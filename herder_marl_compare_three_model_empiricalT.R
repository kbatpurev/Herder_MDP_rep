# ============================================================
# Herder MARL: Action-based PES vs Outcome-based PES vs No PES
# Flat sensitivity grid with persistent worker-local model environments
# ============================================================
#
# Architecture:
#   1. Build one flat model x scenario x replicate task grid.
#   2. Start a persistent PSOCK worker pool once.
#   3. Split the flat grid into contiguous worker chunks.
#   4. Within each worker chunk, instantiate each required model at most once.
#   5. Reuse that model environment across all simulations handled by that worker,
#      resetting the sensitivity parameters before every trajectory.
#
# This avoids the main inefficiency of the earlier flattened version, which
# rebuilt the full model environment for every replicate. Contiguous chunks also
# preserve model locality, so most workers should instantiate only one model.
# Timesteps within simulate_marl() remain sequential.
#
# GitHub workflow:
#   - keep this script + the three Rmd model files under version control;
#   - add outputs/ to .gitignore;
#   - each execution writes its results to one timestamped outputs/sens_.../
#     folder, including metadata and frozen copies of the model Rmds.
# ============================================================
setwd("C:/LocalData/batkhor/OneDrive - University of Helsinki/Papers/Paper 2/mdp_2nd_round/")
library(tidyverse)

if(!requireNamespace("parallelly",quietly=TRUE)){
  stop(
    "This script requires package 'parallelly'. ",
    "Install it once with install.packages('parallelly')."
  )
}

# ============================================================
# 0. USER SETTINGS
# ============================================================

# Exact model filenames expected in the current project directory.
ACTION_MODEL_RMD<-file.path(getwd(),"herder_marl_action_based_PES_csv_transition.Rmd")
OUTCOME_MODEL_RMD<-file.path(getwd(),"herder_marl_outcome_based_PES_csv_transition.Rmd")
NO_PES_MODEL_RMD<-file.path(getwd(),"herder_marl_no_PES_csv_transition.Rmd")

# Full sensitivity run settings.
N_RUNS<-10L
N_STEPS<-300L
SEED_START<-5003L

# The server exposes 128 logical CPUs. Start with 64 persistent workers.
# This flat architecture has hundreds of independent trajectories available,
# so unlike the hierarchical version it can make productive use of >20 workers.
N_WORKERS_TARGET<-6L
N_WORKERS<-as.integer(min(N_WORKERS_TARGET,parallelly::availableCores()))

# Generated files live here. Add "outputs/" to .gitignore.
OUTPUT_ROOT<-file.path(getwd(),"outputs")

# ============================================================
# 1. SCENARIO GRID
# ============================================================
#
# Baseline is run once for each model design.
# Every other scenario changes one model parameter/regime from baseline.
# pes_multiplier scales:
#   Action-based PES  -> pes_payment
#   Outcome-based PES -> state_pes
#   No PES            -> zero PES remains zero (internal control)
# ============================================================

scenario_grid<-tribble(
  ~scenario,~alpha,~gamma,~tau,~memory_window,
  ~drought_prob,~normal_prob,~wet_prob,~pes_multiplier,

  "baseline",
  0.20,0.95,1.00,5L,
  0.30,0.50,0.20,1.00,

  "low_gamma",
  0.20,0.50,1.00,5L,
  0.30,0.50,0.20,1.00,

  "long_memory",
  0.20,0.95,1.00,10L,
  0.30,0.50,0.20,1.00,

  "low_drought",
  0.20,0.95,1.00,5L,
  0.15,0.50,0.35,1.00,

  "high_drought",
  0.20,0.95,1.00,5L,
  0.50,0.35,0.15,1.00,

  "high_PES",
  0.20,0.95,1.00,5L,
  0.30,0.50,0.20,1.50
)

parameter_map<-tribble(
  ~parameter,~variant_scenario,~variant_label,
  "gamma","low_gamma","Low gamma: 0.50",
  "memory_window","long_memory","Long memory: 10",
  "drought_probability_low","low_drought","Low drought: 0.15, normal: 0.50, wet: 0.35",
  "drought_probability","high_drought","High drought: 0.50, normal: 0.35, wet: 0.15",
  "PES_intensity","high_PES","PES multiplier: 1.50"
)

validate_scenario_grid<-function(x){
  required<-c(
    "scenario","alpha","gamma","tau","memory_window",
    "drought_prob","normal_prob","wet_prob","pes_multiplier"
  )
  missing<-setdiff(required,names(x))
  if(length(missing)>0){
    stop("scenario_grid is missing: ",paste(missing,collapse=", "))
  }

  weather_sum<-x$drought_prob+x$normal_prob+x$wet_prob
  if(any(abs(weather_sum-1)>1e-10)){
    stop("Weather probabilities must sum to 1.")
  }
  if(anyDuplicated(x$scenario))stop("Scenario names must be unique.")
  if(any(x$alpha<=0|x$alpha>1))stop("alpha must be in (0,1].")
  if(any(x$gamma<0|x$gamma>=1))stop("gamma must be in [0,1).")
  if(any(x$tau<=0))stop("tau must be > 0.")
  if(any(x$memory_window<1))stop("memory_window must be >= 1.")
  if(any(x$pes_multiplier<0))stop("PES multiplier cannot be negative.")
  invisible(TRUE)
}

validate_scenario_grid(scenario_grid)

# ============================================================
# 2. RUN ID, GIT PROVENANCE, OUTPUT FOLDER, AND FROZEN MODEL FILES
# ============================================================

SENSITIVITY_RUN_ID<-paste0(
  "sens_",
  format(Sys.time(),"%Y%m%d_%H%M%S")
)

SENSITIVITY_RUN_SPEC<-paste0(
  SENSITIVITY_RUN_ID,
  "_runs",N_RUNS,
  "_steps",N_STEPS,
  "_seed",SEED_START
)

RUN_DIR<-file.path(OUTPUT_ROOT,SENSITIVITY_RUN_ID)
PLOT_DIR<-file.path(RUN_DIR,"plots")
RESULT_DIR<-file.path(RUN_DIR,"results")
MODEL_CODE_DIR<-file.path(RUN_DIR,"model_code")

dir.create(PLOT_DIR,recursive=TRUE,showWarnings=FALSE)
dir.create(RESULT_DIR,recursive=TRUE,showWarnings=FALSE)
dir.create(MODEL_CODE_DIR,recursive=TRUE,showWarnings=FALSE)

message("Sensitivity run ID: ",SENSITIVITY_RUN_ID)
message("Output directory: ",RUN_DIR)
message("Available cores reported by parallelly: ",parallelly::availableCores())
message("Workers requested for this run: ",N_WORKERS)

source_model_registry<-tibble(
  model_id=c("action","outcome","none"),
  pes_design=c("Action-based PES","Outcome-based PES","No PES"),
  source_file=c(ACTION_MODEL_RMD,OUTCOME_MODEL_RMD,NO_PES_MODEL_RMD)
)

missing_files<-source_model_registry%>%
  filter(!file.exists(source_file))

if(nrow(missing_files)>0){
  stop(
    "Model file(s) not found:\n",
    paste(
      paste0(missing_files$pes_design,": ",missing_files$source_file),
      collapse="\n"
    ),
    "\n\nWorking directory:\n",getwd(),
    "\n\nFiles currently in working directory:\n",
    paste(list.files(getwd()),collapse="\n")
  )
}

get_git_output<-function(args){
  tryCatch(
    system2("git",args,stdout=TRUE,stderr=FALSE),
    error=function(e)character(0)
  )
}

GIT_COMMIT<-get_git_output(c("rev-parse","HEAD"))
if(length(GIT_COMMIT)==0)GIT_COMMIT<-NA_character_
GIT_COMMIT<-GIT_COMMIT[[1]]

GIT_STATUS<-get_git_output(c("status","--porcelain"))
GIT_DIRTY<-length(GIT_STATUS)>0

# Freeze the exact model code used for this sensitivity run.
walk2(
  source_model_registry$source_file,
  basename(source_model_registry$source_file),
  ~file.copy(.x,file.path(MODEL_CODE_DIR,.y),overwrite=TRUE)
)

model_registry<-source_model_registry%>%
  mutate(
    frozen_file=file.path(MODEL_CODE_DIR,basename(source_file)),
    md5=unname(tools::md5sum(frozen_file)),
    source_modified=as.character(file.info(source_file)$mtime)
  )

write_csv(model_registry,file.path(RUN_DIR,"model_files.csv"))
write_csv(scenario_grid,file.path(RUN_DIR,"scenario_grid.csv"))
write_csv(parameter_map,file.path(RUN_DIR,"parameter_map.csv"))

# ============================================================
# 3. READ AND COMPILE MODEL RMDs ONCE
# ============================================================
#
# The Rmds are parsed once on the main R session. Workers instantiate fresh
# model environments from these frozen expressions. This avoids rereading
# source files hundreds of times and prevents mid-run source edits from
# changing later tasks.
# ============================================================

read_rmd_chunks<-function(path){
  if(!file.exists(path))stop("Model file not found: ",path)
  x<-readLines(path,warn=FALSE)
  starts<-grep("^```\\{r",x)
  chunks<-vector("list",length(starts))

  for(i in seq_along(starts)){
    s<-starts[[i]]
    remaining<-x[(s+1L):length(x)]
    end_candidates<-which(grepl("^```\\s*$",remaining))
    if(length(end_candidates)==0){
      stop("Unclosed R chunk beginning at line ",s," in ",path)
    }
    e<-s+end_candidates[[1]]
    header<-x[[s]]
    inside<-sub("^```\\{r\\s*","",header)
    inside<-sub("\\}\\s*$","",inside)
    label<-trimws(strsplit(inside,",",fixed=TRUE)[[1]][1])
    if(label==""||grepl("=",label,fixed=TRUE))label<-paste0("chunk_",s)
    chunks[[i]]<-list(label=label,code=x[(s+1L):(e-1L)])
  }

  chunks
}

compile_model<-function(path){
  chunks<-read_rmd_chunks(path)
  exclude_pattern<-"^(test-|diagnostic|train-model$|plot-|independent-replicates$)"
  keep<-!str_detect(map_chr(chunks,"label"),exclude_pattern)
  chunks<-chunks[keep]

  map(
    chunks,
    function(ch){
      code<-paste(ch$code,collapse="\n")
      expr<-tryCatch(
        parse(text=code),
        error=function(e){
          stop(
            "Error while parsing chunk '",ch$label,"' in ",path,":\n",
            conditionMessage(e)
          )
        }
      )
      list(label=ch$label,expr=expr)
    }
  )
}

instantiate_model<-function(compiled_model,model_name){
  env<-new.env(parent=globalenv())
  env$model_name<-model_name

  for(ch in compiled_model){
    tryCatch(
      eval(ch$expr,envir=env),
      error=function(e){
        stop(
          "Error while loading ",model_name,
          " from chunk '",ch$label,"':\n",
          conditionMessage(e)
        )
      }
    )
  }

  env
}

compiled_models<-setNames(
  map(model_registry$frozen_file,compile_model),
  model_registry$model_id
)

# Instantiate one baseline copy of each model for validation + metadata.
action_env<-instantiate_model(compiled_models$action,"Action-based PES")
outcome_env<-instantiate_model(compiled_models$outcome,"Outcome-based PES")
no_pes_env<-instantiate_model(compiled_models$none,"No PES")

# ============================================================
# 4. VALIDATE MODEL INTERFACES AND COMPARABILITY
# ============================================================

required_common<-c(
  "simulate_marl",
  "base_transition",
  "climate_transition_source",
  "transition_kernel_version",
  "alpha","gamma","tau","memory_window","ema_alpha",
  "weather_prob","colocation_multiplier",
  "grid_nx","grid_ny",
  "stock_survival","herd_income",
  "designated_cells","initial_actions"
)

validate_model_env<-function(env,pes_design){
  missing<-required_common[
    !vapply(required_common,exists,logical(1),envir=env,inherits=FALSE)
  ]

  if(length(missing)>0){
    stop(pes_design," model is missing: ",paste(missing,collapse=", "))
  }

  if(pes_design=="Action-based PES"&&
     !exists("pes_payment",envir=env,inherits=FALSE)){
    stop("Action-based model must contain pes_payment.")
  }

  if(pes_design=="Outcome-based PES"&&
     !exists("state_pes",envir=env,inherits=FALSE)){
    stop("Outcome-based model must contain state_pes.")
  }

  if(pes_design=="No PES"){
    if(!exists("pes_payment",envir=env,inherits=FALSE)){
      stop("No-PES model must contain pes_payment set to zero.")
    }
    if(any(env$pes_payment!=0)){
      stop("No-PES model has non-zero PES payments.")
    }
  }

  invisible(TRUE)
}

validate_model_env(action_env,"Action-based PES")
validate_model_env(outcome_env,"Outcome-based PES")
validate_model_env(no_pes_env,"No PES")

message("Transition kernel version: ",action_env$transition_kernel_version)

same_transition_kernel<-function(env1,env2){
  identical(env1$transition_kernel_version,env2$transition_kernel_version)&&
    identical(env1$climate_transition_source,env2$climate_transition_source)&&
    identical(env1$base_transition,env2$base_transition)
}

model_structure_check<-tibble(
  feature=c(
    "Grid cells",
    "Agents",
    "Designated-cell rows",
    "Unique designated cells",
    "Emergency cost scale present",
    "Empirical transition kernel identical",
    "Stock survival identical",
    "Herd income identical"
  ),
  action_based=c(
    action_env$grid_nx*action_env$grid_ny,
    length(action_env$initial_actions),
    nrow(action_env$designated_cells),
    n_distinct(action_env$designated_cells$cell_id),
    exists("emergency_cost_scale",envir=action_env,inherits=FALSE),
    TRUE,
    TRUE,
    TRUE
  ),
  outcome_based=c(
    outcome_env$grid_nx*outcome_env$grid_ny,
    length(outcome_env$initial_actions),
    nrow(outcome_env$designated_cells),
    n_distinct(outcome_env$designated_cells$cell_id),
    exists("emergency_cost_scale",envir=outcome_env,inherits=FALSE),
    same_transition_kernel(action_env,outcome_env),
    identical(action_env$stock_survival,outcome_env$stock_survival),
    identical(action_env$herd_income,outcome_env$herd_income)
  ),
  no_pes=c(
    no_pes_env$grid_nx*no_pes_env$grid_ny,
    length(no_pes_env$initial_actions),
    nrow(no_pes_env$designated_cells),
    n_distinct(no_pes_env$designated_cells$cell_id),
    exists("emergency_cost_scale",envir=no_pes_env,inherits=FALSE),
    same_transition_kernel(action_env,no_pes_env),
    identical(action_env$stock_survival,no_pes_env$stock_survival),
    identical(action_env$herd_income,no_pes_env$herd_income)
  )
)

print(model_structure_check)
write_csv(model_structure_check,file.path(RUN_DIR,"model_structure_check.csv"))

structure_same<-
  nrow(action_env$designated_cells)==nrow(outcome_env$designated_cells)&&
  nrow(action_env$designated_cells)==nrow(no_pes_env$designated_cells)&&
  n_distinct(action_env$designated_cells$cell_id)==
    n_distinct(outcome_env$designated_cells$cell_id)&&
  n_distinct(action_env$designated_cells$cell_id)==
    n_distinct(no_pes_env$designated_cells$cell_id)&&
  same_transition_kernel(action_env,outcome_env)&&
  same_transition_kernel(action_env,no_pes_env)&&
  identical(action_env$stock_survival,outcome_env$stock_survival)&&
  identical(action_env$stock_survival,no_pes_env$stock_survival)&&
  identical(action_env$herd_income,outcome_env$herd_income)&&
  identical(action_env$herd_income,no_pes_env$herd_income)

if(!structure_same){
  warning(
    "The three model files differ in shared ecological/economic structure. ",
    "Model-design comparisons may be confounded. Inspect model_structure_check.csv."
  )
}

# ============================================================
# 5. BASELINE SNAPSHOTS + MODEL PARAMETER METADATA
# ============================================================

snapshot_model<-function(env,pes_design){
  out<-list(
    alpha=env$alpha,
    gamma=env$gamma,
    tau=env$tau,
    memory_window=env$memory_window,
    ema_alpha=env$ema_alpha,
    weather_prob=env$weather_prob,
    colocation_multiplier=env$colocation_multiplier
  )

  if(exists("emergency_cost_scale",envir=env,inherits=FALSE)){
    out$emergency_cost_scale<-env$emergency_cost_scale
  }

  if(pes_design=="Action-based PES"){
    out$pes_vector<-env$pes_payment
  }else if(pes_design=="Outcome-based PES"){
    out$pes_vector<-env$state_pes
  }else{
    out$pes_vector<-env$pes_payment
  }

  out
}

action_baseline<-snapshot_model(action_env,"Action-based PES")
outcome_baseline<-snapshot_model(outcome_env,"Outcome-based PES")
no_pes_baseline<-snapshot_model(no_pes_env,"No PES")

baseline_registry<-list(
  action=action_baseline,
  outcome=outcome_baseline,
  none=no_pes_baseline
)

format_parameter<-function(x){
  if(is.null(x))return(NA_character_)
  if(length(x)==1L)return(as.character(x))
  if(!is.null(names(x))){
    return(paste(paste0(names(x),"=",x),collapse="; "))
  }
  paste(x,collapse="; ")
}

metadata_parameters<-c(
  "alpha","gamma","tau","memory_window","ema_alpha",
  "weather_prob","colocation_multiplier","transition_kernel_version",
  "travel_cost_weight","social_cost_weight","emergency_cost_scale",
  "pes_payment","state_pes","herd_income"
)

extract_model_parameters<-function(env,model_name){
  map_dfr(
    metadata_parameters,
    function(parameter){
      if(exists(parameter,envir=env,inherits=FALSE)){
        value<-get(parameter,envir=env,inherits=FALSE)
        tibble(
          pes_design=model_name,
          parameter=parameter,
          value=format_parameter(value)
        )
      }else{
        tibble(
          pes_design=model_name,
          parameter=parameter,
          value=NA_character_
        )
      }
    }
  )
}

model_parameters<-bind_rows(
  extract_model_parameters(action_env,"Action-based PES"),
  extract_model_parameters(outcome_env,"Outcome-based PES"),
  extract_model_parameters(no_pes_env,"No PES")
)

write_csv(model_parameters,file.path(RUN_DIR,"model_parameters.csv"))

# ============================================================
# 6. APPLY ONE SCENARIO TO A FRESH MODEL ENVIRONMENT
# ============================================================

apply_scenario<-function(env,baseline,pes_design,scenario_row){
  if(nrow(scenario_row)!=1L){
    stop("scenario_row must have exactly one row.")
  }

  env$alpha<-scenario_row$alpha[[1]]
  env$gamma<-scenario_row$gamma[[1]]
  env$tau<-scenario_row$tau[[1]]

  env$memory_window<-as.integer(scenario_row$memory_window[[1]])
  env$ema_alpha<-2/(env$memory_window+1)

  env$weather_prob<-c(
    drought=scenario_row$drought_prob[[1]],
    normal=scenario_row$normal_prob[[1]],
    wet=scenario_row$wet_prob[[1]]
  )

  if(pes_design=="Action-based PES"){
    env$pes_payment<-baseline$pes_vector*scenario_row$pes_multiplier[[1]]
  }else if(pes_design=="Outcome-based PES"){
    env$state_pes<-baseline$pes_vector*scenario_row$pes_multiplier[[1]]
  }else{
    # No-PES remains exactly zero under every scenario.
    env$pes_payment<-baseline$pes_vector
  }

  invisible(NULL)
}

# ============================================================
# 7. BUILD FLAT MODEL x SCENARIO x REPLICATE TASK GRID
# ============================================================
#
# Every row below is one independent MARL trajectory. The No-PES x high_PES
# combination is omitted because PES intensity is not defined for No PES.
#
# The grid is explicitly ordered model -> scenario -> run. We later split these
# rows into contiguous worker chunks, which means most workers handle only one
# model and therefore instantiate that model only once.
# ============================================================

scenario_columns<-names(scenario_grid)
model_order<-c("action","outcome","none")

design_grid<-tidyr::expand_grid(
  model_id=model_order,
  scenario=scenario_grid$scenario
)%>%
  filter(!(model_id=="none"&scenario=="high_PES"))%>%
  left_join(
    model_registry%>%select(model_id,pes_design,md5),
    by="model_id"
  )%>%
  left_join(scenario_grid,by="scenario")%>%
  mutate(
    model_index=match(model_id,model_order),
    scenario_index=match(scenario,scenario_grid$scenario)
  )%>%
  arrange(model_index,scenario_index)%>%
  select(-model_index,-scenario_index)

task_grid<-tidyr::expand_grid(
  design_row=seq_len(nrow(design_grid)),
  run=seq_len(N_RUNS)
)%>%
  left_join(
    design_grid%>%mutate(design_row=row_number()),
    by="design_row"
  )%>%
  arrange(design_row,run)%>%
  mutate(
    task_id=row_number(),
    steps=N_STEPS,
    seed=SEED_START+run-1L,
    sensitivity_run_id=SENSITIVITY_RUN_ID,
    sensitivity_run_spec=SENSITIVITY_RUN_SPEC
  )%>%
  select(-design_row)

EXPECTED_SIMULATIONS<-nrow(design_grid)*N_RUNS
stopifnot(nrow(task_grid)==EXPECTED_SIMULATIONS)

# Do not create more workers than there are trajectories.
N_WORKERS<-min(N_WORKERS,nrow(task_grid))

# Contiguous, nearly equal-sized chunks. parallel::splitIndices() preserves
# contiguous task indices, so model locality is retained.
task_chunks<-parallel::splitIndices(
  nx=nrow(task_grid),
  ncl=N_WORKERS
)

chunk_plan<-map_dfr(
  seq_along(task_chunks),
  function(worker_slot){
    idx<-task_chunks[[worker_slot]]
    tibble(
      worker_slot=worker_slot,
      n_tasks=length(idx),
      first_task=min(idx),
      last_task=max(idx),
      models=paste(unique(task_grid$model_id[idx]),collapse=";"),
      scenarios=paste(unique(task_grid$scenario[idx]),collapse=";"),
      model_instantiations=n_distinct(task_grid$model_id[idx])
    )
  }
)

EXPECTED_MODEL_INSTANTIATIONS<-sum(chunk_plan$model_instantiations)

message(
  "Flat task grid: ",nrow(task_grid)," independent simulations across ",
  nrow(design_grid)," model/scenario combinations."
)
message(
  "Persistent worker plan: ",N_WORKERS," workers; ",
  min(chunk_plan$n_tasks),"-",max(chunk_plan$n_tasks)," simulations per worker."
)
message(
  "Expected model instantiations across all workers: ",
  EXPECTED_MODEL_INSTANTIATIONS,
  " (instead of ",nrow(task_grid)," in the earlier per-replicate flat version)."
)

write_csv(task_grid,file.path(RUN_DIR,"task_grid.csv"))
write_csv(chunk_plan,file.path(RUN_DIR,"worker_chunk_plan.csv"))

# ============================================================
# 8. RUN ONE PERSISTENT WORKER CHUNK
# ============================================================
#
# A worker receives several contiguous simulation tasks. Its local cache starts
# empty. The first time a model is needed, that model is instantiated from the
# already-compiled frozen Rmd expressions. The same environment is then reused
# for all later tasks for that model on that worker.
#
# This is safe for the current model implementation because simulate_marl()
# creates agents, land grid, Q tables, weather, and logs locally for each call.
# apply_scenario() resets every sensitivity-varying model parameter before each
# trajectory, including the PES vector from the immutable baseline snapshot.
# ============================================================

run_task_chunk<-function(task_indices){
  model_cache<-new.env(parent=emptyenv())

  one_simulation<-function(i){
    task_row<-task_grid[i,,drop=FALSE]
    model_id<-task_row$model_id[[1]]
    pes_design<-task_row$pes_design[[1]]
    scenario_name<-task_row$scenario[[1]]
    task_id<-task_row$task_id[[1]]
    run_id<-task_row$run[[1]]
    simulation_seed<-task_row$seed[[1]]
    simulation_steps<-task_row$steps[[1]]
    run_stamp<-task_row$sensitivity_run_id[[1]]
    run_spec<-task_row$sensitivity_run_spec[[1]]

    if(!exists(model_id,envir=model_cache,inherits=FALSE)){
      assign(
        model_id,
        instantiate_model(
          compiled_models[[model_id]],
          pes_design
        ),
        envir=model_cache
      )
    }

    env<-get(model_id,envir=model_cache,inherits=FALSE)

    scenario_row<-task_row[,scenario_columns,drop=FALSE]

    apply_scenario(
      env=env,
      baseline=baseline_registry[[model_id]],
      pes_design=pes_design,
      scenario_row=scenario_row
    )

    res<-env$simulate_marl(
      steps=simulation_steps,
      seed=simulation_seed
    )

    landscape<-res$landscape_counts%>%
      mutate(
        sensitivity_run_id=run_stamp,
        sensitivity_run_spec=run_spec,
        task_id=task_id,
        seed=simulation_seed,
        pes_design=pes_design,
        scenario=scenario_name,
        run=run_id,
        prop_degraded=degraded/(env$grid_nx*env$grid_ny)
      )

    cells<-res$cell_ledger%>%
      mutate(
        sensitivity_run_id=run_stamp,
        sensitivity_run_spec=run_spec,
        task_id=task_id,
        seed=simulation_seed,
        pes_design=pes_design,
        scenario=scenario_name,
        run=run_id
      )

    actions<-res$action_log%>%
      mutate(
        sensitivity_run_id=run_stamp,
        sensitivity_run_spec=run_spec,
        task_id=task_id,
        seed=simulation_seed,
        pes_design=pes_design,
        scenario=scenario_name,
        run=run_id
      )

    if(pes_design=="Outcome-based PES"&&
       "state_pes_reward"%in%names(actions)){
      actions$pes_signal<-actions$state_pes_reward
    }else if("pes_reward"%in%names(actions)){
      actions$pes_signal<-actions$pes_reward
    }else{
      actions$pes_signal<-NA_real_
    }

    if(!is.null(res$q_log)){
      q<-res$q_log%>%
        mutate(
          sensitivity_run_id=run_stamp,
          sensitivity_run_spec=run_spec,
          task_id=task_id,
          seed=simulation_seed,
          pes_design=pes_design,
          scenario=scenario_name,
          run=run_id
        )
    }else{
      q<-tibble()
    }

    list(
      landscape=landscape,
      cells=cells,
      actions=actions,
      q=q
    )
  }

  simulation_results<-lapply(task_indices,one_simulation)

  list(
    worker_pid=Sys.getpid(),
    n_tasks=length(task_indices),
    models_instantiated=ls(model_cache,all.names=TRUE),
    landscape=bind_rows(lapply(simulation_results,`[[`,"landscape")),
    cells=bind_rows(lapply(simulation_results,`[[`,"cells")),
    actions=bind_rows(lapply(simulation_results,`[[`,"actions")),
    q=bind_rows(lapply(simulation_results,`[[`,"q"))
  )
}

# ============================================================
# 9. METADATA BEFORE THE EXPENSIVE RUN
# ============================================================

RUN_START<-Sys.time()

metadata_start<-c(
  paste0("Sensitivity run ID: ",SENSITIVITY_RUN_ID),
  paste0("Sensitivity run spec: ",SENSITIVITY_RUN_SPEC),
  paste0("Start time: ",RUN_START),
  paste0("Working directory: ",getwd()),
  paste0("Git commit: ",GIT_COMMIT),
  paste0("Uncommitted Git changes at start: ",ifelse(GIT_DIRTY,"Yes","No")),
  "",
  "SIMULATION SETTINGS",
  paste0("Models: 3"),
  paste0("Active scenarios: ",nrow(scenario_grid)),
  paste0("Model/scenario combinations: ",nrow(design_grid)),
  paste0("Independent runs per model/scenario: ",N_RUNS),
  paste0("Iterations per run: ",N_STEPS),
  paste0("Independent simulations: ",EXPECTED_SIMULATIONS),
  paste0("Persistent workers: ",N_WORKERS),
  paste0("Worker chunks: ",length(task_chunks)),
  paste0("Expected model instantiations across workers: ",EXPECTED_MODEL_INSTANTIATIONS),
  paste0("Seed start: ",SEED_START),
  "",
  "SOFTWARE",
  paste0("R version: ",R.version.string),
  paste0("Platform: ",R.version$platform)
)

writeLines(metadata_start,file.path(RUN_DIR,"metadata.txt"))
writeLines(capture.output(sessionInfo()),file.path(RUN_DIR,"session_info.txt"))

# ============================================================
# 10. START ONE PERSISTENT PSOCK POOL AND RUN ALL CHUNKS
# ============================================================
#
# Unlike future_map(), the cluster is created once and retained for the whole
# simulation phase. Compiled model expressions, baseline snapshots, and the task
# grid are exported once. Because task_chunks has exactly N_WORKERS elements,
# clusterApply() assigns exactly one contiguous chunk to each worker.
# ============================================================

CLUSTER_START<-Sys.time()
message("Starting ",N_WORKERS," persistent PSOCK workers.")

cl<-parallel::makePSOCKcluster(
  N_WORKERS,
  outfile=""
)

worker_pids<-unlist(parallel::clusterCall(cl,Sys.getpid))
stopifnot(length(unique(worker_pids))==N_WORKERS)
message("Persistent worker processes confirmed: ",length(unique(worker_pids)))

# Load packages once per worker. The model setup chunks also call library(),
# but those later calls are then effectively no-ops.
parallel::clusterEvalQ(
  cl,
  {
    suppressPackageStartupMessages(library(tidyverse))
    NULL
  }
)

parallel::clusterExport(
  cl,
  varlist=c(
    "compiled_models",
    "baseline_registry",
    "task_grid",
    "scenario_columns",
    "instantiate_model",
    "apply_scenario"
  ),
  envir=environment()
)

CLUSTER_READY<-Sys.time()
CLUSTER_INIT_MIN<-as.numeric(
  difftime(CLUSTER_READY,CLUSTER_START,units="mins")
)
message(
  "Worker pool initialized in ",
  round(CLUSTER_INIT_MIN,2),
  " minutes."
)

SIMULATION_START<-Sys.time()

chunk_results<-tryCatch(
  parallel::clusterApply(
    cl,
    task_chunks,
    run_task_chunk
  ),
  finally={
    parallel::stopCluster(cl)
  }
)

SIMULATION_END<-Sys.time()
SIMULATION_DURATION_MIN<-as.numeric(
  difftime(SIMULATION_END,SIMULATION_START,units="mins")
)
RUN_END<-SIMULATION_END
RUN_DURATION_MIN<-as.numeric(difftime(RUN_END,RUN_START,units="mins"))

message(
  "Parallel simulation phase finished in ",
  round(SIMULATION_DURATION_MIN,2),
  " minutes."
)

worker_summary<-map_dfr(
  seq_along(chunk_results),
  function(i){
    tibble(
      worker_slot=i,
      worker_pid=chunk_results[[i]]$worker_pid,
      n_tasks=chunk_results[[i]]$n_tasks,
      models_instantiated=paste(
        chunk_results[[i]]$models_instantiated,
        collapse=";"
      ),
      n_models_instantiated=length(
        chunk_results[[i]]$models_instantiated
      )
    )
  }
)

write_csv(worker_summary,file.path(RUN_DIR,"worker_summary.csv"))

# ============================================================
# 11. REASSEMBLE comparison_results + VALIDATE TASK COMPLETENESS
# ============================================================

comparison_results<-list(
  landscape=map_dfr(chunk_results,"landscape"),
  cells=map_dfr(chunk_results,"cells"),
  actions=map_dfr(chunk_results,"actions"),
  q=map_dfr(chunk_results,"q"),
  scenarios=scenario_grid,
  parameter_map=parameter_map,
  structure_check=model_structure_check,
  task_grid=task_grid,
  sensitivity_run_id=SENSITIVITY_RUN_ID,
  sensitivity_run_spec=SENSITIVITY_RUN_SPEC
)

run_validation<-comparison_results$landscape%>%
  summarise(
    tasks=n_distinct(task_id),
    models=n_distinct(pes_design),
    scenarios=n_distinct(scenario),
    simulations=n_distinct(interaction(pes_design,scenario,run,drop=TRUE)),
    runs=n_distinct(run),
    run_ids=n_distinct(sensitivity_run_id),
    max_step=max(step,na.rm=TRUE)
  )

print(run_validation)

stopifnot(
  run_validation$tasks[[1]]==nrow(task_grid),
  run_validation$models[[1]]==3L,
  run_validation$scenarios[[1]]==nrow(scenario_grid),
  run_validation$simulations[[1]]==EXPECTED_SIMULATIONS,
  run_validation$runs[[1]]==N_RUNS,
  run_validation$run_ids[[1]]==1L,
  run_validation$max_step[[1]]==N_STEPS
)

returned_task_ids<-sort(unique(comparison_results$landscape$task_id))
planned_task_ids<-sort(task_grid$task_id)
stopifnot(identical(returned_task_ids,planned_task_ids))

# No-PES intentionally has no high_PES sensitivity run.
stopifnot(
  !any(
    comparison_results$landscape$pes_design=="No PES"&
      comparison_results$landscape$scenario=="high_PES"
  )
)

write_csv(run_validation,file.path(RUN_DIR,"run_validation.csv"))
saveRDS(comparison_results,file.path(RESULT_DIR,"comparison_results.rds"))

metadata_end<-c(
  "",
  "RUN COMPLETED",
  paste0("End time: ",RUN_END),
  paste0("Elapsed minutes: ",round(RUN_DURATION_MIN,2)),
  paste0("Returned simulation tasks: ",run_validation$tasks[[1]]),
  paste0("Returned independent simulations: ",run_validation$simulations[[1]]),
  paste0("Cluster initialization minutes: ",round(CLUSTER_INIT_MIN,2)),
  paste0("Parallel simulation minutes: ",round(SIMULATION_DURATION_MIN,2))
)

write(metadata_end,file.path(RUN_DIR,"metadata.txt"),append=TRUE)

# ============================================================
# 12. PARAMETER LOOKUP + PREPARATION
# ============================================================

get_parameter_scenarios<-function(parameter_name){
  row<-parameter_map%>%filter(parameter==parameter_name)
  if(nrow(row)!=1){
    stop(
      "parameter_name must be one of: ",
      paste(parameter_map$parameter,collapse=", ")
    )
  }
  list(
    scenarios=c("baseline",row$variant_scenario[[1]]),
    variant=row$variant_scenario[[1]],
    variant_label=row$variant_label[[1]]
  )
}

prepare_parameter_data<-function(data,parameter_name){
  info<-get_parameter_scenarios(parameter_name)
  data%>%
    filter(scenario%in%info$scenarios)%>%
    mutate(
      setting=if_else(scenario=="baseline","Baseline",info$variant_label),
      setting=factor(setting,levels=c("Baseline",info$variant_label)),
      pes_design=factor(
        pes_design,
        levels=c("Action-based PES","Outcome-based PES","No PES")
      )
    )
}

# ============================================================
# 13. LANDSCAPE DEGRADATION SENSITIVITY PLOTS
# ============================================================

make_parameter_degradation_plot<-function(results,parameter_name){
  plot_data<-prepare_parameter_data(results$landscape,parameter_name)

  summary_data<-plot_data%>%
    group_by(pes_design,setting,step)%>%
    summarise(
      mean_prop_degraded=mean(prop_degraded),
      lower=quantile(prop_degraded,0.025),
      upper=quantile(prop_degraded,0.975),
      .groups="drop"
    )

  ggplot()+
    geom_line(
      data=plot_data,
      aes(
        x=step,
        y=prop_degraded,
        colour=setting,
        group=interaction(setting,run)
      ),
      alpha=0.15,
      linewidth=0.30
    )+
    geom_ribbon(
      data=summary_data,
      aes(x=step,ymin=lower,ymax=upper,fill=setting),
      alpha=0.10,
      colour=NA
    )+
    geom_line(
      data=summary_data,
      aes(x=step,y=mean_prop_degraded,colour=setting),
      linewidth=1
    )+
    facet_wrap(~pes_design,nrow=1)+
    scale_y_continuous(
      limits=c(0,1),
      labels=scales::percent_format(accuracy=1)
    )+
    labs(
      title=paste("Sensitivity to",parameter_name),
      subtitle=paste0("Sensitivity run: ",results$sensitivity_run_id),
      x="Iteration",
      y="Proportion of degraded cells",
      colour="Setting",
      fill="Setting"
    )+
    theme_minimal()+
    theme(legend.position="bottom")
}

parameter_degradation_plots<-setNames(
  map(
    parameter_map$parameter,
    ~make_parameter_degradation_plot(comparison_results,.x)
  ),
  parameter_map$parameter
)

safe_plot_name<-function(x){
  str_replace_all(x,"[^A-Za-z0-9_-]+","_")
}

walk2(
  parameter_degradation_plots,
  names(parameter_degradation_plots),
  function(plot_object,plot_name){
    ggsave(
      filename=file.path(
        PLOT_DIR,
        paste0("degradation_",safe_plot_name(plot_name),".png")
      ),
      plot=plot_object,
      width=11,
      height=6,
      dpi=300
    )
  }
)

# ============================================================
# 14. OPTIONAL DIAGNOSTIC PLOT FUNCTIONS
# ============================================================

plot_cell_degradation_comparison<-function(
    results,
    parameter_name,
    run_id=1L
){
  plot_data<-prepare_parameter_data(results$cells,parameter_name)%>%
    filter(run==run_id)%>%
    mutate(
      cell_id=factor(
        cell_id,
        levels=rev(sort(unique(cell_id)))
      )
    )

  if(nrow(plot_data)==0)stop("No cell data found for run ",run_id,".")

  ggplot(plot_data,aes(x=step,y=cell_id,fill=ema_after))+
    geom_tile()+
    scale_fill_viridis_c(
      limits=c(0,3),
      breaks=c(0,0.5,1.5,2.5,3)
    )+
    facet_grid(setting~pes_design)+
    labs(
      title=paste("Cell-level EMA sensitivity to",parameter_name),
      subtitle=paste0("Run ",run_id," | sensitivity run: ",results$sensitivity_run_id),
      x="Iteration",
      y="Cell ID",
      fill="EMA degradation"
    )+
    theme_minimal()
}

plot_action_mix_comparison<-function(results,parameter_name){
  action_order<-c("<200","200-400","400-600",">600")

  plot_data<-prepare_parameter_data(results$actions,parameter_name)%>%
    mutate(action=factor(action,levels=action_order))%>%
    count(pes_design,setting,step,action,name="n")%>%
    group_by(pes_design,setting,step)%>%
    mutate(prop=n/sum(n))%>%
    ungroup()

  ggplot(plot_data,aes(step,prop,fill=action))+
    geom_area()+
    facet_grid(setting~pes_design)+
    scale_y_continuous(labels=scales::percent_format(accuracy=1))+
    labs(
      title=paste("Action mix sensitivity to",parameter_name),
      subtitle=paste0("Sensitivity run: ",results$sensitivity_run_id),
      x="Iteration",
      y="Proportion of agent actions",
      fill="Stocking action"
    )+
    theme_minimal()
}

# ============================================================
# 15. FINAL-STEP SUMMARY
# ============================================================

summarise_final_landscape<-function(results){
  last_step<-max(results$landscape$step)
  results$landscape%>%
    filter(step==last_step)%>%
    group_by(
      sensitivity_run_id,
      pes_design,
      scenario
    )%>%
    summarise(
      mean_prop_degraded=mean(prop_degraded),
      sd_prop_degraded=sd(prop_degraded),
      median_prop_degraded=median(prop_degraded),
      mean_ema=mean(mean_ema),
      .groups="drop"
    )%>%
    arrange(scenario,pes_design)
}

final_landscape_summary<-summarise_final_landscape(comparison_results)
print(final_landscape_summary)
write_csv(
  final_landscape_summary,
  file.path(RESULT_DIR,"final_landscape_summary.csv")
)

# ============================================================
# 16. CUMULATIVE LANDSCAPE DEGRADATION (NORMALIZED AUC)
# ============================================================

landscape_auc<-comparison_results$landscape%>%
  arrange(
    sensitivity_run_id,
    pes_design,
    scenario,
    run,
    step
  )%>%
  group_by(
    sensitivity_run_id,
    sensitivity_run_spec,
    pes_design,
    scenario,
    run
  )%>%
  summarise(
    auc=sum(
      diff(step)*
        (head(prop_degraded,-1)+tail(prop_degraded,-1))/2
    ),
    time_span=max(step)-min(step),
    normalized_auc=100*auc/time_span,
    .groups="drop"
  )

stopifnot(
  n_distinct(landscape_auc$sensitivity_run_id)==1L,
  unique(landscape_auc$sensitivity_run_id)==comparison_results$sensitivity_run_id
)

scenario_order<-scenario_grid$scenario
pes_order<-c("Action-based PES","Outcome-based PES","No PES")
pes_labels<-c("Action","Outcome","None")

landscape_auc<-landscape_auc%>%
  mutate(
    scenario=factor(scenario,levels=scenario_order),
    pes_design=factor(
      pes_design,
      levels=pes_order,
      labels=pes_labels
    )
  )

p_degradation<-ggplot(
  landscape_auc,
  aes(x=scenario,y=normalized_auc,fill=pes_design)
)+
  geom_boxplot(
    position=position_dodge(width=0.8),
    width=0.7,
    outlier.shape=NA
  )+
  geom_point(
    aes(group=pes_design),
    position=position_jitterdodge(
      dodge.width=0.8,
      jitter.width=0.08
    ),
    alpha=0.35,
    size=1.2
  )+
  labs(
    title="Cumulative landscape degradation",
    subtitle=paste0(
      "Run: ",SENSITIVITY_RUN_ID,
      " | ",N_RUNS," replicates x ",N_STEPS," steps"
    ),
    x="Sensitivity scenario",
    y="Cumulative landscape degradation (%)",
    fill="PES design"
  )+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45,hjust=1))

print(p_degradation)

ggsave(
  file.path(PLOT_DIR,"cumulative_landscape_degradation.png"),
  p_degradation,
  width=11,
  height=6,
  dpi=300
)

write_csv(landscape_auc,file.path(RESULT_DIR,"landscape_auc.csv"))

# ============================================================
# 17. CUMULATIVE REWARD, LIVESTOCK INCOME, AND PES SIGNAL
# ============================================================

cumulative_economic_results<-comparison_results$actions%>%
  group_by(
    sensitivity_run_id,
    sensitivity_run_spec,
    pes_design,
    scenario,
    run
  )%>%
  summarise(
    cumulative_reward=sum(reward,na.rm=TRUE),
    cumulative_livestock_income=sum(livestock_income,na.rm=TRUE),
    cumulative_pes_signal=sum(pes_signal,na.rm=TRUE),
    .groups="drop"
  )%>%
  mutate(
    scenario=factor(scenario,levels=scenario_order),
    pes_design=factor(
      pes_design,
      levels=pes_order,
      labels=pes_labels
    )
  )

stopifnot(
  n_distinct(cumulative_economic_results$sensitivity_run_id)==1L,
  unique(cumulative_economic_results$sensitivity_run_id)==comparison_results$sensitivity_run_id
)

p_reward<-ggplot(
  cumulative_economic_results,
  aes(x=scenario,y=cumulative_reward,fill=pes_design)
)+
  geom_boxplot(
    position=position_dodge(width=0.8),
    width=0.7,
    outlier.shape=NA
  )+
  geom_point(
    aes(group=pes_design),
    position=position_jitterdodge(
      dodge.width=0.8,
      jitter.width=0.08
    ),
    alpha=0.35,
    size=1.2
  )+
  labs(
    title="Cumulative reward",
    subtitle=paste0("Sensitivity run: ",SENSITIVITY_RUN_ID),
    x="Sensitivity scenario",
    y="Cumulative reward",
    fill="PES design"
  )+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45,hjust=1))

p_income<-ggplot(
  cumulative_economic_results,
  aes(x=scenario,y=cumulative_livestock_income,fill=pes_design)
)+
  geom_boxplot(
    position=position_dodge(width=0.8),
    width=0.7,
    outlier.shape=NA
  )+
  geom_point(
    aes(group=pes_design),
    position=position_jitterdodge(
      dodge.width=0.8,
      jitter.width=0.08
    ),
    alpha=0.35,
    size=1.2
  )+
  labs(
    title="Cumulative livestock income",
    subtitle=paste0("Sensitivity run: ",SENSITIVITY_RUN_ID),
    x="Sensitivity scenario",
    y="Cumulative livestock income",
    fill="PES design"
  )+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45,hjust=1))

print(p_reward)
print(p_income)

ggsave(
  file.path(PLOT_DIR,"cumulative_reward.png"),
  p_reward,
  width=11,
  height=6,
  dpi=300
)

ggsave(
  file.path(PLOT_DIR,"cumulative_livestock_income.png"),
  p_income,
  width=11,
  height=6,
  dpi=300
)

write_csv(
  cumulative_economic_results,
  file.path(RESULT_DIR,"cumulative_reward_income_pes.csv")
)

# ============================================================
# 18. FINISH
# ============================================================

message("Sensitivity analysis complete.")
message("Run ID: ",SENSITIVITY_RUN_ID)
message("Results: ",RESULT_DIR)
message("Plots: ",PLOT_DIR)
message("Full result object: ",file.path(RESULT_DIR,"comparison_results.rds"))
message("Cluster initialization: ",round(CLUSTER_INIT_MIN,2)," min")
message("Parallel simulations: ",round(SIMULATION_DURATION_MIN,2)," min")
message("Run through simulation completion: ",round(RUN_DURATION_MIN,2)," min")

# Useful in-memory objects:
# comparison_results
# comparison_results$landscape
# comparison_results$cells
# comparison_results$actions
# comparison_results$q
# comparison_results$task_grid
parameter_degradation_plots
# final_landscape_summary
# landscape_auc
# cumulative_economic_results
# p_degradation
# p_reward
# p_income
