# ============================================================
# Herder MARL: Action-based PES vs Outcome-based PES vs No PES sensitivity tests
# ============================================================
#
# Purpose
# -------
# Run the SAME parameter scenarios through:
#   1. Action-based PES model
#   2. Outcome-based PES model
#   3. No-PES model
#
# The three model files are loaded into separate R environments so
# functions/parameters with identical names do not overwrite each other.
#
# For each shared parameter, the main plot places:
#   Action-based PES | Outcome-based PES | No PES
# side-by-side in facets, with baseline and altered parameter settings
# shown within each panel.
#
# PES intensity is not defined for the No-PES model, so the PES-intensity
# sensitivity plot contains only the two PES models.
#
# IMPORTANT
# ---------
# This script compares the three model files AS CURRENTLY DEFINED.
# The structural diagnostic below checks that the common movement/ecological
# assumptions remain comparable across models.
#
# Required files in the working directory:
#   herder_marl_action_based_PES_csv_transition.Rmd
#   herder_marl_outcome_based_PES_csv_transition.Rmd
#   herder_marl_no_PES_csv_transition.Rmd
#
# ============================================================

library(tidyverse)

#Exact model filenames expected in the CURRENT working directory
ACTION_MODEL_RMD<-file.path(getwd(),"herder_marl_action_based_PES_csv_transition.Rmd")
OUTCOME_MODEL_RMD<-file.path(getwd(),"herder_marl_outcome_based_PES_csv_transition.Rmd")
NO_PES_MODEL_RMD<-file.path(getwd(),"herder_marl_no_PES_csv_transition.Rmd")

model_files<-c(
  "Action-based PES"=ACTION_MODEL_RMD,
  "Outcome-based PES"=OUTCOME_MODEL_RMD,
  "No PES"=NO_PES_MODEL_RMD
)

missing_files<-model_files[!file.exists(model_files)]
if(length(missing_files)>0){
  stop(
    "Model file(s) not found:\n",
    paste(names(missing_files),missing_files,sep=": ",collapse="\n"),
    "\n\nWorking directory:\n",getwd(),
    "\n\nFiles currently in working directory:\n",
    paste(list.files(getwd()),collapse="\n")
  )
}

message("Working directory: ",getwd())
walk2(names(model_files),model_files,~message(.x,": ",basename(.y)))

# ============================================================
# 1. LOAD EACH RMD INTO AN ISOLATED MODEL ENVIRONMENT
# ============================================================

read_rmd_chunks<-function(path){
  if(!file.exists(path))stop("Model file not found: ",path)
  x<-readLines(path,warn=FALSE)
  starts<-grep("^```\\{r",x)
  chunks<-vector("list",length(starts))
  for(i in seq_along(starts)){
    s<-starts[[i]]
    end_rel<-which(grepl("^```\\s*$",x[(s+1L):length(x)]))[[1]]
    e<-s+end_rel
    header<-x[[s]]
    inside<-sub("^```\\{r\\s*","",header)
    inside<-sub("\\}\\s*$","",inside)
    label<-trimws(strsplit(inside,",",fixed=TRUE)[[1]][1])
    if(label==""||grepl("=",label,fixed=TRUE))label<-paste0("chunk_",s)
    chunks[[i]]<-list(label=label,code=x[(s+1L):(e-1L)])
  }
  chunks
}

load_model_env<-function(path,model_name){
  chunks<-read_rmd_chunks(path)
  exclude_pattern<-"^(test-|diagnostic|train-model$|plot-|independent-replicates$)"
  keep<-!str_detect(map_chr(chunks,"label"),exclude_pattern)
  chunks<-chunks[keep]
  env<-new.env(parent=globalenv())
  env$model_name<-model_name
  for(ch in chunks){
    code<-paste(ch$code,collapse="\n")
    tryCatch(
      eval(parse(text=code),envir=env),
      error=function(e){
        stop(
          "Error while loading ",model_name,
          " from chunk '",ch$label,"':\n",conditionMessage(e)
        )
      }
    )
  }
  env
}

action_env<-load_model_env(ACTION_MODEL_RMD,"Action-based PES")
outcome_env<-load_model_env(OUTCOME_MODEL_RMD,"Outcome-based PES")
no_pes_env<-load_model_env(NO_PES_MODEL_RMD,"No PES")

# ============================================================
# 2. VALIDATE MODEL INTERFACES AND COMPARABILITY
# ============================================================

required_common<-c(
  "simulate_marl","make_transition_matrix","base_transition",
  "climate_transition_source",
  "ecological_actions","alpha","gamma","tau","memory_window","ema_alpha",
  "weather_prob","colocation_multiplier","grid_nx","grid_ny",
  "actions","r_states","stock_survival","herd_income",
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

same_transition_kernel<-function(env1,env2){
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
    "CSV transition kernel identical",
    "Stock survival identical",
    "Herd income identical"
  ),
  action_based=c(
    action_env$grid_nx*action_env$grid_ny,
    length(action_env$initial_actions),
    nrow(action_env$designated_cells),
    n_distinct(action_env$designated_cells$cell_id),
    exists("emergency_cost_scale",envir=action_env,inherits=FALSE),
    TRUE,TRUE,TRUE
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

structure_same<-
  nrow(action_env$designated_cells)==nrow(outcome_env$designated_cells)&&
  nrow(action_env$designated_cells)==nrow(no_pes_env$designated_cells)&&
  n_distinct(action_env$designated_cells$cell_id)==
    n_distinct(outcome_env$designated_cells$cell_id)&&
  n_distinct(action_env$designated_cells$cell_id)==
    n_distinct(no_pes_env$designated_cells$cell_id)&&
  exists("emergency_cost_scale",envir=action_env,inherits=FALSE)==
    exists("emergency_cost_scale",envir=outcome_env,inherits=FALSE)&&
  exists("emergency_cost_scale",envir=action_env,inherits=FALSE)==
    exists("emergency_cost_scale",envir=no_pes_env,inherits=FALSE)&&
  same_transition_kernel(action_env,outcome_env)&&
  same_transition_kernel(action_env,no_pes_env)

if(!structure_same){
  warning(
    "The three model files differ in movement, emergency-cost, or ",
    "CSV-transition structure. Model-design comparisons may be confounded."
  )
}

# ============================================================
# 3. SCENARIO GRID
# ============================================================
#
# Baseline is run once for each model design.
# Each remaining scenario changes ONE parameter from baseline.
#
# pes_multiplier scales:
#   action model  -> pes_payment
#   outcome model -> state_pes
#   No-PES model  -> not applicable; high_PES is not run for this model
#
# ============================================================

scenario_grid<-tribble(
  ~scenario,~alpha,~gamma,~tau,~memory_window,
  ~drought_prob,~normal_prob,~rainy_prob,~pes_multiplier,

  "baseline",
  0.20,0.95,1.00,5L,
  0.30,0.50,0.20,1.00,

  "high_alpha",
  0.40,0.95,1.00,5L,
  0.30,0.50,0.20,1.00,

  "low_gamma",
  0.20,0.50,1.00,5L,
  0.30,0.50,0.20,1.00,

  "high_tau",
  0.20,0.95,2.00,5L,
  0.30,0.50,0.20,1.00,

  "long_memory",
  0.20,0.95,1.00,10L,
  0.30,0.50,0.20,1.00,

  "high_drought",
  0.20,0.95,1.00,5L,
  0.50,0.35,0.15,1.00,

  "high_PES",
  0.20,0.95,1.00,5L,
  0.30,0.50,0.20,1.50
)

parameter_map<-tribble(
  ~parameter,~variant_scenario,~variant_label,
  "alpha","high_alpha","High alpha: 0.40",
  "gamma","low_gamma","Low gamma: 0.50",
  "tau","high_tau","High tau: 2.00",
  "memory_window","long_memory","Long memory: 10",
  "drought_probability","high_drought","Drought probability: 0.30 -> 0.50",
  "PES_intensity","high_PES","PES multiplier: 1.50"
)

validate_scenario_grid<-function(x){
  required<-c(
    "scenario","alpha","gamma","tau","memory_window",
    "drought_prob","normal_prob","rainy_prob","pes_multiplier"
  )
  missing<-setdiff(required,names(x))
  if(length(missing)>0){
    stop("scenario_grid is missing: ",paste(missing,collapse=", "))
  }

  weather_sum<-x$drought_prob+x$normal_prob+x$rainy_prob
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

# NOTE:
# The old high_grazing_impact and high_weather_impact scenarios were removed.
# They depended on action_deterioration and weather_deterioration, which no
# longer exist in the CSV-derived ecological transition model. Sensitivity
# tests on the transition kernel should be defined directly on the empirical
# transition probabilities rather than recreating those retired parameters.

# ============================================================
# 4. STORE BASELINE PARAMETERS FOR EACH MODEL
# ============================================================

snapshot_model<-function(env,pes_design){
  out<-list(
    alpha=env$alpha,
    gamma=env$gamma,
    tau=env$tau,
    memory_window=env$memory_window,
    ema_alpha=env$ema_alpha,
    weather_prob=env$weather_prob,
    base_transition=env$base_transition,
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

# ============================================================
# 5. RESTORE AND APPLY ONE SCENARIO
# ============================================================

restore_model<-function(env,baseline,pes_design){
  env$alpha<-baseline$alpha
  env$gamma<-baseline$gamma
  env$tau<-baseline$tau
  env$memory_window<-baseline$memory_window
  env$ema_alpha<-baseline$ema_alpha
  env$weather_prob<-baseline$weather_prob
  env$base_transition<-baseline$base_transition
  env$colocation_multiplier<-baseline$colocation_multiplier

  if(!is.null(baseline$emergency_cost_scale)){
    env$emergency_cost_scale<-baseline$emergency_cost_scale
  }

  if(pes_design=="Action-based PES"){
    env$pes_payment<-baseline$pes_vector
  }else if(pes_design=="Outcome-based PES"){
    env$state_pes<-baseline$pes_vector
  }else{
    env$pes_payment<-baseline$pes_vector
  }

  invisible(NULL)
}

apply_scenario<-function(env,baseline,pes_design,scenario_row){
  if(nrow(scenario_row)!=1L){
    stop("scenario_row must have exactly one row.")
  }

  restore_model(env,baseline,pes_design)

  env$alpha<-scenario_row$alpha[[1]]
  env$gamma<-scenario_row$gamma[[1]]
  env$tau<-scenario_row$tau[[1]]

  env$memory_window<-as.integer(scenario_row$memory_window[[1]])
  env$ema_alpha<-2/(env$memory_window+1)

  env$weather_prob<-c(
    drought=scenario_row$drought_prob[[1]],
    normal=scenario_row$normal_prob[[1]],
    rainy=scenario_row$rainy_prob[[1]]
  )

  if(pes_design=="Action-based PES"){
    env$pes_payment<-baseline$pes_vector*scenario_row$pes_multiplier[[1]]
  }else if(pes_design=="Outcome-based PES"){
    env$state_pes<-baseline$pes_vector*scenario_row$pes_multiplier[[1]]
  }else{
    #No-PES remains exactly zero under every scenario.
    env$pes_payment<-baseline$pes_vector
  }

  invisible(NULL)
}

# ============================================================
# 6. RUN ONE SCENARIO FOR ONE PES DESIGN
# ============================================================

run_scenario<-function(
    env,
    baseline,
    pes_design,
    scenario_row,
    n_runs=20L,
    steps=500L,
    seed_start=50000L
){
  scenario_name<-scenario_row$scenario[[1]]
  message("\n",pes_design," | ",scenario_name)
  apply_scenario(env,baseline,pes_design,scenario_row)

  runs<-map(seq_len(n_runs),function(run_id){
    message("  Run ",run_id," / ",n_runs)
    #Same seed sequence is used across all model designs and scenarios.
    res<-env$simulate_marl(
      steps=steps,
      seed=seed_start+run_id-1L
    )

    landscape<-res$landscape_counts%>%
      mutate(
        pes_design=pes_design,
        scenario=scenario_name,
        run=run_id,
        prop_degraded=degraded/(env$grid_nx*env$grid_ny)
      )

    cells<-res$cell_ledger%>%
      mutate(
        pes_design=pes_design,
        scenario=scenario_name,
        run=run_id
      )

    actions<-res$action_log%>%
      mutate(
        pes_design=pes_design,
        scenario=scenario_name,
        run=run_id
      )

    if("pes_reward"%in%names(actions)){
      actions$pes_signal<-actions$pes_reward
    }else if("state_pes_reward"%in%names(actions)){
      actions$pes_signal<-actions$state_pes_reward
    }else{
      actions$pes_signal<-NA_real_
    }

    if(!is.null(res$q_log)){
      q<-res$q_log%>%
        mutate(
          pes_design=pes_design,
          scenario=scenario_name,
          run=run_id
        )
    }else{
      q<-tibble()
    }

    list(landscape=landscape,cells=cells,actions=actions,q=q)
  })

  list(
    landscape=map_dfr(runs,"landscape"),
    cells=map_dfr(runs,"cells"),
    actions=map_dfr(runs,"actions"),
    q=map_dfr(runs,"q")
  )
}

# ============================================================
# 7. RUN ALL SCENARIOS FOR ONE PES DESIGN
# ============================================================

run_model_grid<-function(
    env,
    baseline,
    pes_design,
    scenario_grid,
    n_runs=20L,
    steps=500L,
    seed_start=50000L
){
  validate_scenario_grid(scenario_grid)
  on.exit(restore_model(env,baseline,pes_design),add=TRUE)

  scenario_results<-map(seq_len(nrow(scenario_grid)),function(i){
    run_scenario(
      env=env,
      baseline=baseline,
      pes_design=pes_design,
      scenario_row=scenario_grid[i,,drop=FALSE],
      n_runs=n_runs,
      steps=steps,
      seed_start=seed_start
    )
  })

  list(
    landscape=map_dfr(scenario_results,"landscape"),
    cells=map_dfr(scenario_results,"cells"),
    actions=map_dfr(scenario_results,"actions"),
    q=map_dfr(scenario_results,"q")
  )
}

# ============================================================
# 8. RUN ALL THREE MODEL DESIGNS
# ============================================================

run_three_model_designs<-function(
    scenario_grid,
    n_runs=20L,
    steps=500L,
    seed_start=50000L
){
  action_results<-run_model_grid(
    env=action_env,
    baseline=action_baseline,
    pes_design="Action-based PES",
    scenario_grid=scenario_grid,
    n_runs=n_runs,
    steps=steps,
    seed_start=seed_start
  )

  outcome_results<-run_model_grid(
    env=outcome_env,
    baseline=outcome_baseline,
    pes_design="Outcome-based PES",
    scenario_grid=scenario_grid,
    n_runs=n_runs,
    steps=steps,
    seed_start=seed_start
  )

  # PES intensity has no meaning in a model with no PES.
  no_pes_grid<-scenario_grid%>%filter(scenario!="high_PES")

  no_pes_results<-run_model_grid(
    env=no_pes_env,
    baseline=no_pes_baseline,
    pes_design="No PES",
    scenario_grid=no_pes_grid,
    n_runs=n_runs,
    steps=steps,
    seed_start=seed_start
  )

  list(
    landscape=bind_rows(
      action_results$landscape,
      outcome_results$landscape,
      no_pes_results$landscape
    ),
    cells=bind_rows(
      action_results$cells,
      outcome_results$cells,
      no_pes_results$cells
    ),
    actions=bind_rows(
      action_results$actions,
      outcome_results$actions,
      no_pes_results$actions
    ),
    q=bind_rows(
      action_results$q,
      outcome_results$q,
      no_pes_results$q
    ),
    scenarios=scenario_grid,
    parameter_map=parameter_map,
    structure_check=model_structure_check
  )
}

# ============================================================
# 9. TEST RUN
# ============================================================
#
# Use this first.
#
#comparison_test<-run_three_model_designs(
#   scenario_grid=scenario_grid,
#   n_runs=3L,
#   steps=100L,
#   seed_start=50000L
# )
#
# stopifnot(n_distinct(comparison_test$landscape$pes_design)==3)
# stopifnot(n_distinct(comparison_test$landscape$scenario)==nrow(scenario_grid))
# stopifnot(n_distinct(comparison_test$landscape$run)==3)

# ============================================================
# 10. FULL RUN
# ============================================================

N_RUNS<-20L
N_STEPS<-500L
SEED_START<-50000L

comparison_results<-run_three_model_designs(
  scenario_grid=scenario_grid,
  n_runs=N_RUNS,
  steps=N_STEPS,
  seed_start=SEED_START
)

# ============================================================
# 11. PARAMETER LOOKUP
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
# 12. VISUAL A: LANDSCAPE DEGRADATION
#     ACTION vs OUTCOME vs NO PES SIDE-BY-SIDE
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
    #geom_line(
      #data=plot_data,
      #aes(
        #x=step,
        #y=prop_degraded,
        #colour=setting,
        #group=interaction(setting,run)
      #),
      #alpha=0.15,
      #linewidth=0.30
    #)+
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
      subtitle=if(
        parameter_name=="PES_intensity"
      ){
        "PES payment increased by 50% under two designs"
      }else{
        ""
      },
      x="Iteration",
      y="Proportion of degraded cells",
      colour="Setting",
      fill="Setting"
    )+
    theme_minimal()+
    theme(legend.position="bottom")
}

#Examples:
make_parameter_degradation_plot(comparison_results,"drought_probability")
make_parameter_degradation_plot(comparison_results,"PES_intensity")
make_parameter_degradation_plot(comparison_results,"tau")

# ============================================================
# 13. PRINT ONE FACET PLOT FOR EVERY PARAMETER
# ============================================================

parameter_degradation_plots<-setNames(
  map(
    parameter_map$parameter,
    ~make_parameter_degradation_plot(comparison_results,.x)
  ),
  parameter_map$parameter
)

walk(parameter_degradation_plots,print)

#Individual plots are available as:
parameter_degradation_plots$alpha
parameter_degradation_plots$gamma
parameter_degradation_plots$tau
parameter_degradation_plots$memory_window
parameter_degradation_plots$drought_probability
#parameter_degradation_plots$grazing_impact
#parameter_degradation_plots$weather_impact
parameter_degradation_plots$PES_intensity

# ============================================================
# 14. VISUAL B: CELL-LEVEL EMA
#     ACTION vs OUTCOME vs NO PES SIDE-BY-SIDE
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

  ggplot(
    plot_data,
    aes(x=step,y=cell_id,fill=ema_after)
  )+
    geom_tile()+
    scale_fill_viridis_c(
      limits=c(0,3),
      breaks=c(0,0.5,1.5,2.5,3)
    )+
    facet_grid(setting~pes_design)+
    labs(
      title=paste("Cell-level EMA sensitivity to",parameter_name),
      subtitle=paste("Run",run_id,"; columns compare model design"),
      x="Iteration",
      y="Cell ID",
      fill="EMA degradation"
    )+
    theme_minimal()
}

#Example:
plot_cell_degradation_comparison(
  comparison_results,
  parameter_name="memory_window",
  run_id=1
 )

# ============================================================
# 15. VISUAL C: ACTION MIX
#     ACTION vs OUTCOME vs NO PES SIDE-BY-SIDE
# ============================================================

plot_action_mix_comparison<-function(results,parameter_name){
  plot_data<-prepare_parameter_data(results$actions,parameter_name)%>%
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
      subtitle="Columns compare model design; rows compare baseline with the altered parameter",
      x="Iteration",
      y="Proportion of agent actions",
      fill="Stocking action"
    )+
    theme_minimal()
}

#Example:
plot_action_mix_comparison(comparison_results,"tau")

# ============================================================
# 16. FINAL-STEP SUMMARY TABLE
# ============================================================

summarise_final_landscape<-function(results){
  last_step<-max(results$landscape$step)
  results$landscape%>%
    filter(step==last_step)%>%
    group_by(pes_design,scenario)%>%
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

# ============================================================
# 17. USEFUL RESULT OBJECTS
# ============================================================
#
comparison_results$landscape
# comparison_results$cells
comparison_results$actions
# comparison_results$q
# comparison_results$scenarios
# comparison_results$parameter_map
# comparison_results$structure_check
#
# Main plots:
# parameter_degradation_plots
# plot_cell_degradation_comparison(...)
# plot_action_mix_comparison(...)
#
# ============================================================
# ============================================================
# 17. LOOKING AT PARAMETERS SEPARATELY for cumulative impact
# ============================================================

landscape_auc<-landscape_auc%>%
  mutate(
    scenario=factor(
      scenario,
      levels=c(
        "baseline",
        "high_alpha",
        "low_gamma",
        "high_tau",
        "long_memory",
        "high_drought",
        "high_PES"
      )
    ),
    pes_design=factor(
      pes_design,
      levels=c(
        "Action-based PES",
        "Outcome-based PES",
        "No PES"
      ),
      labels=c(
        "Action",
        "Outcome",
        "None"
      )
    )
  )

landscape_auc%>%
  ggplot(
    aes(
      x=scenario,
      y=normalized_auc,
      fill=pes_design
    )
  )+
  geom_boxplot(
    position=position_dodge(width=0.8),
    width=0.7,
    outlier.shape=NA
  )+
  geom_point(
    aes(
      group=pes_design
    ),
    position=position_jitterdodge(
      dodge.width=0.8,
      jitter.width=0.08
    ),
    alpha=0.45,
    size=1.5
  )+
  labs(
    title="Cumulative landscape degradation across sensitivity scenarios",
    x="Scenario",
    y="Normalized AUC (%)",
    fill="PES design"
  )+
  theme_minimal()+
  theme(
    axis.text.x=element_text(
      angle=45,
      hjust=1
    )
  )