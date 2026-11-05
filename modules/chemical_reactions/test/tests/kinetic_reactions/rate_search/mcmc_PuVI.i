# =====================================================================
# Bayesian search, with parallel Markov chain Monte Carlo (MCMC), of one rate
# constant of the plutonium radiolysis network (forward_PuVI.i), that of
#   PuO22p + H -> PuO2p + Hp
# from the observed final log10([PuO2p]).
#
# The rate constant is sampled as log10(k) with a uniform prior over about one
# decade on either side of the reference value k = 2e8 dm^3 mol^-1 s^-1 (log10(k) = 8.301)
# used to compute the synthetic observation. The posterior samples are written
# to the JSON output (mcmc_reporter/inputs).
# =====================================================================

[StochasticTools]
[]

[Distributions]
  [prior_logk]
    type = Uniform
    # A wider prior lets the initial walkers spread so far apart that one can
    # stall: the differential-evolution proposals are built from differences
    # between the other walkers, which shrink once those have converged
    lower_bound = 7
    upper_bound = 9.5
  []
[]

[Likelihood]
  [gaussian]
    type = Gaussian
    noise = 'noise/noise'
    file_name = 'observed_PuVI.csv'
    log_likelihood = true
  []
[]

[Samplers]
  [sample]
    type = AffineInvariantDES
    prior_distributions = 'prior_logk'
    num_parallel_proposals = 5
    # Experimental configuration: the observation time [s]
    file_name = 'config_PuVI.csv'
    initial_values = '8'
    previous_state = 'mcmc_reporter/inputs'
    previous_state_var = 'mcmc_reporter/variance'
    execute_on = 'PRE_MULTIAPP_SETUP TIMESTEP_END'
    seed = 2547
    # The time step indexes the random seeds: needs at least num_steps + 1
    num_random_seeds = 51
  []
[]

[MultiApps]
  [forward]
    type = SamplerFullSolveMultiApp
    input_files = forward_PuVI.i
    sampler = sample
  []
[]

[Controls]
  [cmdline]
    type = MultiAppSamplerControl
    multi_app = forward
    sampler = sample
    param_names = 'logk Executioner/end_time'
  []
[]

[Transfers]
  [observed]
    type = SamplerReporterTransfer
    from_multi_app = forward
    sampler = sample
    stochastic_reporter = 'constant'
    from_reporter = 'log10_PuO2p/value'
  []
[]

[Reporters]
  [constant]
    type = StochasticReporter
    outputs = none
  []
  [noise]
    type = ConstantReporter
    real_names = 'noise'
    # Standard deviation of the measured log10([PuO2p]), i.e. a ~5% relative error
    real_values = '0.02'
    outputs = none
  []
  [mcmc_reporter]
    type = AffineInvariantDifferentialDecision
    output_value = constant/observed:log10_PuO2p:value
    sampler = sample
    likelihoods = 'gaussian'
    execute_on = 'TIMESTEP_BEGIN'
  []
[]

[Executioner]
  type = Transient
  # Number of MCMC steps of each of the 5 parallel chains
  num_steps = 50
[]

[Outputs]
  [out]
    type = JSON
    execute_system_information_on = NONE
  []
[]
