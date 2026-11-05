# =====================================================================
# Bayesian search, with parallel Markov chain Monte Carlo (MCMC), of
# two rate constants (k1 of A -> B, k2 of B -> C) from the observed final [C]
# in the network A -> B -> C (forward_stochastic.i).
#
# One final concentration cannot determine two rate constants, so [C] is
# observed in two experiments (config_several_rates.csv): starting from [A] = 1 M,
# [C] depends on k1 and k2; starting from [B] = 1 M, it depends on k2 only.
#
# Rate constants are sampled as log10(k) with uniform priors over four decades;
# proposals outside the priors have zero prior density and are rejected.
# The synthetic observations were computed with k1 = 2 /s and k2 = 0.5 /s,
# i.e. log10(k1) = 0.301, log10(k2) = -0.301. The posterior samples are written
# to the JSON output (mcmc_reporter/inputs).
# =====================================================================

[StochasticTools]
[]

[Distributions]
  [prior_logk1]
    type = Uniform
    lower_bound = -2
    upper_bound = 2
  []
  [prior_logk2]
    type = Uniform
    lower_bound = -2
    upper_bound = 2
  []
[]

[Likelihood]
  [gaussian]
    type = Gaussian
    noise = 'noise/noise'
    file_name = 'observed_several_rates.csv'
    log_likelihood = true
  []
[]

[Samplers]
  [sample]
    type = AffineInvariantDES
    prior_distributions = 'prior_logk1 prior_logk2'
    num_parallel_proposals = 5
    file_name = 'config_several_rates.csv'
    num_columns = 2
    initial_values = '0 0'
    previous_state = 'mcmc_reporter/inputs'
    previous_state_var = 'mcmc_reporter/variance'
    execute_on = 'PRE_MULTIAPP_SETUP TIMESTEP_END'
    seed = 2547
    # The time step indexes the random seeds: needs at least num_steps + 1
    num_random_seeds = 101
  []
[]

[MultiApps]
  [forward]
    type = SamplerFullSolveMultiApp
    input_files = forward_stochastic.i
    sampler = sample
  []
[]

[Controls]
  [cmdline]
    type = MultiAppSamplerControl
    multi_app = forward
    sampler = sample
    param_names = 'logk1 logk2 experiment species'
  []
[]

[Transfers]
  [observed]
    type = SamplerReporterTransfer
    from_multi_app = forward
    sampler = sample
    stochastic_reporter = 'constant'
    from_reporter = 'observed/value'
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
    # Standard deviation of the measurement error [mol dm^-3]. Concentration misfits
    # are below 1 M, so with this noise the Gaussian log-likelihood stays above
    # -1250 and does not underflow to -inf for any sample drawn from the priors.
    real_values = '0.02'
    outputs = none
  []
  [mcmc_reporter]
    type = AffineInvariantDifferentialDecision
    output_value = constant/observed:observed:value
    sampler = sample
    likelihoods = 'gaussian'
    execute_on = 'TIMESTEP_BEGIN'
  []
[]

[Executioner]
  type = Transient
  # Number of MCMC steps of each of the 5 parallel chains
  num_steps = 100
[]

[Outputs]
  [out]
    type = JSON
    execute_system_information_on = NONE
  []
[]
