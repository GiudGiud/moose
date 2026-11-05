# =====================================================================
# Forward model for the stochastic-tools rate searches: the consecutive
# first-order network A -> B -> C in a single 0-D element, run to t = 1 s.
#
# The rate constants enter as log10(k) so that a search can span several
# decades. The synthetic "observations" were generated from this input with
# the reference values k1 = 2 /s and k2 = 0.5 /s (the defaults below).
#
# Top-level variables set by MultiAppSamplerControl:
#   logk1, logk2 : log10 of the rate constants of A -> B and B -> C
#   species      : observed species, 0 = A, 1 = B, 2 = C
#   experiment   : 0 starts from [A] = 1 M, 1 starts from [B] = 1 M
# =====================================================================

logk1 = 0.301029995663981 # log10(2)
logk2 = -0.301029995663981 # log10(0.5)
species = 0
experiment = 0

[Mesh]
  # Single element (0-D homogeneous kinetics)
  [single]
    type = ElementGenerator
    nodal_positions = '0 0 0
                       1 0 0'
    element_connectivity = '0 1'
    elem_type = EDGE2
  []
[]

[Functions]
  [k1]
    type = ParsedFunction
    expression = '10^${logk1}'
  []
  [k2]
    type = ParsedFunction
    expression = '10^${logk2}'
  []
[]

[Physics]
  [KineticReactions/network]
    species = 'A B C'
    diffusivities = '0 0 0'
    initial_conditions = '${fparse 1 - experiment} ${experiment} 0'
    reactions = 'A -> B [k=k1]
                 B -> C [k=k2]'
  []
[]

[Postprocessors]
  [A]
    type = ElementAverageValue
    variable = A
    outputs = none
  []
  [B]
    type = ElementAverageValue
    variable = B
    outputs = none
  []
  [C]
    type = ElementAverageValue
    variable = C
    outputs = none
  []
  # Final concentration of the observed species, returned to the sampler
  [observed]
    type = ParsedPostprocessor
    expression = 'if(species < 0.5, A, if(species < 1.5, B, C))'
    pp_names = 'A B C'
    constant_names = 'species'
    constant_expressions = '${species}'
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = NEWTON
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  dt = 0.01
  end_time = 1
[]

[Outputs]
  console = false
[]
