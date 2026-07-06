# =====================================================================
# Burns, Sims & Goodall (1984) spherical-spur water-radiolysis case WITH
# diffusing scavengers (Appleby & Schwarz conditions, guide reactions 11-13),
# assembled with the KineticReactionsPhysics.
#
# The scavengers (and their products) are untracked bulk species held at a
# fixed uniform concentration, so each scavenging reaction
#     spur_species + S -> product
# reduces to a pseudo-first-order sink on the spur species with effective rate
#     k_eff = k_bimolecular * [S].
#   e_aq:  0.87e10 * 7e-4 = 6.09e6 1/s   ([S1] = 7e-4 mol/dm^3)
#   H:     0.05e10 * 1e-2 = 5.0e6  1/s   ([S2] = 1e-2 mol/dm^3)
#   OH:    0.13e10 * 1e-2 = 1.3e7  1/s   ([S3] = 1e-2 mol/dm^3)
#
# Units: mol, dm, s.  r = x in dm (1 nm = 1e-8 dm).
# =====================================================================

# TODO: find diffusivity of H2

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    xmin = 0
    xmax = 6e-6      # 600 nm
    nx = 300
    bias_x = 1.03
  []
  coord_type = RSPHERICAL
[]

# Gaussian initial conditions C(r,0) = (N/Na)/(pi*2*r0^2)^1.5 * exp(-r^2/(2*r0^2))
[Functions]
  [ic_e_aq]
    type = ParsedFunction
    expression = '(N/Na)/(pi*2*r0^2)^1.5 * exp(-x^2/(2*r0^2))'
    symbol_names = 'N Na r0'
    symbol_values = '2.9875 6.022e23 2.458e-8'
  []
  [ic_Hp]
    type = ParsedFunction
    expression = '(N/Na)/(pi*2*r0^2)^1.5 * exp(-x^2/(2*r0^2))'
    symbol_names = 'N Na r0'
    symbol_values = '2.9875 6.022e23 1.145e-8'
  []
  [ic_H]
    type = ParsedFunction
    expression = '(N/Na)/(pi*2*r0^2)^1.5 * exp(-x^2/(2*r0^2))'
    symbol_names = 'N Na r0'
    symbol_values = '0.3875 6.022e23 1.145e-8'
  []
  [ic_OH]
    type = ParsedFunction
    expression = '(N/Na)/(pi*2*r0^2)^1.5 * exp(-x^2/(2*r0^2))'
    symbol_names = 'N Na r0'
    symbol_values = '3.5625 6.022e23 1.145e-8'
  []
  [ic_H2]
    type = ParsedFunction
    expression = '(N/Na)/(pi*2*r0^2)^1.5 * exp(-x^2/(2*r0^2))'
    symbol_names = 'N Na r0'
    symbol_values = '0.09375 6.022e23 1.145e-8'
  []
[]

[Physics]
  [KineticReactions/network]
    species = 'e_aq Hp H OH OHm H2 H2O2'
    # OHm and H2O2 start at zero (constant-0 functor)
    initial_conditions = 'ic_e_aq ic_Hp ic_H ic_OH 0 ic_H2 0'
    diffusivities = '4.5e-7 9.0e-7 7.0e-7 2.8e-7 5.0e-7 5.0e-7 2.2e-7'
    order = FIRST

    # Schwarz rate constants [dm^3 mol^-1 s^-1]. H2O is the (untracked) solvent.
    # The factor-of-2 convention for like-species reactions is automatic.
    # The last three are the scavenging reactions: the untracked scavenger and
    # product give a pseudo-first-order sink with k = k_bimolecular * [scavenger].
    reactions = '2e_aq -> H2 + 2OHm   [k=5.5e9]
                 H + e_aq -> H2 + OHm [k=2.5e10]
                 2H -> H2             [k=1.0e10]
                 e_aq + Hp -> H       [k=1.7e10]
                 e_aq + OH -> OHm     [k=2.5e10]
                 e_aq + H2O2 -> OHm + OH [k=1.3e10]
                 H + H2O2 -> OH + H2O [k=1.0e8]
                 Hp + OHm -> H2O      [k=1.0e11]
                 2OH -> H2O2          [k=6.0e9]
                 H + OH -> H2O        [k=2.0e10]
                 e_aq -> Scavenged1   [k=8.7e9]
                 H -> Scavenged2      [k=5.0e8]
                 OH -> Scavenged3     [k=1.3e9]'
  []
[]

# G value of species s = Na * integral(C_s 4 pi r^2 dr) * 100 / E_spur.
# In RSPHERICAL the integral already carries the 4 pi r^2 weight.
[Postprocessors]
  [int_e_aq]
    type = ElementIntegralVariablePostprocessor
    variable = e_aq
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [int_OH]
    type = ElementIntegralVariablePostprocessor
    variable = OH
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [int_H]
    type = ElementIntegralVariablePostprocessor
    variable = H
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [int_H2]
    type = ElementIntegralVariablePostprocessor
    variable = H2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [int_H2O2]
    type = ElementIntegralVariablePostprocessor
    variable = H2O2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [G_e_aq]
    type = ParsedPostprocessor
    pp_names = 'int_e_aq'
    constant_names =       'Na       Espur'
    constant_expressions = '6.022e23 62.5'
    expression = 'Na*int_e_aq*100/Espur'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [G_OH]
    type = ParsedPostprocessor
    pp_names = 'int_OH'
    constant_names =       'Na       Espur'
    constant_expressions = '6.022e23 62.5'
    expression = 'Na*int_OH*100/Espur'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [G_H]
    type = ParsedPostprocessor
    pp_names = 'int_H'
    constant_names =       'Na       Espur'
    constant_expressions = '6.022e23 62.5'
    expression = 'Na*int_H*100/Espur'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [G_H2]
    type = ParsedPostprocessor
    pp_names = 'int_H2'
    constant_names =       'Na       Espur'
    constant_expressions = '6.022e23 62.5'
    expression = 'Na*int_H2*100/Espur'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [G_H2O2]
    type = ParsedPostprocessor
    pp_names = 'int_H2O2'
    constant_names =       'Na       Espur'
    constant_expressions = '6.022e23 62.5'
    expression = 'Na*int_H2O2*100/Espur'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  # scheme = bdf2
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'
  automatic_scaling = true
  # nl_abs_tol ~ 0 so the relative tolerance binds (dm/nm units make the
  # absolute reaction residual tiny; otherwise the solver "converges" at
  # iteration 0 without applying the reactions).
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-30
  l_max_its = 50
  end_time = 1e-6
  dtmin = 1e-18
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-15
    growth_factor = 1.2
    cutback_factor = 0.5
    optimal_iterations = 8
  []
  dtmax = 1e-8
[]

[Outputs]
  csv = true
[]
