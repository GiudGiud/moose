# =====================================================================
# Burns, Sims & Goodall (1984), Radiat. Phys. Chem. 23(1/2), 143-180
# Spherical spur, Schwarz parameters, Gaussian initial conditions.
# Reproduces the FACSIMILE numerical case (Figs. 2-7).
#
# Diffusion-reaction of 7 radiolytic species:
#   e_aq, Hp(=H+), H, OH, OHm(=OH-), H2, H2O2
#
# Consistent units: mol, dm, s.
#   concentration [mol dm^-3]; D [dm^2 s^-1]; k [dm^3 mol^-1 s^-1];
#   radius r=x [dm]  (1 nm = 1e-8 dm)
#
# Implemented entirely with existing framework objects:
#   ADTimeDerivative + ADMatDiffusion (PDE), ADMatBodyForce fed by an
#   ADParsedMaterial holding each species' net reaction rate (guide S.12).
# Zero-flux at r=0 and r=r_max are the natural BC, so no [BCs] block.
# =====================================================================

[Mesh]
  # 1D radial mesh, r in [0, 600 nm] = [0, 6e-6 dm].
  # Domain must be >> the diffusion length at t_end: sigma = sqrt(2*D*t)
  # ~ 95 nm for e_aq at 1e-6 s, so the cloud spans a few hundred nm. A too
  # small domain (e.g. 30 nm) confines the species against the reflecting
  # outer wall, causing spurious over-recombination so G never plateaus.
  # Graded mesh (bias_x=1.03) resolves the Gaussian core (dx ~ 0.003 nm at
  # r=0, where the Gaussian shell r ~ sqrt(2)*r0 ~ 3.5 nm lives) while
  # cheaply reaching 600 nm (dx ~ 17 nm in the far field). This is the role
  # of the paper's expanding K=sqrt(1.5) grid; we need more zones because
  # MOOSE uses nodal FE interpolation rather than the paper's cell-averaged
  # finite volumes. Recovers the Table-2 G values at t=0 to <0.1%.
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    xmin = 0
    xmax = 6e-6
    nx = 300
    bias_x = 1.03
  []
  coord_type = RSPHERICAL
[]

[Variables]
  [e_aq]
  []
  [Hp]
  []
  [H]
  []
  [OH]
  []
  [OHm]
  []
  [H2]
  []
  [H2O2]
  []
[]

# ---------------------------------------------------------------------
# Initial conditions: Gaussian per species (guide eq. 1, t=0)
#   C(r,0) = (N/Na) / (pi*2*r0^2)^(3/2) * exp(-r^2/(2*r0^2))
# OHm and H2O2 start at zero (produced only by reactions).
# ---------------------------------------------------------------------
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

[ICs]
  [e_aq]
    type = FunctionIC
    variable = e_aq
    function = ic_e_aq
  []
  [Hp]
    type = FunctionIC
    variable = Hp
    function = ic_Hp
  []
  [H]
    type = FunctionIC
    variable = H
    function = ic_H
  []
  [OH]
    type = FunctionIC
    variable = OH
    function = ic_OH
  []
  [H2]
    type = FunctionIC
    variable = H2
    function = ic_H2
  []
  # OHm, H2O2 default to 0
[]

# ---------------------------------------------------------------------
# Kernels: dC/dt = D * laplacian(C) + S(C)    [spherical, via RSPHERICAL]
#   ADMatBodyForce adds the residual -S*test, where the material S_x is
#   the net production rate of species x (guide S.12, corrected S_OHm).
# ---------------------------------------------------------------------
[Kernels]
  # ---- e_aq ----
  [e_aq_time]
    type = ADTimeDerivative
    variable = e_aq
  []
  [e_aq_diff]
    type = ADMatDiffusion
    variable = e_aq
    diffusivity = D_e_aq
  []
  [e_aq_rxn]
    type = ADMatBodyForce
    variable = e_aq
    material_property = S_e_aq
  []

  # ---- Hp ----
  [Hp_time]
    type = ADTimeDerivative
    variable = Hp
  []
  [Hp_diff]
    type = ADMatDiffusion
    variable = Hp
    diffusivity = D_Hp
  []
  [Hp_rxn]
    type = ADMatBodyForce
    variable = Hp
    material_property = S_Hp
  []

  # ---- H ----
  [H_time]
    type = ADTimeDerivative
    variable = H
  []
  [H_diff]
    type = ADMatDiffusion
    variable = H
    diffusivity = D_H
  []
  [H_rxn]
    type = ADMatBodyForce
    variable = H
    material_property = S_H
  []

  # ---- OH ----
  [OH_time]
    type = ADTimeDerivative
    variable = OH
  []
  [OH_diff]
    type = ADMatDiffusion
    variable = OH
    diffusivity = D_OH
  []
  [OH_rxn]
    type = ADMatBodyForce
    variable = OH
    material_property = S_OH
  []

  # ---- OHm ----
  [OHm_time]
    type = ADTimeDerivative
    variable = OHm
  []
  [OHm_diff]
    type = ADMatDiffusion
    variable = OHm
    diffusivity = D_OHm
  []
  [OHm_rxn]
    type = ADMatBodyForce
    variable = OHm
    material_property = S_OHm
  []

  # ---- H2 ----
  [H2_time]
    type = ADTimeDerivative
    variable = H2
  []
  [H2_diff]
    type = ADMatDiffusion
    variable = H2
    diffusivity = D_H2
  []
  [H2_rxn]
    type = ADMatBodyForce
    variable = H2
    material_property = S_H2
  []

  # ---- H2O2 ----
  [H2O2_time]
    type = ADTimeDerivative
    variable = H2O2
  []
  [H2O2_diff]
    type = ADMatDiffusion
    variable = H2O2
    diffusivity = D_H2O2
  []
  [H2O2_rxn]
    type = ADMatBodyForce
    variable = H2O2
    material_property = S_H2O2
  []
[]

# ---------------------------------------------------------------------
# Materials
# ---------------------------------------------------------------------
[Materials]
  # Diffusion coefficients (Schwarz, Table 1), units dm^2 s^-1.
  # D(H2) is not given for Schwarz; use 5.0e-7 per guide S.11.
  [diffusivities]
    type = ADGenericConstantMaterial
    prop_names =  'D_e_aq  D_Hp    D_H     D_OH    D_OHm   D_H2    D_H2O2'
    prop_values = '4.5e-7  9.0e-7  7.0e-7  2.8e-7  5.0e-7  5.0e-7  2.2e-7'
  []

  # Net reaction rates S_x (guide S.12). Rate constants [dm^3 mol^-1 s^-1]:
  #   k1=5.5e9 k2=2.5e10 k3=1.0e10 k4=1.7e10 k5=2.5e10
  #   k6=1.3e10 k7=1.0e8 k8=1.0e11 k9=6.0e9 k10=2.0e10
  # Factor of 2 on like-species reactions (1, 3, 9) per Table 3 note a.
  [S_e_aq]
    type = ADParsedMaterial
    property_name = S_e_aq
    coupled_variables = 'e_aq Hp H OH H2O2'
    constant_names =       'k1    k2     k4     k5     k6'
    constant_expressions = '5.5e9 2.5e10 1.7e10 2.5e10 1.3e10'
    expression = '-2*k1*e_aq*e_aq - k2*H*e_aq - k4*e_aq*Hp - k5*e_aq*OH - k6*e_aq*H2O2'
  []
  [S_Hp]
    type = ADParsedMaterial
    property_name = S_Hp
    coupled_variables = 'e_aq Hp OHm'
    constant_names =       'k4     k8'
    constant_expressions = '1.7e10 1.0e11'
    expression = '-k4*e_aq*Hp - k8*Hp*OHm'
  []
  [S_H]
    type = ADParsedMaterial
    property_name = S_H
    coupled_variables = 'e_aq Hp H OH H2O2'
    constant_names =       'k2     k3     k4     k7    k10'
    constant_expressions = '2.5e10 1.0e10 1.7e10 1.0e8 2.0e10'
    expression = 'k4*e_aq*Hp - 2*k3*H*H - k2*H*e_aq - k7*H*H2O2 - k10*H*OH'
  []
  [S_OH]
    type = ADParsedMaterial
    property_name = S_OH
    coupled_variables = 'e_aq H OH H2O2'
    constant_names =       'k5     k6     k7    k9    k10'
    constant_expressions = '2.5e10 1.3e10 1.0e8 6.0e9 2.0e10'
    expression = 'k6*e_aq*H2O2 + k7*H*H2O2 - k5*e_aq*OH - 2*k9*OH*OH - k10*H*OH'
  []
  [S_OHm]
    type = ADParsedMaterial
    property_name = S_OHm
    coupled_variables = 'e_aq Hp H OH OHm H2O2'
    constant_names =       'k1    k2     k5     k6     k8'
    constant_expressions = '5.5e9 2.5e10 2.5e10 1.3e10 1.0e11'
    expression = '2*k1*e_aq*e_aq + k2*H*e_aq + k5*e_aq*OH + k6*e_aq*H2O2 - k8*Hp*OHm'
  []
  [S_H2]
    type = ADParsedMaterial
    property_name = S_H2
    coupled_variables = 'e_aq H'
    constant_names =       'k1    k2     k3'
    constant_expressions = '5.5e9 2.5e10 1.0e10'
    expression = 'k1*e_aq*e_aq + k2*H*e_aq + k3*H*H'
  []
  [S_H2O2]
    type = ADParsedMaterial
    property_name = S_H2O2
    coupled_variables = 'e_aq H OH H2O2'
    constant_names =       'k6     k7    k9'
    constant_expressions = '1.3e10 1.0e8 6.0e9'
    expression = 'k9*OH*OH - k6*e_aq*H2O2 - k7*H*H2O2'
  []
[]

# ---------------------------------------------------------------------
# Postprocessors: G value of species i = N_i(t) * 100 / E_spur, with
#   N_i = Na * integral C_i dV.
# In RSPHERICAL 1D, ElementIntegralVariablePostprocessor integrates over
# the full spherical volume element (it carries the 4*pi*r^2 Jacobian),
# so it returns integral(C 4*pi*r^2 dr) and N_i = Na * int (no extra 4*pi).
# E_spur = 62.5 eV. (Sanity: at t=0, G recovers Table 2 values.)
# ---------------------------------------------------------------------
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
  [int_Hp]
    type = ElementIntegralVariablePostprocessor
    variable = Hp
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [int_OHm]
    type = ElementIntegralVariablePostprocessor
    variable = OHm
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
  [G_Hp]
    type = ParsedPostprocessor
    pp_names = 'int_Hp'
    constant_names =       'Na       Espur'
    constant_expressions = '6.022e23 62.5'
    expression = 'Na*int_Hp*100/Espur'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [G_OHm]
    type = ParsedPostprocessor
    pp_names = 'int_OHm'
    constant_names =       'Na       Espur'
    constant_expressions = '6.022e23 62.5'
    expression = 'Na*int_OHm*100/Espur'
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
  scheme = bdf2

  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'

  # Stiff system (k8 = 1e11); concentrations span many orders -> autoscale.
  automatic_scaling = true

  # IMPORTANT: with dm/nm units the absolute residuals are tiny, and
  # automatic_scaling shrinks the reaction residual ~ dt (the 1/dt time
  # derivative dominates the scaling). A normal nl_abs_tol would be met at
  # iteration 0, so the solver "converges" without ever applying the
  # reaction (worse at smaller dt). Drive nl_abs_tol to ~0 so the relative
  # tolerance binds and the reaction is always integrated.
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-30
  l_max_its = 50

  end_time = 1e-6
  dtmin = 1e-18
  # IterationAdaptiveDT grows dt from the iteration count, not a temporal
  # error estimate; cap dt (dtmax) so the fast early recombination
  # (t < 1e-8 s) stays time-resolved. Result is dt-converged at these
  # settings (matches the paper's error-controlled Gear integration).
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-14
    growth_factor = 1.2
    cutback_factor = 0.5
    optimal_iterations = 8
  []
  dtmax = 1e-8
[]

[Outputs]
  exodus = true   # concentration profiles C(r,t) for Figs. 2-5
  csv = true      # G(t) curves for Figs. 6-7
[]
