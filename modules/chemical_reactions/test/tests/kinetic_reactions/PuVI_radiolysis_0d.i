# =====================================================================
# Plutonium radiolysis in nitric acid (Conrad, Donoclift, Holmbeck & Kynman 2026),
# reaction network extracted from PuVI_2026.fac section 6.3, assembled with the
# KineticReactionsPhysics on a single-element (0-D) mesh -- no diffusion.
#
#   - Rate constants evaluated at T = 20 C (PuVIGamma1M_2025 data file), gamma 1M HNO3.
#   - H2O is the untracked solvent: dropped from reactants (its concentration
#     [H2O] = 52.067 M folded into k), kept as an untracked product token.
#   - FACSIMILE rate-tracking pseudo-variables (+R0..R242, +Y0..Y30) stripped.
#   - The dose-driven primary-species generation (sec 6.2) is NOT a constant-rate
#     reaction; instead the primary radicals are seeded as initial conditions from
#     the file's gamma G-values x a 100 Gy pulse, then the network relaxes.
#   - Charged names mapped to valid identifiers: '+' -> p, '-' -> m
#     (e.g. H+ -> Hp, NO3- -> NO3m, Eaq- -> Eaqm, PuO22+ -> PuO22p).
#
# Units: mol, dm^3, s.  (47 species, 157 reactions)
# =====================================================================

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

[Physics]
  [KineticReactions/PuVI]
    species = 'AmII AmIII OHm H Eaqm OH AmIV NO3 NO3m PuIV PuIII PuO2p Hp PuO22p Hm H2O2 N2O4 NO2 NO2m HO2 HNO2 PuO23p HNO3 O O2 Eprem NO32m HNO3m O2m H2NO3 OONOOm HOONOO ONOOH NO N2O3 N2O5 HNO2m H2NO2 Om ONOOm NO22m N2O6 N2O N2 H2 HO2m O3m'
    # no diffusion: all diffusivities zero (Physics omits the diffusion kernel)
    diffusivities = '0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0'
    order = FIRST

    # Initial conditions [mol dm^-3]: bulk acid/Pu/O2 from the data file; primary
    # radicals seeded from gamma G-values (sec 6.2) x 100 Gy. Others start at 0.
    initial_conditions = '0 0 1e-07 9.28794e-07 0 1.50908e-05 0 1.61139e-05 0 0 0 0 1e-07 0.000359 0 6.84052e-06 0 0 0 4.35372e-08 0 0 1 1.82442e-07 0.00025 0 3.73953e-05 2.68894e-06 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.17136e-06 0 0'

    # Reaction network (PuVI_2026.fac sec 6.3), forward rate constants at 20 C.
    reactions = 'AmII -> AmIII + OHm + H [k=4.89431e6]
                 AmIII + Eaqm -> AmII [k=1.6e8]
                 AmIII + OH -> AmIV + OHm [k=3e8]
                 AmIII + NO3 -> AmIV + NO3m [k=1.32e8]
                 2PuIV -> PuIII + PuO2p + 4Hp [k=1e-6]
                 PuIII + PuO2p + 4Hp -> 2PuIV [k=1]
                 PuIV + PuO2p -> PuIII + PuO22p [k=13.228]
                 PuIII + PuO22p -> PuIV + PuO2p [k=1]
                 2PuO2p -> PuIV + PuO22p [k=0.01]
                 PuIV + PuO22p -> 2PuO2p + 4Hp [k=2.71098e9]
                 PuIII + H -> PuIV + Hm [k=1e6]
                 PuIII + OH -> PuIV + OHm [k=4.2e8]
                 PuIII + NO3 -> PuIV + NO3m [k=2.5e8]
                 PuIII + H2O2 -> PuIV + OH + OHm [k=0.055]
                 PuIII + N2O4 -> PuIV + NO2 + NO2m [k=0.001]
                 PuIV + H -> PuIII + Hp [k=2e7]
                 PuIV + Eaqm -> PuIII [k=3e10]
                 PuIV + NO3 -> PuO2p + NO3m + 4Hp [k=10000]
                 PuIV + H2O2 -> PuIII + HO2 + Hp [k=0.1]
                 PuIV + HNO2 -> PuIII + NO2 + Hp [k=0.033]
                 PuO2p + NO3 -> PuO22p + NO3m [k=1e7]
                 PuO2p + NO2 + Hp -> PuO22p + HNO2 [k=19000]
                 PuO22p + H -> PuO2p + Hp [k=2e8]
                 PuO22p + Eaqm -> PuO2p [k=6.1e10]
                 PuO22p + H2O2 -> PuO2p + HO2 + Hp [k=0.0062]
                 PuO22p + HNO2 -> PuO2p + NO2 + Hp [k=0.03]
                 PuO22p + NO3 -> PuO23p + NO3m [k=5.97e8]
                 PuO23p + Eaqm -> PuO22p [k=3.5e10]
                 PuO23p + HNO3 -> PuO22p + NO3 + Hp [k=10000]
                 PuO23p -> PuO22p + OH + Hp [k=520671]
                 NO3m + O -> NO2m + O2 [k=1e10]
                 HNO3 -> Hp + NO3m [k=1.46e10]
                 HNO3 + OH -> NO3 + H2O [k=1.9e7]
                 NO3m + Eprem -> NO32m [k=1e13]
                 NO3m + Eaqm -> NO32m [k=9.7e9]
                 NO3m + H -> HNO3m [k=3.99186e6]
                 NO3m + Hp -> HNO3 [k=6e8]
                 NO32m -> NO2 + 2OHm [k=52067.1]
                 NO32m + OH -> NO3m + OHm [k=3e9]
                 NO32m + H2O2 -> NO3m + OH + OHm [k=1e8]
                 NO32m + Hp -> HNO3m [k=5e10]
                 NO32m + O2 -> NO3m + O2m [k=2.4e8]
                 HNO3m + Hp -> H2NO3 [k=5e10]
                 HNO3m -> NO32m + Hp [k=1600]
                 HNO3m -> NO2 + OHm [k=200000]
                 HNO3m + O2 -> NO3m + O2m + Hp [k=1e8]
                 H2NO3 -> NO2 + H2O [k=700000]
                 H2NO3 -> HNO3m + Hp [k=1600]
                 H2NO3 + O2 -> HO2 + HNO3 [k=1e6]
                 2NO2 -> N2O4 [k=4.5e8]
                 NO2 + O2m -> OONOOm [k=4.5e9]
                 NO2 + HO2 -> HOONOO [k=1.8e9]
                 NO2 + OH -> ONOOH [k=1e10]
                 NO2 + H -> HNO2 [k=1e10]
                 NO2 + NO -> N2O3 [k=1.1e9]
                 NO2 + Eaqm -> NO2m [k=1e10]
                 NO2 + NO3 -> N2O5 [k=1.1e9]
                 N2O5 -> 2HNO3 [k=1874.42]
                 N2O4 -> HNO2 + HNO3 [k=937.208]
                 N2O4 -> 2NO2 [k=6000]
                 HNO2 -> NO2m + Hp [k=4.6e6]
                 HNO2 + OH -> NO2 + H2O [k=2.9e9]
                 HNO2 + H2O2 -> ONOOH + H2O [k=717000]
                 HNO2 + HNO3 -> 2NO2 + H2O [k=0.045]
                 2HNO2 -> NO2 + NO + H2O [k=1.34]
                 HNO2 + Eaqm -> HNO2m [k=4e9]
                 HNO2 + H -> H2NO2 [k=3.32364e8]
                 NO2m + Hp -> HNO2 [k=1e10]
                 NO2m + OH -> NO2 + OHm [k=6e9]
                 NO2m + H -> NO + OHm [k=1.45221e9]
                 NO2m + Om -> ONOOm [k=3.5e9]
                 NO2m + Eaqm -> NO22m [k=4.1e9]
                 2NO3 -> N2O6 [k=790000]
                 NO3 + HNO2 -> HNO3 + NO2 [k=2e8]
                 NO3 + H2O2 -> HNO3 + HO2 [k=7.1e6]
                 NO3 + NO2m -> NO3m + NO2 [k=1.2e9]
                 NO3 + HO2 -> HNO3 + O2 [k=3e9]
                 NO3 + OHm -> NO3m + OH [k=8.2e7]
                 NO3 + OH -> NO2 + HO2 [k=1e10]
                 NO3 -> HNO3 + OH [k=9372.08]
                 ONOOH -> HNO2 + H2O2 [k=15620.1]
                 ONOOH -> HNO3 [k=0.9]
                 ONOOH -> NO2 + OH [k=0.35]
                 HOONOO -> NO2 + HO2 [k=0.026]
                 HOONOO -> HNO2 + O2 [k=0.0007]
                 HOONOO + HNO2 -> 2HNO3 [k=12]
                 HOONOO -> OONOOm + Hp [k=71000]
                 OONOOm -> O2m + NO2 [k=1.05]
                 OONOOm -> O2 + NO2m [k=1.35]
                 OONOOm + Hp -> HOONOO [k=5e10]
                 NO22m -> NO + 2OHm [k=8.33074e7]
                 HNO2m -> NO + OHm [k=5000]
                 NO + OH -> NO2m + Hp [k=1e10]
                 2NO -> 2NO2 [k=5.9e6]
                 NO + HO2 -> ONOOH [k=3.2e9]
                 NO + O2m -> ONOOm [k=5e9]
                 ONOOm -> NO + O2m [k=0.02]
                 ONOOH + Hp -> HNO2 + H2O2 [k=6.3]
                 ONOOH -> HNO3 [k=4.3]
                 ONOOm + OH -> NO + O2 + OHm [k=4.8e9]
                 N2O3 -> NO2 + NO [k=84000]
                 ONOOm + N2O3 -> 2NO2 + NO2m [k=3.1e8]
                 N2O3 -> 2NO2m + 2Hp [k=2000]
                 ONOOm + Hp -> ONOOH [k=5e10]
                 NO2 + ONOOm -> NO2m + NO3 [k=24000]
                 Eaqm + N2O -> Om + N2 [k=9.6e9]
                 H + N2O -> N2 + OH [k=90000]
                 Eprem -> Eaqm [k=2.4e12]
                 2Eaqm -> H2 + 2OHm [k=6.20063e9]
                 Eaqm + H -> H2 + OHm [k=2.49241e10]
                 Eaqm + OH -> OHm [k=3.32388e10]
                 Eaqm + Om -> 2OHm [k=1.14037e12]
                 Eaqm + H2O2 -> OHm + OH [k=1.22207e10]
                 Eaqm + HO2 -> HO2m [k=1.18887e10]
                 Eaqm + O2 -> O2m [k=2.11404e10]
                 Eaqm + O2m -> OHm + HO2m [k=1.18887e10]
                 2H -> H2 [k=4.62071e9]
                 H + OH -> H2O [k=1.02751e10]
                 H + H2O2 -> OH + H2O [k=3.15782e7]
                 H + O2 -> HO2 [k=1.2016e10]
                 H + O2m -> HO2m [k=1.02561e10]
                 H + HO2 -> H2O2 [k=1.02561e10]
                 H -> H2 + OH [k=109341]
                 2OH -> H2O2 [k=4.53602e9]
                 OH + H2O2 -> HO2 + H2O [k=2.65484e7]
                 OH + HO2 -> O2 + H2O [k=8.4449e9]
                 OH + O2m -> OHm + O2 [k=1.01835e10]
                 OH + H2 -> H + H2O [k=3.68021e7]
                 HO2 + O2m -> O2 + HO2m [k=9.47446e7]
                 2HO2 -> H2O2 + O2 [k=731489]
                 2H2O2 -> O2 + 2H2O [k=4.97541e-7]
                 Om + H2 -> H + OHm [k=1.10162e8]
                 Eaqm + HO2m -> Om + OHm [k=3.15432e9]
                 Om + OH -> HO2m [k=7.21743e9]
                 OH + HO2m -> O2m + H2O [k=7.47208e9]
                 Om -> HO2 + OHm [k=2.59366e10]
                 Om + O2 -> O3m [k=3.43357e9]
                 O3m -> Om + O2 [k=1977.1]
                 2HO2m -> 2OHm + O2 [k=4.97541e-6]
                 Om -> OHm + OH [k=4.89431e9]
                 Hp + OHm -> H2O [k=1.00458e11]
                 H2O2 -> Hp + HO2m [k=0.0541918]
                 Hp + HO2m -> H2O2 [k=4.40826e10]
                 H2O2 + OHm -> HO2m + H2O [k=1.13592e10]
                 HO2m -> H2O2 + OHm [k=6.15792e7]
                 OH -> Hp + Om [k=0.0541918]
                 Hp + Om -> OH [k=4.40826e10]
                 OH + OHm -> Om + H2O [k=1.13592e10]
                 Om -> OH + OHm [k=6.15792e7]
                 HO2 -> Hp + O2m [k=633138]
                 Hp + O2m -> HO2 [k=4.40826e10]
                 HO2 + OHm -> O2m + H2O [k=1.13592e10]
                 O2m -> HO2 + OHm [k=5.27071]
                 H -> Hp + Eaqm [k=3.90449]
                 Hp + Eaqm -> H [k=2.07423e10]
                 H + OHm -> Eaqm + H2O [k=1.92709e7]
                 Eaqm -> H + OHm [k=682.257]'
  []
[]

[Postprocessors]
  [PuO22p]
    type = ElementAverageValue
    variable = PuO22p
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [PuO2p]
    type = ElementAverageValue
    variable = PuO2p
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [PuO23p]
    type = ElementAverageValue
    variable = PuO23p
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [PuIV]
    type = ElementAverageValue
    variable = PuIV
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [PuIII]
    type = ElementAverageValue
    variable = PuIII
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [H2O2]
    type = ElementAverageValue
    variable = H2O2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [H2]
    type = ElementAverageValue
    variable = H2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [HNO2]
    type = ElementAverageValue
    variable = HNO2
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
  line_search = none   # full Newton steps; bt line search diverges on this stiff network
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'
  automatic_scaling = true
  compute_scaling_once = false   # recompute scaling every step: dt grows ~10 orders, so
                                 # factors fixed at dt0 go stale and freeze the slow channels
  # With automatic_scaling the residual is normalized, so an absolute tolerance is
  # meaningful: here the 1/dt term dominates the Jacobian diagonal (scaling factors
  # ~= dt for every species), so the scaled residual is small and nl_abs_tol binds.
  # A relative-only criterion (nl_abs_tol ~ 0) is unreachable on this ~1e17-stiff
  # network and diverges (DIVERGED_MAX_IT); a scaled nl_abs_tol converges robustly.
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
  l_max_its = 100
  nl_max_its = 30
  end_time = 1.0
  dtmin = 1e-15
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-12        # resolve the fast HNO3 dissociation / radical chemistry
    growth_factor = 1.2
    cutback_factor = 0.5
    optimal_iterations = 8
  []
  dtmax = 1e-2
[]

[Outputs]
  csv = true
[]
