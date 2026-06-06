# KineticReactionsPhysics

!syntax description /Physics/KineticReactions/KineticReactionsPhysics

## Overview

This `Physics` assembles and solves a kinetic (mass-action) diffusion-reaction
network. For each species $s$ it forms a continuous Galerkin equation

\begin{equation}
\frac{\partial C_s}{\partial t} = \nabla \cdot \left( D_s \nabla C_s \right) + S_s(C),
\end{equation}

where the net production is summed over every reaction $j$,

\begin{equation}
S_s = \sum_j \nu_{s,j}\, k_j \prod_r C_r^{a_{r,j}},
\end{equation}

with $\nu_{s,j}$ the net stoichiometry of $s$ in reaction $j$ (products minus
reactants), $k_j$ the forward rate constant, and $a_{r,j}$ the stoichiometry of
reactant $r$. Because a reactant's stoichiometry appears both as the rate-law
exponent and as the net-stoichiometry multiplier, the factor-of-two convention
for like-species reactions (e.g. $2A \rightarrow P$ giving $dA/dt = -2k[A]^2$)
is handled automatically.

The reactions are specified as a multi-line list; each reaction carries its
forward rate constant in metadata, for example

```
reactions = '2e_aq -> H2 + 2OHm [k=5.5e9]
             e_aq + Hp -> H     [k=1.7e10]'
```

Species must be named to match the [!param](/Physics/KineticReactions/KineticReactionsPhysics/species)
parameter (no `+`/`-` charge notation). Products that are not listed as species
(such as the solvent) are treated as untracked sinks and receive no equation.

The equations are built entirely from existing framework objects:
[ADTimeDerivative.md], [ADMatDiffusion.md] (which obeys the mesh `coord_type`),
an [ADParsedMaterial.md] holding each species' net rate, applied with
[ADMatBodyForce.md], and an [ADGenericConstantMaterial.md] for the diffusivities.

## Example Input File Syntax

!listing modules/chemical_reactions/test/tests/kinetic_reactions/burns1984_spherical_spur.i block=Physics

!syntax parameters /Physics/KineticReactions/KineticReactionsPhysics

!syntax inputs /Physics/KineticReactions/KineticReactionsPhysics

!syntax children /Physics/KineticReactions/KineticReactionsPhysics
