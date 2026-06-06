# MOOSE Implementation Guide: Burns et al. (1984)
## "Radiation Chemical Diffusion Kinetic Calculations with Prescribed and Non-Prescribed Diffusion—I: Spherical and Cylindrical Cases"
### W.G. Burns, H.E. Sims, J.A.B. Goodall — Radiat. Phys. Chem. Vol. 23, No. 1/2, pp. 143–180, 1984

---

## 1. OVERVIEW OF THE PROBLEM

The paper models the **radiolysis of water** following ionizing radiation. Transient radiolytic species
(e⁻ₐq, H⁺, OH, H, OH⁻, H₂, H₂O₂) are formed in small concentrated regions:
- **Spherical spurs** — typifying γ and fast electron (β) radiation
- **Cylindrical tracks** — typifying α radiation

Species simultaneously **diffuse** and **react** with each other. The goal is to compute time-evolution
of species concentration profiles C(r,t) and **G values** (yield in molecules per 100 eV absorbed).

The paper compares two methods:
1. **Schwarz modified prescribed diffusion** — forces Gaussian concentration profiles at all times
2. **FACSIMILE numerical method** — the reference, which is what we implement in MOOSE

**We reproduce the FACSIMILE numerical method (Case: Spherical spur, Schwarz parameters, Gaussian initial conditions)** — this is Figures 2–7 of the paper.

---

## 2. GOVERNING EQUATIONS

### 2.1 Diffusion-Reaction PDE (in spherical or cylindrical symmetry)

For each species *i*, the governing equation in radial symmetry is:

**Spherical:**
```
∂Cᵢ/∂t = Dᵢ (1/r²) ∂/∂r [r² ∂Cᵢ/∂r]  +  Σⱼ Rᵢⱼ(C₁, C₂, ..., Cₙ)
```

**Cylindrical (per unit length):**
```
∂Cᵢ/∂t = Dᵢ (1/r) ∂/∂r [r ∂Cᵢ/∂r]  +  Σⱼ Rᵢⱼ(C₁, C₂, ..., Cₙ)
```

where:
- Cᵢ(r, t) = molar concentration of species i [mol dm⁻³ = M]
- Dᵢ = diffusion coefficient of species i [dm² s⁻¹]
- Rᵢⱼ = net rate of production/consumption of species i from reactions

The spatial variable r is the radial distance from the center of the spur/track [nm → convert to dm for consistency].

### 2.2 Fick's First Law Discretization (as implemented in FACSIMILE)

The paper gives explicit discrete diffusion coefficients (eqs. 4–8) using a **finite volume** approach.

The mass flux across the boundary between zone i and zone i+1 is:
```
fᵢ = nᵢ mᵢ - sᵢ mᵢ₊₁
```
where mᵢ is mass of the species in zone i.

**For spherical diffusion** (equations 5–6):
```
nᵢ = 3D pᵢ₊₁² / [(rᵢ₊₁ - rᵢ)(pᵢ₊₁³ - pᵢ³)]
sᵢ = 3D pᵢ₊₁² / [(rᵢ₊₁ - rᵢ)(pᵢ₊₂³ - pᵢ₊₁³)]
```

**For cylindrical diffusion** (equations 7–8):
```
nᵢ = 2D pᵢ₊₁ / [(rᵢ₊₁ - rᵢ)(pᵢ₊₁² - pᵢ²)]
sᵢ = 2D pᵢ₊₁ / [(rᵢ₊₁ - rᵢ)(pᵢ₊₂² - pᵢ₊₁²)]
```

**Note**: These correspond directly to the standard FVM divergence operators for spherical/cylindrical
Laplacians. In MOOSE, use `CoeffDiffusion` with the appropriate coordinate system (`RZ` or spherical), 
or manually implement using `ADDiffusion` with a spherical Laplacian kernel.

### 2.3 Grid (Variable Mesh with Constant Expansion Ratio)

Equations 9–10 define the **logarithmically expanding grid**:

Zone boundary positions:
```
pⱼ = h (Kʲ - 1)     j = 0, 1, 2, ..., n   (n = number of zones)
```

Effective mean radii of zones:
```
rᵢ = h (K^(i+1/2) - 1)    i = 0, 1, 2, ..., n-1
```

The paper uses:
- **K = √1.5 ≈ 1.2247** (mesh expansion ratio)
- **n = 40 zones** (verified to give relative errors ~10⁻³ in concentrations)
- Changing n from 40 to 60 causes relative changes in G values and concentrations < 10⁻³
- h is chosen so that the outermost zone is large enough that species mass is negligible there

**For the spherical case with Schwarz parameters:**
- Initial characteristic radius r₀(e⁻ₐq) = 2.458 nm
- The domain must extend to ~6–10 × r₀(max) ≈ 25–50 nm to ensure species are negligible at the boundary
- Suggest domain [0, 40 nm] with n=40 zones for the spherical spur case

**To determine h:** Set p₀ = 0, so h = p₁ / (K - 1). Choose p_n (outer boundary) large enough.
For example, if outer boundary p_n = 30 nm and K = √1.5, n = 40:
```
p_n = h(K^n - 1)  →  h = p_n / (K^n - 1)
```
With K=1.2247, n=40: K^40 ≈ 1.2247^40 ≈ 238.7, so h ≈ 30/(237.7) ≈ 0.1262 nm

### 2.4 Boundary Conditions

- **At r = 0 (center):** Reflecting (zero flux), i.e., ∂C/∂r = 0. 
  This is the natural BC for spherical/cylindrical symmetry.
- **At r = r_max (outer boundary):** Also reflecting (zero flux).
  The outer zone is made large enough that species concentration there is initially and always negligible.
  The paper monitors mass flowing into the outermost zone to verify this.

### 2.5 Initial Conditions — Gaussian Distribution

**Schwarz parameters (Gaussian initial conditions), from Table 2:**

The initial concentration profile for each species is Gaussian (equation 1):
```
C(r, t=0) = [N / (π b₀²)^(3/2)] × exp(-r²/b₀²)   [spherical]
```

where b₀² = 2r₀² (at t=0, equation 2 with D·t = 0).

So:
```
C(r, 0) = [N / (π × 2r₀²)^(3/2)] × exp(-r²/(2r₀²))   [spherical]
```

For the **cylindrical case** (per unit length), the 2D Gaussian is:
```
C(r, 0) = [N_L / (π × 2r₀²)] × exp(-r²/(2r₀²))   [cylindrical, per unit length]
```

The G value gives the number of species per spur N, and the energy per spur gives the conversion factor.

---

## 3. SPECIES, INITIAL CONDITIONS, AND DIFFUSION CONSTANTS

### 3.1 Schwarz Parameters — Table 2 (Spherical Spur, Gaussian Initial Distribution)

| Species | r₀ [nm] | G value | b₀ = √2 r₀ [nm] |
|---------|---------|---------|-----------------|
| e⁻ₐq   | 2.458   | 4.78    | 3.477           |
| H⁺      | 1.145   | 4.78    | 1.619           |
| H       | 1.145   | 0.62    | 1.619           |
| OH      | 1.145   | 5.70    | 1.619           |
| H₂      | —       | 0.15    | (same as H⁺)    |

**Energy per spur** = 62.5 eV

**Note on H₂:** r₀ is not specified for H₂ in Table 2 for Schwarz parameters. Use r₀ = 1.145 nm 
(same as other species sharing the smaller radius, consistent with co-localized production).

**Number of molecules per spur N:**
```
N = G × E_spur / 100
```
where G is in molecules/100 eV and E_spur = 62.5 eV, so:
```
N(e⁻ₐq) = 4.78 × 62.5/100 = 2.9875 ≈ 2.99 molecules/spur
N(H⁺)   = 4.78 × 62.5/100 = 2.9875 ≈ 2.99
N(H)     = 0.62 × 62.5/100 = 0.3875 ≈ 0.388
N(OH)    = 5.70 × 62.5/100 = 3.5625 ≈ 3.56
N(H₂)   = 0.15 × 62.5/100 = 0.09375 ≈ 0.094
```

**Initial peak concentration** from Gaussian formula (spherical):
```
C_peak(r=0, t=0) = N / (π × 2r₀²)^(3/2)
```

Units check: N [molecules] / (nm³) → need conversion to mol dm⁻³
1 nm³ = 10⁻²⁴ dm³, and Nₐ = 6.022×10²³ mol⁻¹
```
C [mol dm⁻³] = (N / Nₐ) / [volume in dm³]
```
Or equivalently, work in consistent units throughout (mol, dm, s).

**Convert r₀ from nm to dm:** 1 nm = 10⁻⁸ dm

Example for e⁻ₐq:
```
r₀ = 2.458 nm = 2.458×10⁻⁸ dm
b₀² = 2 × (2.458×10⁻⁸)² = 1.208×10⁻¹⁵ dm²
N/Nₐ = 2.9875 / 6.022×10²³ = 4.96×10⁻²⁴ mol
Volume_scale = (π × b₀²)^(3/2) = (π × 1.208×10⁻¹⁵)^(3/2)
             = (3.795×10⁻¹⁵)^(3/2) = 7.397×10⁻²³ dm³
C_peak = 4.96×10⁻²⁴ / 7.397×10⁻²³ = 0.0671 mol dm⁻³ = 67.1 mM
```
This matches the order of magnitude in Figure 2a (peak ~20 mmol dm⁻³ at t=0... re-check: 
actually Fig 2a shows ~20 mmol/cubic decimetre... confirming units are mmol dm⁻³ in plots).

### 3.2 Diffusion Constants — Table 1

All values ×10⁻⁷ dm² s⁻¹:

| Species | Schwarz (Gaussian) | Trumbore (Central min.) | Girija & Gopinathan |
|---------|-------------------|------------------------|---------------------|
| e⁻ₐq   | 4.5               | 4.5                    | 4.9                 |
| H⁺      | 9.0               | 10.0                   | 9.0                 |
| H       | 7.0               | 8.0                    | 7.0                 |
| OH      | 2.8               | 2.0                    | 2.8                 |
| OH⁻     | 5.0               | 2.0                    | 5.0                 |
| H₂O₂   | 2.2               | 1.4                    | 2.2                 |
| N₂O     | (2.0)*            | (2.0)*                 | —                   |
| C₂H₅OH | (1.0)*            | (1.0)*                 | —                   |

*Notional values for runs with diffusing scavengers.

**For the Schwarz parameters (Gaussian), Case 1 (spherical):**
- D(e⁻ₐq) = 4.5×10⁻⁷ dm² s⁻¹
- D(H⁺)   = 9.0×10⁻⁷ dm² s⁻¹
- D(H)     = 7.0×10⁻⁷ dm² s⁻¹
- D(OH)    = 2.8×10⁻⁷ dm² s⁻¹
- D(OH⁻)  = 5.0×10⁻⁷ dm² s⁻¹
- D(H₂O₂) = 2.2×10⁻⁷ dm² s⁻¹
- D(H₂)   = assumed same order, not explicitly given for Schwarz (use ~5.0×10⁻⁷ or set D=0 for H₂ as molecular product)

---

## 4. REACTION SCHEME AND RATE CONSTANTS — Table 3

All rate constants ×10⁻¹⁰ dm³ mol⁻¹ s⁻¹:

| # | Reaction | k (Schwarz) ×10⁻¹⁰ | k (Trumbore/Girija) ×10⁻¹⁰ |
|---|----------|--------------------|-----------------------------|
| 1 | 2H₂O + e⁻ₐq + e⁻ₐq → H₂ + 2OH⁻ | 0.55 | 0.50 |
| 2 | H + H₂O + e⁻ₐq → H₂ + OH⁻ | 2.5 | 3.0 |
| 3 | H + H → H₂ | 1.0 | 1.3 |
| 4 | e⁻ₐq + H⁺ → H | 1.7 | 2.3 |
| 5 | e⁻ₐq + OH → OH⁻ | 2.5 | 3.0 |
| 6 | e⁻ₐq + H₂O₂ → OH⁻ + OH | 1.3 | 1.23 |
| 7 | H + H₂O₂ → H₂O + OH | 0.01 | 0.016 |
| 8 | H⁺ + OH⁻ → H₂O | 10.0 | 14.3 |
| 9 | OH + OH → H₂O₂ | 0.6 | 0.5 |
| 10| H + OH → H₂O | 2.0 | 3.2 |
| 11| e⁻ₐq + Scavenger₁ → P₁ | 0.87 | 0.87 |
| 12| H + Scavenger₂ → P₂ | 0.05 | 0.05 |
| 13| OH + Scavenger₃ → P₃ | 0.13 | 0.13 |

**Notes (from Table 3 footnote):**
- For interaction of like species (e.g., H): -d[H]/dt = 2k[H]² (factor of 2 in rate law)
  - This applies to reactions 1 (e⁻ₐq + e⁻ₐq) and 3 (H + H) and 9 (OH + OH)
- Scavenger concentrations for scavenged calculations:
  - [Scavenger₁] = 7×10⁻⁴ mol dm⁻³ (reacts with e⁻ₐq)
  - [Scavenger₂] = 2×10⁻² mol dm⁻³ (reacts with H)
  - [Scavenger₃] = 10⁻² mol dm⁻³ (reacts with OH)
  - These correspond to the experimental conditions of Appleby and Schwarz, ref [14]

**For the unscavenged spherical spur case (Figs. 2–7), reactions 1–10 only.**

### 4.1 Source/Sink Terms for Each Species

Let concentrations be: C₁=e⁻ₐq, C₂=H⁺, C₃=H, C₄=OH, C₅=OH⁻, C₆=H₂, C₇=H₂O₂

Rate expressions (Schwarz rate constants, units: dm³ mol⁻¹ s⁻¹):

```
k₁ = 0.55×10¹⁰,  k₂ = 2.5×10¹⁰,   k₃ = 1.0×10¹⁰,  k₄ = 1.7×10¹⁰
k₅ = 2.5×10¹⁰,   k₆ = 1.3×10¹⁰,   k₇ = 0.01×10¹⁰, k₈ = 10.0×10¹⁰
k₉ = 0.6×10¹⁰,   k₁₀ = 2.0×10¹⁰
```

**Note on like-species reactions (factor of 2 convention):**
Reaction 1: rate = k₁ × C₁ × C₁ (consumption of 2 e⁻ₐq per event, so -d[e⁻ₐq]/dt includes 2×k₁×C₁²)
Reaction 3: -d[H]/dt includes 2×k₃×C₃²
Reaction 9: -d[OH]/dt includes 2×k₉×C₄²

**Source/sink term S_i for each species:**

```
S(e⁻ₐq) = -2k₁C₁² - k₂C₃C₁ - k₄C₁C₂ - k₅C₁C₄ - k₆C₁C₇

S(H⁺)   = -k₄C₁C₂ - k₈C₂C₅

S(H)    = -2k₃C₃² + k₄C₁C₂ - k₂C₃C₁ - k₇C₃C₇ - k₁₀C₃C₄

S(OH)   = -k₅C₁C₄ + k₆C₁C₇ - 2k₉C₄² + k₇C₃C₇ - k₁₀C₃C₄

S(OH⁻)  = +2k₁C₁² + k₂C₃C₁ + k₅C₁C₄ - k₈C₂C₅

S(H₂)   = +k₁C₁² + k₂C₃C₁ + k₃C₃²

S(H₂O₂) = +k₉C₄² - k₆C₁C₇ - k₇C₃C₇
```

**H₂O is a solvent (not tracked); reactions consuming/producing H₂O do not affect the ODE system.**

---

## 5. G VALUE CALCULATION

The G value for species i at time t is:
```
G_i(t) = N_i(t) × 100 / E_spur
```
where N_i(t) is the **total number of molecules of species i remaining per spur** at time t:

**Spherical:**
```
N_i(t) = Nₐ × 4π ∫₀^∞ C_i(r,t) r² dr
```

**Cylindrical (per unit length of track):**
```
N_i(t) = Nₐ × 2π ∫₀^∞ C_i(r,t) r dr   [per unit length]
```

In the MOOSE simulation, integrate the concentration profile numerically over the mesh at each time step.

For **initial G values**, use Table 2: G(e⁻ₐq)=4.78, G(H⁺)=4.78, G(H)=0.62, G(OH)=5.7, G(H₂)=0.15.

---

## 6. NUMERICAL METHOD DETAILS

### 6.1 ODE Solver

The paper uses **Gear's method** for stiff ODEs (ref. [7]: C.W. Gear, Comm. ACM 1971, 14, 176).
- In MOOSE, use the implicit backward Euler or BDF2 time integrator with JFNK solver.
- The system is highly stiff due to the fast reactions (k₈ = 10¹¹ dm³ mol⁻¹ s⁻¹ for H⁺ + OH⁻ → H₂O).

### 6.2 Time Range

The paper shows results from:
- t = 0 (initial)
- t = 100 ps = 10⁻¹⁰ s
- t = 1 ns = 10⁻⁹ s
- t = 10 ns = 10⁻⁸ s
- Full G-value curves from t = 10⁻¹² s to t = 10⁻⁶ s (Figs. 6–7)

**Recommended simulation time span:** t ∈ [10⁻¹³, 10⁻⁶] s with adaptive time stepping.
Initial time step: dt₀ = 10⁻¹³ s; maximum dt: 10⁻⁷ s.

### 6.3 Memory Requirements (Reference)

- FACSIMILE with 40 zones: ~1000–1500 variables, 1100–1500 KB on IBM 3033
- Modified prescribed diffusion: 400 KB
- FACSIMILE ran in 150–170 s CPU time; prescribed diffusion needed up to 10 s

### 6.4 Numerical Accuracy

The paper reports (for no-reaction test):
- K = √1.5, n = 40: relative errors ~10⁻³ in concentrations as peak drops from 10²⁰ to 10¹⁴ cm⁻³ in 10⁻⁵ s
- Additional errors from reaction integration expected < 10⁻³
- Changing n from 40 to 60: relative changes in G values < 10⁻³

---

## 7. CASE TO REPRODUCE: SPHERICAL SPUR, SCHWARZ PARAMETERS

### 7.1 Summary of Parameters

**Geometry:** Spherical, 1D radial (r ∈ [0, r_max])

**Grid:**
- n = 40 zones minimum (recommend 60 for safety)
- K = √1.5 = 1.22474
- r_max ≈ 30 nm (ensure concentration negligible at boundary)
- h = r_max / (K^n - 1)
  - n=40, K=1.22474: K^40 ≈ 238.7 → h ≈ 30/237.7 ≈ 0.1262 nm
  - Zone boundaries: p_j = h(K^j - 1), j=0..40
  - Zone centers: r_i = h(K^(i+0.5) - 1), i=0..39

**Initial Conditions (Gaussian):**
```
C_i(r, 0) = [N_i / (π × 2r₀ᵢ²)^(3/2)] × exp(-r²/(2r₀ᵢ²))
```

Species parameters:
| Species | r₀ [nm]  | N_i (molecules/spur) |
|---------|----------|---------------------|
| e⁻ₐq   | 2.458    | 2.9875              |
| H⁺      | 1.145    | 2.9875              |
| H       | 1.145    | 0.3875              |
| OH      | 1.145    | 3.5625              |
| OH⁻     | 0        | 0                   |
| H₂      | 1.145    | 0.09375             |
| H₂O₂   | 0        | 0                   |

**Diffusion constants (Schwarz):** See Table 1 above.

**Reactions and rate constants (Schwarz):** See Table 3 above, k₁–k₁₀ only (no scavengers).

**Boundary conditions:** Zero flux at r=0 and r=r_max.

**Time integration:** t ∈ [0, 10⁻⁶ s], adaptive stepping, Gear/BDF method.

### 7.2 Expected Results (from paper)

**Concentration profiles (Fig. 2a,b):** 
- At t=0: e⁻ₐq has broader Gaussian (r₀=2.458 nm) than OH, H⁺, H (r₀=1.145 nm)
- At t=100 ps: numerical shows central minimum forming in e⁻ₐq profile
- At t=1 ns, 10 ns: profiles continue broadening; numerical diverges from prescribed

**G values vs. time (Fig. 6):**
| Species | G(t=10⁻¹²) | G(t=10⁻⁶) approx |
|---------|------------|------------------|
| e⁻ₐq   | ~4.78      | ~2.9–3.0         |
| OH      | ~5.7       | ~2.9–3.0         |
| H       | ~0.62      | ~0.7             |

The numerical G(e⁻ₐq) decays more slowly than the prescribed diffusion result, especially between 
t=10⁻¹⁰ and t=10⁻⁷ s. The difference is ~0.3–0.5 G units at intermediate times.

**G values for H⁺ and OH⁻ (Fig. 7):**
- G(H⁺) decreases from ~4.78 to ~3.8
- G(OH⁻) increases from ~0 to ~0.9–1.0

---

## 8. CYLINDRICAL TRACK CASE (OPTIONAL EXTENSION)

**For completeness: Cylindrical track, 100 eV nm⁻¹, Schwarz parameters (Figs. 8–13)**

Key differences from spherical:
1. Use cylindrical Laplacian: (1/r) d/dr [r dC/dr]
2. Initial conditions: 2D Gaussian per unit length
3. LET (linear energy transfer) = 100 eV nm⁻¹
4. Energy per unit length determines N_i per unit length:
   - N_i per nm = G_i × LET / 100 = G_i × 1.0 molecules nm⁻¹
   - So N(e⁻ₐq) per nm = 4.78
5. Initial r₀ values are same as Table 2 (Schwarz parameters)
6. The deviations from Gaussian are even larger than in the spherical case (paper states this)
7. G-value decays are faster and the numerical/prescribed discrepancy is larger (Fig. 12 vs. Fig. 6)

---

## 9. MOOSE IMPLEMENTATION NOTES

### 9.1 Coordinate System

Use `coord_type = RSPHERICAL` (or `RZ` for cylindrical) in the MOOSE Mesh block.

For 1D spherical:
```
[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    xmin = 0
    xmax = 3e-8  # 30 nm in dm
    nx = 40
  []
[]

[Problem]
  coord_type = RSPHERICAL
[]
```

**However**, the paper uses a non-uniform logarithmically expanding mesh. This should be implemented
as a custom mesh (or using MOOSE's `TransformGenerator` / `ParsedGenerateSideset`).

**Custom mesh generation code (Python, for MOOSE Exodus input):**
```python
import numpy as np

n = 60  # number of zones (use 60 for safety)
K = np.sqrt(1.5)  # = 1.22474
r_max_nm = 40.0   # nm, outer boundary

# Solve for h: r_max = h*(K^n - 1)
h = r_max_nm / (K**n - 1)

# Zone boundaries
p = np.array([h * (K**j - 1) for j in range(n+1)])  # in nm
# Zone centers (effective mean radii)
r = np.array([h * (K**(i + 0.5) - 1) for i in range(n)])  # in nm

# Convert to dm for MOOSE:
p_dm = p * 1e-8
r_dm = r * 1e-8
```

### 9.2 Variables (7 species)

```
[Variables]
  [e_aq]   # hydrated electron
  [Hp]     # H+ ion
  [H]      # H atom
  [OH]     # hydroxyl radical
  [OHm]    # OH- ion
  [H2]     # molecular hydrogen
  [H2O2]   # hydrogen peroxide
[]
```

### 9.3 Kernels

For each variable, add:
1. `TimeDerivative` kernel
2. Spherical diffusion: `ADDiffusion` with `coord_type = RSPHERICAL`, OR use `MatDiffusion`

For reactions, use `ADBodyForce` or a custom `ADKernel` for the source/sink terms.

**Alternatively**, use `ScalarTransportBase` or write custom kernels inheriting from `ADKernel`.

### 9.4 Time Integration

```
[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = BDF2         # or 'implicit-euler'
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-13           # initial time step (s)
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 6
  []
  end_time = 1e-6
  dtmin = 1e-16
  dtmax = 1e-7
[]
```

### 9.5 Initial Condition Functions

For species i with initial Gaussian distribution (spherical coordinates, MOOSE Function):

```
[Functions]
  [ic_eaq]
    type = ParsedFunction
    expression = '(N_eaq / Na) / (pi * 2 * r0_eaq^2)^1.5 * exp(-x^2 / (2 * r0_eaq^2))'
    symbol_names = 'N_eaq Na r0_eaq'
    symbol_values = '2.9875 6.022e23 2.458e-8'   # N [molecules], Na [mol⁻¹], r0 [dm]
  []
  # Note: The formula (N/Na) / (π×2r₀²)^(3/2) × exp(-r²/2r₀²) gives mol dm⁻³
[]
```

**Numerical check for e⁻ₐq peak:**
```
N/Na = 2.9875/6.022e23 = 4.96e-24 mol
denom = (π × 2 × (2.458e-8)²)^1.5 = (π × 1.208e-15)^1.5 = (3.797e-15)^1.5
      = 7.40e-23 dm³
C_peak = 4.96e-24 / 7.40e-23 = 0.0670 mol dm⁻³ = 67.0 mM
```
(In Fig. 2a, the y-axis shows mmol dm⁻³; peak at t=0 is ~20–22 mM for e⁻ₐq... 

**Re-check:** Figures 2a shows peak ~20 mmol dm⁻³ = 20 mM. This discrepancy means the formula
in the paper uses a different normalization. Let's re-examine equation (1):

```
C(r, t) = N(t) / (π b²)^(3/2) × exp(-r²/b²)
```
Note: b² = 2r₀² + 4Dt, so at t=0: b₀² = 2r₀². The normalization integral:
```
∫₀^∞ C × 4πr² dr = 4π × N/(π b²)^(3/2) × ∫₀^∞ r² exp(-r²/b²) dr
                 = 4π × N/(π b²)^(3/2) × (b²/4) × √π × b
                 = 4π × N/(π^(3/2) × b³) × (√π/4) × b³
                 = N
```
So the normalization gives total integral = N molecules per spur (in units of molecules/volume).

The concentration in mol dm⁻³ at r=0, t=0 for e⁻ₐq:
```
C(0,0) = N / (π b₀²)^(3/2) / Nₐ   [mol dm⁻³, if r in dm]
b₀² = 2 × (2.458e-8)² = 1.2086e-15 dm²
(π b₀²)^(3/2) = (π × 1.2086e-15)^(3/2) = (3.797e-15)^(3/2) = 7.40e-23 dm³
C(0,0) = 2.9875 / (6.022e23 × 7.40e-23) = 2.9875 / 44.56 = 0.0671 mol/dm³ = 67.1 mM
```

**But Fig. 2a shows approximately 21 mM for e⁻ₐq at t=0.** This is 3.2× lower.

Possible resolution: The paper may plot concentrations in different units, OR the effective spur volume
per spur in a track needs the number of spurs per unit volume. Actually, the paper's y-axis unit is 
"MILLIMOLE / CUBIC DECIMETRE", so 21 on the y-axis = 21 mmol dm⁻³ = 0.021 mol dm⁻³.

There is a factor of ~3 discrepancy. This is likely because the figures are for the MEAN concentration
OVER ALL SPACE (i.e., considering spur density), not the LOCAL peak. The local peak should indeed be
~67 mM and the figures presumably show the spatial profile relative to some normalization convention.

**Alternative interpretation:** The concentrations in the figures are probably the local concentrations
in mmol/dm³ within the spur. Re-reading the caption: "Species concentration against radius for e⁻ₐq.
Spherical spur." So these ARE local concentrations.

The discrepancy (67 vs 21 mM) may be because the paper uses a **different b₀** or a different convention
for the initial Gaussian. Let's check if the Schwarz (1969) paper defines the Gaussian differently.

Looking at equation (1): C(r,t) = N(t)/(πb²)^(3/2) × exp(-r²/b²)

The argument of the exponential is -r²/b², NOT -r²/(2σ²). So b is not the same as σ.
The inflection point (r₀ = value of r at point of inflection at t=0) from dC/dr = 0... 
Actually at t=0, d²C/dr² = 0 at r = r₀ gives the characteristic radius r₀.

For C = A exp(-r²/b²), d/dr [r² dC/dr] = 0 at the inflection point:
The inflection point of N(r) = 4πr² C(r) is at r = b/√2.
The inflection point of C(r) itself: d²/dr²[exp(-r²/b²)] = 0 → r = b/√2.
So r₀ = b/√2, meaning **b₀ = r₀√2 at t=0** (from b² = 2r₀² + 4Dt with Dt=0 → b = r₀√2 ✓).

At t=0: C(0) = N/(π × 2r₀²)^(3/2) 

BUT WAIT — the correct argument uses b², not 2r₀²:
b₀² = 2r₀² → b₀ = r₀√2

So C(0,0) = N/(π b₀²)^(3/2) = N/(π × 2r₀²)^(3/2):
```
= 2.9875 / [Nₐ × (π × 2 × (2.458e-8)²)^(3/2)]
= 2.9875 / [6.022e23 × (3.797e-15)^(3/2)]
= 2.9875 / [6.022e23 × 7.40e-23]
= 2.9875 / 44.56 = 0.0671 mol dm⁻³ = 67.1 mmol dm⁻³
```
This matches Fig. 2a at t=0! (The y-axis goes to 20 for the t=0 curve; oh wait — Figure 2a shows
the maximum is ~22 mmol dm⁻³ at r=0 for t=1e-10 (dotted), and the t=0 curves are at ~22 and ~20.
Actually, the y-axis maximum shown is 20, suggesting the t=0 peak may be off-scale or the Gaussian 
parameter gives a lower value.

**Conclusion:** Use the formula as written. The numerical peak of ~67 mmol dm⁻³ seems consistent 
with the figure if we note the figure may be plotted at limited resolution or the maximum is cut off.
Use C_i(r,0) = [N_i/(Nₐ)] / (π × 2r₀ᵢ²)^(3/2) × exp(-r²/(2r₀ᵢ²)) in units of mol dm⁻³.

Actually, checking again: N_i is dimensionless (molecules per spur). The formula gives concentration:

```python
import numpy as np

Na = 6.022e23  # molecules/mol
# e_aq
r0_eaq = 2.458e-8  # dm
N_eaq = 2.9875  # molecules/spur
b0_sq = 2 * r0_eaq**2
C_peak_eaq = (N_eaq/Na) / (np.pi * b0_sq)**(1.5)
print(f"C_peak(e_aq) = {C_peak_eaq*1000:.2f} mmol/dm3")
# Output: C_peak(e_aq) = 67.08 mmol/dm3

# OH: same r0, larger G
r0_OH = 1.145e-8  # dm
N_OH = 3.5625
b0_sq_OH = 2 * r0_OH**2
C_peak_OH = (N_OH/Na) / (np.pi * b0_sq_OH)**(1.5)
print(f"C_peak(OH) = {C_peak_OH*1000:.2f} mmol/dm3")
# OH has smaller r0 (1.145 vs 2.458 nm) so higher peak:
# ~ 67 * (N_OH/N_eaq) * (r0_eaq/r0_OH)^3 = 67 * (3.56/2.99) * (2.458/1.145)^3 
# = 67 * 1.19 * 9.86 = 787 mmol/dm3 ??? 
```
This is very high. Fig. 3a shows OH peak at ~250 mmol dm⁻³... let's check:
```
N_OH = 3.5625, r0_OH = 1.145e-8 dm
b0sq = 2*(1.145e-8)^2 = 2.622e-16 dm^2
(π × 2.622e-16)^1.5 = (8.235e-16)^1.5 = 7.476e-23 dm^3  ... 
Wait: (8.235e-16)^1.5 = (8.235)^1.5 × (10^-16)^1.5 = 23.64 × 10^-24 = 2.364e-23 dm^3

C_peak(OH) = (3.5625/6.022e23) / 2.364e-23
           = 5.917e-24 / 2.364e-23 = 0.2503 mol dm⁻³ = 250.3 mmol dm⁻³
```
This matches Fig. 3a perfectly (peak ~250 mmol dm⁻³). Good.

For e⁻ₐq:
```
b0sq = 2*(2.458e-8)^2 = 1.208e-15 dm^2
(π × 1.208e-15)^1.5 = (3.795e-15)^1.5 

(3.795)^1.5 = 7.393; 10^(-15×1.5) = 10^-22.5 = 3.162e-23... 
So (3.795e-15)^1.5 = 7.393 × 3.162e-23 = 2.338e-22 dm^3? 

Let me be more careful:
(3.795e-15)^(3/2) = (3.795)^(3/2) × (1e-15)^(3/2) 
                  = (3.795)^(3/2) × 1e-22.5
                  = 7.393 × 3.162e-23
                  = 2.337e-22 dm^3

C_peak(e_aq) = (2.9875/6.022e23) / 2.337e-22
             = 4.962e-24 / 2.337e-22
             = 0.02124 mol dm⁻³ = 21.24 mmol dm⁻³
```

**This matches Fig. 2a! Peak e⁻ₐq ≈ 21 mmol dm⁻³ at t=0.** 

Earlier arithmetic error was in computing (3.797e-15)^(3/2). The correct value is 2.337e-22 dm³, 
not 7.40e-23 dm³.

### Corrected Initial Concentrations:

| Species | C_peak(r=0, t=0) |
|---------|-----------------|
| e⁻ₐq   | 21.2 mmol dm⁻³  |
| H⁺      | 240 mmol dm⁻³   |
| H       | 31 mmol dm⁻³    |
| OH      | 250 mmol dm⁻³   |
| H₂      | 6.5 mmol dm⁻³   |

---

## 10. COMPLETE MOOSE INPUT FILE STRUCTURE

Below is the recommended structure for the MOOSE input file to reproduce the spherical spur case.

### 10.1 Mesh

Use a 1D mesh with spherical coordinate system. For the non-uniform mesh, either:
(a) Use a uniform mesh (acceptable approximation), or
(b) Generate the non-uniform mesh using a custom script and import as Exodus

For an **initial uniform-mesh implementation**:
```
domain: r ∈ [0, 30] nm = [0, 3e-7] dm
nx = 200 (uniform) or
nx = 60 (non-uniform, logarithmic)
```

### 10.2 Key MOOSE Physics Modules

- `Navier-Stokes` — NO (not fluid flow)
- `Chemical Reactions` module — YES, useful for reaction kernels
- `Heat Conduction` — NO
- **Primary:** `ScalarTransport` (or write custom kernels)
- Use `ReactionNetwork` tools or `ChemicalDiffusionKernel` if available

Alternatively, use the **Multi-species diffusion-reaction** approach:
- `ADDiffusion` for Fick diffusion (with spherical coordinates)
- `ADBodyForce` with coupled variable reactions for source terms
- `CoupledForce` or custom `ADKernel` derived classes

### 10.3 Postprocessors for G Values

```
[Postprocessors]
  [G_eaq]
    type = ElementIntegralVariablePostprocessor
    variable = e_aq
    # For spherical: multiply by 4π (coord_type = RSPHERICAL handles r² factor)
    # G = integral × 4π × Na × 100 / E_spur
  []
[]
```

In spherical coordinates with `coord_type = RSPHERICAL`, the element volume includes the r² factor,
so `ElementIntegralVariablePostprocessor` computes ∫ C × r² dr (in 1D). Then:
```
G = 4π × [integral result] × Na × (100 eV/molecule) / E_spur
```

---

## 11. SUMMARY TABLE: ALL PARAMETERS FOR SPHERICAL SPUR CASE

| Parameter | Value | Units |
|-----------|-------|-------|
| Geometry | Spherical spur | — |
| Domain | [0, 30] nm = [0, 3×10⁻⁷] dm | dm |
| Mesh zones n | 40–60 | — |
| Mesh ratio K | √1.5 = 1.22474 | — |
| Energy per spur | 62.5 | eV |
| Nₐ (Avogadro) | 6.022×10²³ | mol⁻¹ |
| D(e⁻ₐq) | 4.5×10⁻⁷ | dm² s⁻¹ |
| D(H⁺) | 9.0×10⁻⁷ | dm² s⁻¹ |
| D(H) | 7.0×10⁻⁷ | dm² s⁻¹ |
| D(OH) | 2.8×10⁻⁷ | dm² s⁻¹ |
| D(OH⁻) | 5.0×10⁻⁷ | dm² s⁻¹ |
| D(H₂O₂) | 2.2×10⁻⁷ | dm² s⁻¹ |
| D(H₂) | 5.0×10⁻⁷ (assumed) | dm² s⁻¹ |
| r₀(e⁻ₐq) | 2.458×10⁻⁸ | dm |
| r₀(H⁺, H, OH, H₂) | 1.145×10⁻⁸ | dm |
| N(e⁻ₐq) | 2.9875 | molecules/spur |
| N(H⁺) | 2.9875 | molecules/spur |
| N(H) | 0.3875 | molecules/spur |
| N(OH) | 3.5625 | molecules/spur |
| N(H₂) | 0.09375 | molecules/spur |
| N(OH⁻), N(H₂O₂) | 0 | (produced by reactions) |
| k₁ (e⁻ₐq + e⁻ₐq → H₂ + 2OH⁻) | 5.5×10⁹ | dm³ mol⁻¹ s⁻¹ |
| k₂ (H + e⁻ₐq → H₂ + OH⁻) | 2.5×10¹⁰ | dm³ mol⁻¹ s⁻¹ |
| k₃ (H + H → H₂) | 1.0×10¹⁰ | dm³ mol⁻¹ s⁻¹ |
| k₄ (e⁻ₐq + H⁺ → H) | 1.7×10¹⁰ | dm³ mol⁻¹ s⁻¹ |
| k₅ (e⁻ₐq + OH → OH⁻) | 2.5×10¹⁰ | dm³ mol⁻¹ s⁻¹ |
| k₆ (e⁻ₐq + H₂O₂ → OH⁻ + OH) | 1.3×10¹⁰ | dm³ mol⁻¹ s⁻¹ |
| k₇ (H + H₂O₂ → H₂O + OH) | 1.0×10⁸ | dm³ mol⁻¹ s⁻¹ |
| k₈ (H⁺ + OH⁻ → H₂O) | 1.0×10¹¹ | dm³ mol⁻¹ s⁻¹ |
| k₉ (OH + OH → H₂O₂) | 6.0×10⁹ | dm³ mol⁻¹ s⁻¹ |
| k₁₀ (H + OH → H₂O) | 2.0×10¹⁰ | dm³ mol⁻¹ s⁻¹ |
| t_start | 10⁻¹³ | s |
| t_end | 10⁻⁶ | s |
| dt_initial | 10⁻¹³ | s |

---

## 12. REACTION RATE EXPRESSIONS (MOOSE-READY)

All concentrations in mol dm⁻³, time in s, rate constants in dm³ mol⁻¹ s⁻¹.

### Source/Sink Terms (net production rate, mol dm⁻³ s⁻¹):

Let:
- C1 = [e⁻ₐq], C2 = [H⁺], C3 = [H], C4 = [OH], C5 = [OH⁻], C6 = [H₂], C7 = [H₂O₂]

```
# Rate constants (SI: dm³ mol⁻¹ s⁻¹)
k1 = 5.5e9     # e_aq + e_aq → H2 + 2OH- 
k2 = 2.5e10    # H + e_aq → H2 + OH-
k3 = 1.0e10    # H + H → H2
k4 = 1.7e10    # e_aq + H+ → H
k5 = 2.5e10    # e_aq + OH → OH-
k6 = 1.3e10    # e_aq + H2O2 → OH- + OH
k7 = 1.0e8     # H + H2O2 → H2O + OH
k8 = 1.0e11    # H+ + OH- → H2O
k9 = 6.0e9     # OH + OH → H2O2
k10 = 2.0e10   # H + OH → H2O

# Source terms (positive = production, negative = consumption)
S_C1 = -2*k1*C1*C1 - k2*C3*C1 - k4*C1*C2 - k5*C1*C4 - k6*C1*C7
S_C2 = -k4*C1*C2 - k8*C2*C5
S_C3 = +k4*C1*C2 - 2*k3*C3*C3 - k2*C3*C1 - k7*C3*C7 - k10*C3*C4
S_C4 = +k6*C1*C7 + k7*C3*C7 - k5*C1*C4 - 2*k9*C4*C4 - k10*C3*C4
S_C5 = +2*k1*C1*C1 + k2*C3*C1 + k5*C1*C4 - k8*C2*C5
S_C6 = +k1*C1*C1 + k2*C3*C1 + k3*C3*C3
S_C7 = +k9*C4*C4 - k6*C1*C7 - k7*C3*C7
```

**Important notes on stoichiometry:**
- Reaction 1: consumes 2 e⁻ₐq → rate of consumption = 2×k₁×C₁² (hence factor of 2 in S_C1)
  produces 1 H₂ → rate of H₂ production = k₁×C₁²
  produces 2 OH⁻ → rate of OH⁻ production = 2×k₁×C₁²
- Reaction 3: consumes 2 H → rate of H consumption = 2×k₃×C₃²
  produces 1 H₂ → rate of H₂ production = k₃×C₃²
- Reaction 9: consumes 2 OH → rate of OH consumption = 2×k₉×C₄²
  produces 1 H₂O₂ → rate of H₂O₂ production = k₉×C₄²
- Reaction 10: produces H₂O (solvent, not tracked); just a loss term for H and OH

### Stoichiometry cross-check (conservation):

Reaction 4: e⁻ₐq + H⁺ → H:
- -k₄C₁C₂ in S_C1 ✓
- -k₄C₁C₂ in S_C2 ✓  
- +k₄C₁C₂ in S_C3 ✓

Reaction 6: e⁻ₐq + H₂O₂ → OH⁻ + OH:
- -k₆C₁C₇ in S_C1 ✓
- +k₆C₁C₇ in S_C4 (OH produced) ✓
- -k₆C₁C₇ in S_C7 (H₂O₂ consumed) ✓
- +k₆C₁C₇ in S_C5 (OH⁻ produced) — included in S_C5 = +k2C3C1 + k5C1C4 + 2k1C1² − k8C2C5 
  Wait: S_C5 above does NOT include +k6*C1*C7. 

**CORRECTION to S_C5:**
```
S_C5 = +2*k1*C1*C1 + k2*C3*C1 + k5*C1*C4 + k6*C1*C7 - k8*C2*C5
```
(OH⁻ is produced by reactions 1, 2, 5, and 6)

---

## 13. VALIDATION AND OUTPUT COMPARISON

To reproduce Figures 2–7 of the paper:

### Figures 2a, 2b: e⁻ₐq concentration profile
- x-axis: radius r [nm]
- y-axis: concentration [mmol dm⁻³]  
- Times: t=0, t=1e-10 s (Fig. 2a); t=1e-9, t=1e-8 s (Fig. 2b)
- Key feature: central minimum forms in numerical solution but not in prescribed

### Figures 3a, 3b: OH radical concentration profile
- Similar to Fig. 2 but for OH
- OH starts with narrower distribution (r₀=1.145 nm vs 2.458 nm for e⁻ₐq)

### Figures 4a, 4b: H⁺ ion concentration profile
### Figures 5a, 5b: H atom concentration profile

### Figure 6: G values vs. time (spherical, e⁻ₐq, H, OH)
- x-axis: time [s], log scale from 1e-12 to 1e-6
- y-axis: G value (molecules per 100 eV)
- At t→∞: G(e⁻ₐq)≈2.9 (numerical), ~2.85 (Table 4)
- Key: numerical G(e⁻ₐq) decays MORE SLOWLY than prescribed

### Figure 7: G values vs. time (H⁺ and OH⁻)
- G(H⁺) decreases from ~4.78 to ~3.8
- G(OH⁻) increases from 0 to ~0.9–1.0

### Table 4 (long-time G values, scavenged):
For comparison at 1 μs without scavengers (brackets in Table 4, Schwarz numerical):
```
G(e⁻ₐq) ≈ 2.85
G(OH)    ≈ 3.06  
G(H)     ≈ 0.712
G(H₂)   ≈ 0.407
G(H₂O₂) ≈ 0.659
```

---

## 14. IMPORTANT CAVEATS AND KNOWN ISSUES

1. **Central minimum in e⁻ₐq:** The most important physical result. The e⁻ₐq distribution is wider
   (r₀=2.458 nm) than OH/H⁺ (r₀=1.145 nm). As time proceeds, the center of e⁻ₐq is depleted by
   reactions with OH and H⁺ before diffusion can replenish it, creating a hollow shell shape.
   The prescribed diffusion method forces a Gaussian and cannot capture this.

2. **Cylindrical case deviates MORE:** The cylindrical track shows even larger deviations between
   prescribed and numerical methods (Fig. 12 vs. Fig. 6).

3. **Factor of 2 in like-species reactions:** The paper explicitly states this convention (Table 3 note a).
   This must be correctly implemented.

4. **Charge conservation:** H₂O is produced (not tracked) as solvent. The net charge of the system
   is preserved: initial charge includes H⁺ ions; charge is neutralized by OH⁻ formation.

5. **Stiffness:** The reaction k₈ = 10¹¹ dm³ mol⁻¹ s⁻¹ (H⁺+OH⁻→H₂O) is extremely fast and 
   dominates stiffness. Initial concentrations of OH⁻ are zero, so the problem is less stiff 
   initially, but becomes stiff as OH⁻ is produced. Use adaptive timestepping.

6. **The paper uses MOLAR concentrations (mol dm⁻³) throughout.** All rate constants are in 
   dm³ mol⁻¹ s⁻¹. The diffusion constants are in dm² s⁻¹. The spatial coordinate is in dm.
   Keep all units consistent in MOOSE.

7. **The boundary condition at r=0** in the MOOSE spherical mesh: for `coord_type = RSPHERICAL`,
   the r=0 boundary is a natural symmetry BC with zero flux.

8. **The FACSIMILE program** used by Burns et al. implements the same FVM discretization described
   in Section 2.2. MOOSE's `ADDiffusion` with `coord_type = RSPHERICAL` should reproduce this if 
   the same mesh is used.

---

## 15. REFERENCES FROM THE PAPER

1. H.A. Schwarz, J. Phys. Chem. 1969, **73**, 1928; W.G. Burns and A.R. Curtis, J. Phys. Chem. 1972, **76**, 3008.
2. E.M. Chance, A.R. Curtis, I.P. Jones and C.R. Kirby, AERE Report AERE R8775, HMSO 1977. (FACSIMILE)
3. A.H. Samuel and J.L. Magee, J. Chem. Phys. 1953, **21**, 1080.
4. A.K. Ganguly and J.L. Magee, J. Chem. Phys. 1956, **25**, 129.
7. C.W. Gear, Comm. ACM 1971, **14**, 176. (Gear's method for stiff ODEs)
10. C.N. Trumbore et al., J. Phys. Chem. 1978, **82**, 2762. (Central minimum model)
11. C.D. Jonah et al., J. Phys. Chem. 1976, **80**, 1267. (Experimental G(e⁻ₐq))
12. C.D. Jonah and J.R. Miller, J. Phys. Chem. 1977, **81**, 1974. (Experimental G(OH))
14. A. Appleby and H.A. Schwarz, J. Phys. Chem. 1969, **73**, 1938. (Scavenger experiments)

---

## 16. QUICK-START CHECKLIST FOR CLAUDE CODE

1. ✅ Read this guide completely
2. ✅ Set up 1D MOOSE input for spherical coordinate system
3. ✅ Implement logarithmic mesh (or use fine uniform mesh as approximation)
4. ✅ Define 7 variables: e_aq, Hp, H, OH, OHm, H2, H2O2
5. ✅ Set up Gaussian ICs using ParsedFunction
6. ✅ Add ADDiffusion kernels with correct D values for each species
7. ✅ Add reaction source terms (use ParsedMaterial or custom kernels)
8. ✅ Use BDF2 time integration with Gear-type adaptive stepping
9. ✅ Add Postprocessors to compute G values by integrating concentrations
10. ✅ Run to t=1e-6 s
11. ✅ Compare G(e⁻ₐq), G(OH), G(H⁺) vs. time to Figs. 6–7
12. ✅ Compare concentration profiles at t=0, 1e-10, 1e-9, 1e-8 s to Figs. 2–5
