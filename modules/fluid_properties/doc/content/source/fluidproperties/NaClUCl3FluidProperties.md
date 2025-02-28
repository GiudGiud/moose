# NaClUCl3FluidProperties

!syntax description /FluidProperties/NaClUCl3FluidProperties

These properties are based on experiments and compilations reported in the MCRE design documentation \citep{MCRE2022a,MCRE2022b}. The salt used is a NaCl–UCl₃ eutectic with a composition of 66.7 mol% NaCl and 33.3 mol% UCl₃. Its melting temperature is 523 °C (i.e. 796.15 K). Most properties depend only on temperature; the fluid is treated as incompressible. The fluid properties are summarized in Table \ref{tab:naclucl3}, which lists the equations and corresponding literature references.

!table id=tab:naclucl3 caption=Table of properties, equations, and references for the NaCl–UCl₃ salt fluid properties.
| Properties                             | Equations | Reference |
| :------------------------------------- | :-------- | :-------- |
| Melting point, \(T_{mo}\) (K)            | \(796.15\) | \(\citep{MCRE2022a}\) |
| Density, \(\rho\) (kg/m\(^3\))             | \(\displaystyle 4212.6 - 1.0686\,T\) | \(\citep{MCRE2022b}\) |
| Specific heat capacity, \(c_p\) (J/kg-K)   | \(\displaystyle 8.900439\times10^{3} - 1.377936\times10^{1}\,T + 6.400369\times10^{-3}\,T^2 - \frac{8.443758\times10^{8}}{T^2}\) | \(\citep{MCRE2022b}\) |
| Thermal Conductivity, \(k\) (W/m-K)      | \(\displaystyle 5.6820 - 8.7832\times10^{-3}\,T + 4.0967\times10^{-6}\,T^2 - \frac{5.7642\times10^{5}}{T^2}\) | \(\citep{MCRE2022b}\) |
| Dynamic Viscosity, \(\mu\) (Pa\(\cdot\)s)           | \(\displaystyle 1.505\times10^{-4}\,\exp\!\left(\frac{2.666\times10^{4}}{8.314\,T}\right)\) | \(\citep{MCRE2022b}\) |
| Specific Enthalpy, \(h\) (J/kg)            | \(\displaystyle h(T)=8.900439\times10^{3}(T-T_{mo})-\frac{1.377936\times10^{1}}{2}(T^2-T_{mo}^2)+\frac{6.400369\times10^{-3}}{3}(T^3-T_{mo}^3)-\left(-8.443758\times10^{8}\right)\left(\frac{1}{T}-\frac{1}{T_{mo}}\right)\) | \(\citep{MCRE2022b}\) |
| Internal Energy, \(e\) (J/kg)     | \(\displaystyle e=h-\frac{p}{\rho}\) | \(\citep{MCRE2022b}\) |

## Range of Validity

The properties defined in `NaClUCl3FluidProperties` are valid for the following conditions:
- **Temperature:**
  Density is valid for \(800\,\mathrm{K} \le T \le 1100\,\mathrm{K}\);
  Specific heat capacity and thermal conductivity are valid for \(800\,\mathrm{K} \le T \le 1000\,\mathrm{K}\);
  Dynamic viscosity is valid for \(800\,\mathrm{K} \le T \le 1100\,\mathrm{K}\).
- **Pressure:**
  The salt is assumed to be incompressible over pressures near atmospheric up to a few MPa.

## Uncertainties of NaCl–UCl₃ Salt Fluid Properties

Based on the MCRE documentation \citep{MCRE2022a,MCRE2022b}, the estimated uncertainties for the NaCl–UCl₃ salt properties are approximately:
\[
\begin{array}{|l|c|}
\hline
\textbf{Property} & \textbf{Uncertainty (\%)} \\
\hline
\text{Density} & \approx 1\% \\
\text{Specific Heat Capacity} & \approx 5\% \\
\text{Thermal Conductivity} & \approx 10\% \\
\text{Dynamic Viscosity} & \approx 15\% \\
\hline
\end{array}
\]

!syntax parameters /FluidProperties/NaClUCl3FluidProperties

!syntax inputs /FluidProperties/NaClUCl3FluidProperties

!syntax children /FluidProperties/NaClUCl3FluidProperties

!bibtex bibliography
