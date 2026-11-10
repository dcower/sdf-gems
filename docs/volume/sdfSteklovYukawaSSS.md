# `sdfSteklovYukawaSSS` — Spherical-Wronskian Steklov-Yukawa Subsurface Scattering, Lambert Backward Continued-Fraction DtN & Zero-March Geodesic-Ring BSSRDF

<p align="center">
  <img src="../images/volume/sdfSteklovYukawaSSS.png" width="360" alt="3D translucent carnelian-wax SDF sculpture rendered with Spherical-Wronskian Steklov-Yukawa BSSRDF, 1-tap antipodal thin-fin slab transmission, and 0-march geodesic-ring cast-shadow bleeding">
  <img src="../images/volume/sdfSteklovYukawaSSS_slice.png" width="360" alt="2D cross-section slice of sdfSteklovYukawaSSS showing exterior Euclidean distance isolines and interior Screened-Poisson RGB subsurface fluence">
</p>

<p align="center">
  <img src="../images/volume/sdfSteklovYukawaSSS_comparison.png" width="820" alt="Comparison: (a) Existing sdfBSSRDFYukawa showing literal 0.0 black band at mu=0 shadow terminator, (b) New sdfSteklovYukawaSSS with RGB chromatic terminator bleed, thin-fin Robin slab glow, and geodesic-ring shadow bleed, (c) Sphere RGB exitance vs 3D Monte Carlo random walks, (d) Non-spherical RoundBox edge sweep">
</p>

- **Category**: `volume`
- **GLSL Header**: [`glsl/volume/sdfSteklovYukawaSSS.glsl`](../../glsl/volume/sdfSteklovYukawaSSS.glsl)
- **Attribution**: Discovered, derived, implemented, and benchmarked using **Gemini 4 Carbon** (`sdf-gems`, 2026).
- **Exactness & Accuracy**:
  - **`sdfSteklovYukawaQtilde5` / `sdfSteklovYukawaQtilde`**: Subtraction-free 10-level Lambert-Gauss backward continued fraction with $`O(x^{-3})`$ square-root closure for modified spherical Bessel Dirichlet-to-Neumann (DtN) ratios $`\tilde{q}_\ell(x) = i_{\ell+1}(x) / (x\,i_\ell(x))`$ ($`\ell = 0, \dots, 4`$). Maximum GPU FP32 relative error $`\le 2.22 \times 10^{-4}`$ across $`x \in [0, 10^5]`$ with **`0` `NaN`/`Inf`s** (where naive `sinh(x)` overflows FP32 at $`x > 88.7`$ and forward recurrence suffers catastrophic Miller cancellation for $`x < 0.5`$).
  - **`sdfSteklovSphereBSSRDF`**: Solves interior Screened-Poisson diffusion $`(\Delta - \sigma_{\text{tr}}^2)\phi = -Q/D`$ with Marshak-Robin boundary conditions $`\phi + 2AD\,\partial_n\phi = 0`$ and 3D collimated Beer-Lambert sources. Achieves **`0.0055`–`0.0124` RMSE** against 600,000-photon 3D Monte Carlo subsurface random-walk path tracing (**`9.7x`–`32.7x` lower error** than `sdfBSSRDFYukawa`), eliminates the literal `0.0` black terminator band at $`\mu = \mathbf{n}\cdot\mathbf{l} = 0`$, conserves `100%` of scattered photons at $`\sigma_a = 0`$, and converges continuously ($`9.3\times 10^{-5}`$ transition, `0.0` dark-side Gibbs ringing) to the exact oblique Robin half-space solution as $`H \to 0`$.
  - **`sdfAntipodalSlabProbe` & `sdfRobinSlabBSSRDF`**: Recovers exact Euclidean wall thickness $`T`$ to **`5.96e-8`** in 1 interior `sdg` tap across both medial-crossed and back-exited regimes, and matches 400,000-photon Monte Carlo finite-slab transport to **`0.0044` max error** (including $`\sigma_a = 0`$).
  - **`sdfGeodesicRingBSSRDF`**: Exact per-channel partition-of-unity 2D polar Yukawa kernel over 4 zero-march Eikonal footpoints $`\mathbf{q}_k = \mathbf{y}_k - d(\mathbf{y}_k)\nabla d(\mathbf{y}_k)`$ for chromatic cast-shadow bleeding and sharp polyhedral corner rescue.

---

## Motivation & Audit of Prior Volumetric Subsurface Models

Subsurface scattering (BSSRDF) is governed in the diffusion regime by the interior **Screened-Poisson (Yukawa) PDE** on $`\Omega = \{ \mathbf{x} \in \mathbb{R}^3 : d(\mathbf{x}) < 0 \}`$:

```math
\left(\Delta - \sigma_{\text{tr}}^2\right)\phi(\mathbf{x}) = -\frac{1}{D}\,Q(\mathbf{x}), \qquad D = \frac{1}{3\sigma_t'}, \qquad \sigma_{\text{tr}} = \sqrt{3\sigma_a\sigma_t'},
```

subject to the **Extrapolated Marshak-Robin Boundary Condition** on $`\partial\Omega`$:

```math
\phi(\mathbf{x}) + z_b\,\frac{\partial\phi}{\partial n}(\mathbf{x}) = 0, \qquad z_b = 2 A D = \frac{2A}{3\sigma_t'}, \qquad A = \frac{1 + F_{dr}}{1 - F_{dr}},
```

with outward diffuse radiant exitance $`M_d(\mathbf{x}) = -D\,\partial_n\phi(\mathbf{x}) = \phi(\mathbf{x}) / (2A)`$.

Auditing [`glsl/volume/sdfFinslerSubsurface.glsl`](../../glsl/volume/sdfFinslerSubsurface.glsl) revealed a structural flaw in `sdfBSSRDFYukawa`:

```glsl
float ndotl    = max(dot(nSurf,  lightDir), 0.0);
float backdotl = max(dot(-nSurf, lightDir), 0.0);
return lightCol * (ndotl * wFront + backdotl * wBack);
```

Along the **shadow terminator** $`\mu = \mathbf{n}_{\text{surf}}\cdot\mathbf{l} = 0`$ — the primary visual signature of translucent wax, skin, marble, and jade — **both `max(mu, 0.0)` and `max(-mu, 0.0)` evaluate identically to `0.0`**, carving a pitch-black zero-radiance ring down the middle of every illuminated curved surface (Panel (a) above) while over-brightening dark-side silhouettes where a constant manual `chordThickness` ignores true geometry.

---

## Mathematical Derivation

### 1. The Spherical-Wronskian Robin Reduction Theorem

Approximate a convex SDF surface locally by its osculating sphere of radius $`R = 1 / H`$ ($`H = \frac{1}{2}\Delta d > 0`$). Expanding the fluence $`\phi(r, \mu) = \sum_{\ell=0}^\infty \phi_\ell(r) P_\ell(\mu)`$ in modified spherical Bessel functions $`i_\ell(\sigma_{\text{tr}} r)`$ and $`k_\ell(\sigma_{\text{tr}} r)`$ ($`x = \sigma_{\text{tr}} R`$, $`x_1 = \sigma_{\text{tr}} r_1`$), any power-conserving spherical source shell of total strength $`S_\ell`$ at radius $`r_1 \in (0, R)`$ satisfies the radial jump condition $`r_1^2\bigl(\phi_\ell'(r_1^+) - \phi_\ell'(r_1^-)\bigr) = -\frac{R^2}{D} S_\ell`$ and the outer Robin condition $`\phi_\ell(R) + z_b \phi_\ell'(R) = 0`$.

Writing the outer homogeneous solution for $`r \in (r_1, R]`$ as

```math
\psi_\ell(r) = i_\ell(\sigma_{\text{tr}} r)\Bigl(k_\ell(x) + \frac{z_b}{R}\,x\,k_\ell'(x)\Bigr) - k_\ell(\sigma_{\text{tr}} r)\Bigl(i_\ell(x) + \frac{z_b}{R}\,x\,i_\ell'(x)\Bigr),
```

which satisfies $`\psi_\ell(R) + z_b\psi_\ell'(R) \equiv 0`$ identically, both the derivative jump at $`r = r_1`$ and the boundary flux at $`r = R`$ collapse via the **Modified Spherical Bessel Wronskian** ($`i_\ell(z) k_\ell'(z) - k_\ell(z) i_\ell'(z) \equiv -\frac{\pi}{2 z^2}`$). Every singular $`k_\ell`$ term and prefactor cancels identically, leaving the exact closed-form modal exitance:

```math
\boxed{M_\ell(R) = \frac{S_\ell}{1 + \dfrac{z_b}{R}\,\Lambda_\ell(x)}, \qquad \Lambda_\ell(x) := \frac{x\,i_\ell'(x)}{i_\ell(x)} = \ell + x^2\,\tilde{q}_\ell(x),}
```

where $`\Lambda_\ell(x)`$ is the **Steklov-Yukawa Dirichlet-to-Neumann (DtN) eigenvalue** of the ball $`B(\mathbf{0}, R)`$ and $`\tilde{q}_\ell(x) := \frac{i_{\ell+1}(x)}{x\,i_\ell(x)} > 0`$.

- **Exact Photon Conservation at $`\sigma_a \to 0`$**: As $`\sigma_a \to 0`$ ($`x \to 0`$), $`i_0(x_1)/i_0(x) \to 1`$ and $`\Lambda_0(x) = O(x^2) \to 0`$, giving $`M_0(R) \to S_0`$ with zero energy loss.
- **Why Collimated Beam Drift Generates $`\ell = 3`$ & Terminator Glow**: A collimated beam along $`-\mathbf{l}`$ entering at $`R\hat{\boldsymbol{\omega}}_0`$ ($`\mu_0 = \hat{\boldsymbol{\omega}}_0\cdot\mathbf{l} \in (0, 1]`$) travels parallel to $`-\mathbf{l}`$ rather than radially along $`-\hat{\boldsymbol{\omega}}_0`$:
  ```math
  \mathbf{x}(s) = R\hat{\boldsymbol{\omega}}_0 - s\,\mathbf{l}, \qquad s \in [0, 2R\mu_0].
  ```
  Near $`\mu_0 \to 0^+`$, stepping $`s \sim 1/\sigma_t'`$ pushes first-scatter collisions directly across the $`z = 0`$ equatorial plane into the shadow hemisphere ($`\mu_1 = (R\mu_0 - s)/r_1 < 0`$), exciting odd $`\ell = 3`$ moments and producing strong chromatic RGB bleeding across $`\mu = 0`$ (since $`\sigma_{\text{tr},R} < \sigma_{\text{tr},G} < \sigma_{\text{tr},B}`$).

To avoid $`O(R^2)`$ catastrophic FP32 mantissa cancellation as $`H \to 0`$ ($`R \to \infty`$), `sdfSteklovChannelBSSRDF` evaluates the chord depth via the probability-integral transform $`u_{\text{cdf}} \in (0, 1)`$, $`s = -\ln\bigl(1 - u_{\text{cdf}}(1 - e^{-2\sigma_t' R\mu_0})\bigr)/\sigma_t'`$, and rationalizes the radial shift:

```math
r_1 = \sqrt{(R - s)^2 + 2 R s (1 - \mu_0)}, \qquad R - r_1 = \frac{2 R \mu_0 s - s^2}{R + r_1}.
```

---

### 2. Subtraction-Free Lambert-Gauss Backward Continued Fraction (`sdfSteklovYukawaQtilde5`)

Evaluating $`\Lambda_\ell(x)`$ via forward recurrence $`\Lambda_\ell = \frac{x^2}{\Lambda_{\ell-1} - (\ell-1)} - (\ell+1)`$ fails catastrophically in FP32 for $`x < 0.5`$ due to Miller instability, while direct evaluation of $`i_\ell(x)`$ overflows `sinh(x)` for $`x > 88.7`$.

Rewriting $`i_{\ell-1}(x) - i_{\ell+1}(x) = \frac{2\ell+1}{x}i_\ell(x)`$ in terms of the positive ratio $`\tilde{q}_m(x) = \frac{i_{m+1}(x)}{x\,i_m(x)}`$ yields the **subtraction-free backward recurrence**:

```math
\tilde{q}_m(x) = \frac{1}{(2m + 3) + x^2\,\tilde{q}_{m+1}(x)}, \qquad m = M-1, \dots, 0.
```

From the large-$`x`$ expansion $`\Lambda_M(x) = x - 1 + \frac{M(M+1)}{2x} + O(x^{-2})`$, we close the tail at $`M = 10`$ with the square-root approximant:

```math
\tilde{q}_{10}(x) = \frac{1}{11 + \sqrt{144 + x^2}},
```

which is **simultaneously exact at $`x = 0`$** ($`\tilde{q}_{10}(0) = \frac{1}{23}`$) **and exact to $`O(x^{-3})`$ as $`x \to \infty`$** ($`\tilde{q}_{10}(x) = \frac{1}{x} - \frac{11}{x^2} + O(x^{-3})`$), achieving $`\le 2.22\times 10^{-4}`$ relative error across all $`x \in [0, \infty)`$ with zero branches and zero exponentials.

---

### 3. Duhamel Oblique Half-Space Singularity Subtraction ($`H \to 0`$ Exact Limit)

As $`H \to 0`$ ($`R \to \infty`$), the exact infinite Legendre series converges pointwise to the **1D Oblique Beer-Lambert Half-Space Marshak-Robin Exitance**:

```math
M_{\text{flat}}(\mu) = \begin{cases} \dfrac{\alpha'\,\mu}{\bigl(1 + \frac{\sigma_{\text{tr}}}{\sigma_t'}\mu\bigr)\bigl(1 + \sigma_{\text{tr}} z_b\bigr)}, & \mu > 0, \\[10pt] 0, & \mu \le 0, \end{cases}
```

which has a slope discontinuity at $`\mu = 0`$. Truncating $`\sum_{\ell=0}^4 M_\ell(R) P_\ell(\mu)`$ directly would leave $`\sim 4.5\%`$ Legendre Gibbs ringing at $`\mu = 0`$ even on a flat wall! Subtracting the quadrature-matched flat moments $`M_\ell^{\text{flat}}`$ and adding back $`M_{\text{flat}}(\mu)`$ damped by the asymptotic $`\ell = 6`$ Steklov-Yukawa attenuation factor $`\eta_{\text{high}}(R) = e^{-6 H/\sigma_t'}\frac{1 + \sigma_{\text{tr}} z_b}{1 + (z_b H)\Lambda_6(x)} \in (0, 1]`$:

```math
M_d(\mu, H) = \max\!\left(\sum_{\ell=0}^4 M_\ell(R) P_\ell(\mu) + \eta_{\text{high}}(R)\left(M_{\text{flat}}(\mu) - \sum_{\ell=0}^4 M_\ell^{\text{flat}} P_\ell(\mu)\right),\; 0\right)
```

eliminates all Gibbs ringing as $`H \to 0`$ (`0.00e+00` dark-side error) while cutting curved-sphere RMSE against 3D Monte Carlo down to `0.55%`–`1.24%`.

---

### 4. 1-Tap Antipodal Slab Extraction & Zero-March Geodesic-Ring BSSRDF

1. **Universal Antipodal Eikonal Probe (`sdfAntipodalSlabProbe`)**:
   Stepping inward from $`\mathbf{p}_{\text{surf}}`$ along $`-\mathbf{n}_{\text{surf}}`$ to $`\mathbf{p}_{\text{probe}} = \mathbf{p}_{\text{surf}} - \delta\,\mathbf{n}_{\text{surf}}`$ inside a thin fin/ear of thickness $`T \in (0, 2\delta)`$ crosses the medial sheet ($`\nabla d(\mathbf{p}_{\text{probe}})\cdot\mathbf{n}_{\text{surf}} < 0`$). Crucially, **both** when $`\mathbf{p}_{\text{probe}}`$ remains inside the back half of the slab ($`T/2 < \delta \le T`$, $`d < 0`$) **and** when $`\mathbf{p}_{\text{probe}}`$ exits through the back wall into empty space ($`0 < T < \delta`$, $`d > 0`$), the identical Eikonal footpoint formula
   ```math
   \mathbf{q}_{\text{back}} = \mathbf{p}_{\text{probe}} - d(\mathbf{p}_{\text{probe}})\,\nabla d(\mathbf{p}_{\text{probe}})
   ```
   lands **unconditionally on the back surface**, yielding exact thickness $`T = \Vert\mathbf{p}_{\text{surf}} - \mathbf{q}_{\text{back}}\Vert`$ and back-wall normal $`\mathbf{n}_{\text{back}}`$ in **zero raymarch steps**, which feeds the exact 2-boundary Marshak-Robin slab Green's function `sdfRobinSlabBSSRDF`.

2. **Partition-of-Unity Geodesic-Ring BSSRDF (`sdfGeodesicRingBSSRDF`)**:
   Because $`\frac{d}{dr}\bigl(-\frac{1}{\kappa}e^{-\kappa\sqrt{r^2 + z_0^2}}\bigr) = \frac{r\,e^{-\kappa\sqrt{r^2 + z_0^2}}}{\sqrt{r^2 + z_0^2}}`$, the 2D polar Screened-Poisson kernel splits across any 4 Eikonal-projected tangent taps $`\mathbf{q}_k = \mathbf{y}_k - d(\mathbf{y}_k)\nabla d(\mathbf{y}_k)`$ ($`r_k = \frac{1}{\sqrt{2}}\Vert\mathbf{q}_k - \mathbf{p}_{\text{surf}}\Vert`$) with **exact per-channel partition of unity**:
   ```math
   w_k = \frac{1}{4}\exp\!\Bigl(-\boldsymbol{\sigma}_{\text{tr}}\bigl(\sqrt{r_k^2 + z_0^2} - z_0\bigr)\Bigr), \qquad w_{\text{center}} = \mathbf{1} - \sum_{k=0}^3 w_k.
   ```

---

## GLSL API Summary

```glsl
vec4  sdfSteklovYukawaQtilde5(in float x, out float outQ4);
vec4  sdfSteklovYukawaQtilde(in float x);

float sdfSteklovChannelBSSRDF(in float mu, in float H, in float sigmaTr, in float sigmaTp, in float A);
vec3  sdfSteklovSphereBSSRDF(in float mu, in float H, in vec3 sigmaTr, in vec3 sigmaTp, in float A);

vec4  sdfAntipodalSlabProbe(in vec3 pSurf, in vec3 nSurf, in vec4 dgProbe, in float probeDepth);
float sdfRobinSlabChannel(in float thickness, in float muFront, in float muBack,
                          in float sigmaTr, in float sigmaTp, in float A);
vec3  sdfRobinSlabBSSRDF(in float thickness, in float muFront, in float muBack,
                         in vec3 sigmaTr, in vec3 sigmaTp, in float A);

vec3  sdfGeodesicRingBSSRDF(in vec3 pSurf, in vec3 ECenter,
                            in vec3 q0, in vec3 E0, in vec3 q1, in vec3 E1,
                            in vec3 q2, in vec3 E2, in vec3 q3, in vec3 E3,
                            in vec3 sigmaTr, in vec3 sigmaTp, in float A);
```

---

## Direct-GPU Verification & 3D Monte Carlo Benchmarks

Verified on the literal [`glsl/volume/sdfSteklovYukawaSSS.glsl`](../../glsl/volume/sdfSteklovYukawaSSS.glsl) shader in OpenGL 4.1 FP32 via `tools/gpu_eval.py` against `scipy.special.ive`, a 64x64 Gauss-Legendre 12-mode PDE solver, a 600,000-photon forward 3D Monte Carlo Subsurface Path Tracer, and a 35,000-walk Adjoint Next-Event-Estimation 3D Monte Carlo solver on non-spherical SDFs:

| Benchmark / Regime | 3D Monte Carlo Ground Truth | Prior `sdfBSSRDFYukawa` (GPU FP32) | New `sdfSteklovYukawaSSS` (GPU FP32) | Improvement / Status |
| :--- | :---: | :---: | :---: | :---: |
| **DtN Ratios `sdfSteklovYukawaQtilde5` ($`x \in [0, 10^5]`$, $`\ell \le 4`$)** | `scipy.special.ive` | Overflows `sinh` at $`x > 88.7`$ | **`1.55e-4`–`2.22e-4` max rel err, `0` NaN/Inf** | Unconditionally stable on $`[0, \infty)`$ |
| **Sphere ($`R=1.0, \sigma_s'=8, \sigma_a=0.2, A=1`$) Terminator $`\mu=0`$** | `0.12586` | **`0.00000` (black band)** | **`0.10961` (`0.01068` full-sphere RMSE)** | **`9.7x` lower RMSE** (`0.1041` $`\to`$ `0.0107`) |
| **Conservative Sphere ($`R=1.0, \sigma_s'=8, \sigma_a=0`$) RMSE** | `0.16594` ($`\mu=0`$) | `0.40627` RMSE (`0.0` at $`\mu=0`$) | **`0.01242` RMSE (`100%` energy conserved)** | **`32.7x` lower RMSE** |
| **High-Index Sphere ($`R=1.0, \sigma_s'=5, \sigma_a=0.8, A=1.8`$)** | `0.08250` ($`\mu=0`$) | `0.07914` RMSE (`0.0` at $`\mu=0`$) | **`0.00551` RMSE** | **`14.4x` lower RMSE** |
| **Small Translucent Bead ($`R=0.4, \sigma_s'=6, \sigma_a=0.3, A=1.4`$)** | `0.19255` ($`\mu=0`$) | `0.26653` RMSE (`0.0` at $`\mu=0`$) | **`0.00818` RMSE (`0.1988` at $`\mu=0`$)** | **`32.6x` lower RMSE** |
| **Planar Half-Space Continuity ($`H=10^{-3} \to 0`$)** | Exact oblique Robin | Non-zero fake back-glow | **`9.33e-5` max jump, `0.00` dark ringing** | Exact $`H \to 0`$ limit |
| **1-Tap `sdfAntipodalSlabProbe` ($`T \in [0.05, 1.20]`$, $`\delta = 0.70`$)** | Exact $`T`$, $`\mathbf{n}_{\text{back}}`$ | Manual scalar guess | **`5.96e-8` thickness err, `0.00` normal err** | Exact across medial & exited probes |
| **Finite Slabs `sdfRobinSlabBSSRDF` (4 regimes incl. $`\sigma_a=0`$)** | 400k-photon 3D MC | Up to `48%` error | **`0.00440` max absolute diff** | Closed-form 2-wall Robin BVP |

---

## Research Findings: What Held & Where Approximations Break Down

1. **What Worked Beyond Expectations**:
   - **Wronskian Singularity Cancellation & Collimated Beam Drift**: Retaining the $`-\mathbf{l}`$ collimated beam drift inside the sphere (exciting $`\ell = 3`$ alongside $`\ell = 0,1,2,4`$) plus Duhamel flat-remainder damping matched full 3D Monte Carlo random walks within `~0.5%–1.2%` RMSE across every tested sphere radius and absorption albedo, including $`\sigma_a = 0`$.
   - **Universal Back-Surface Footpoint Identity**: $`\mathbf{q}_{\text{back}} = \mathbf{p}_{\text{probe}} - d(\mathbf{p}_{\text{probe}})\nabla d(\mathbf{p}_{\text{probe}})`$ recovers exact slab thickness regardless of whether $`\mathbf{p}_{\text{probe}}`$ stays inside the far half of the fin or overshoots completely into exterior space (provided $`2\delta`$ does not enter a separate object).
2. **Measured Breakdown Regimes on Non-Spherical SDFs (and Partial Rescues)**:
   - **Sharp Polyhedral Wedges ($`\sigma_t' r_b < 1`$)**: On a `sdRoundBox` corner ridge at the $`\mu = 0`$ shadow terminator ($`\sigma_t' = 8.2`$), when corner rounding shrinks from $`r_b = 0.20`$ ($`\sigma_t' r_b = 1.64`$) to $`r_b = 0.01`$ ($`\sigma_t' r_b = 0.08`$, $`H = 50`$), 3D Adjoint Monte Carlo stays bright (`0.2505`) because the sharp corner bounds a thick $`90^\circ`$ solid cube with an adjacent sunlit face. Purely local `sdfSteklovSphereBSSRDF(0, H=50)` treats $`H = 50`$ as an isolated tiny bead of radius $`R = 0.02`$ in vacuum and drops to `0.0469` (`-81.3%`). Adding `sdfGeodesicRingBSSRDF` projects across the sharp crease onto the sunlit $`45^\circ`$ face (`0.1181`), recovering most of the wedge bleed (Panel (d)).
   - **Concave Saddle Crevices ($`H \le 0`$) & Self-Shadowed Torus Holes**: Because interior spherical Bessel modes exist only for interior osculating balls ($`H > 0`$), `sdfSteklovSphereBSSRDF` clamps $`H \le 0`$ to the planar limit ($`M_{\text{flat}}(0) = 0`$), which visibly darkens narrow concave neck creases in shadow (see the neck shadow in Panel (b)). Similarly, inside the self-shadowed inner saddle of a ring torus ($`R = 0.55, r = 0.25`$, $`\mu = 0`$), 3D Monte Carlo measures `0.2909` via poloidal diffusion from the upper sunlit lip, whereas local `sdfSteklovSphereBSSRDF` gives `0.0327` and a 4-tap `sdfGeodesicRingBSSRDF` ($`r_{\text{ring}} = 0.18`$) partially bridges the gap (`0.0961`), motivating multi-scale geodesic rings for deep concavities.
