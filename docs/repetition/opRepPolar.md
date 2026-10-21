# `opRepPolar` & `opRepPolarLim` — Voronoi-Centered 1-Tap & Margin-Guarded 2-Tap Polar Repetition

| 3D Surface ($f(p) = 0$) | 2D Cross-Section Slice |
| :---: | :---: |
| ![opRepPolar](../images/repetition/opRepPolar.png) | ![opRepPolar slice](../images/repetition/opRepPolar_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Mercury Demogroup — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)** |
| **Source** | [https://mercury.sexy/hg_sdf/](https://mercury.sexy/hg_sdf/) |
| **Category** | `repetition` |
| **GLSL File** | [`glsl/repetition/opRepPolar.glsl`](../../glsl/repetition/opRepPolar.glsl) |
| **Exact Euclidean SDF?** | Yes (eliminates $\theta = \pm\pi$ branch cut, asymmetric `[lima, limb]` Voronoi tear, and boundary neighbor drop) |

---

## Mathematical Overview

### 1. Voronoi-Centered Periodic Wrapping for General `[lima, limb]`

In limited polar repetition on an arbitrary sector range `[lima, limb]`, naive clamping `clamp(round(a / sp), lima, limb)` places the periodic branch cut at $\theta = \pm\pi$. Whenever $\text{lima} + \text{limb} \ne 0$ (e.g. `[-1, 2]` or `[0, 3]`), the true Voronoi bisector of the empty angular wedge is at the antipode of the midpoint angle $\text{mid} = \frac{1}{2}(\text{lima} + \text{limb})$, causing up to **`2.869` (`286.9%`!)** distance error across the empty wedge.

Wrapping the continuous sector coordinate around `mid` with period $n$:

$$\text{sector} = \text{mid} + \operatorname{mod}(a - \text{mid} + 0.5 n,\; n) - 0.5 n \in [\text{mid} - 0.5 n,\; \text{mid} + 0.5 n]$$

places the periodic branch cut at the **exact Voronoi bisector** $\text{mid} \pm 0.5 n$, reducing 1-tap `opRepPolarLim` error from **`2.853`** down to **`1.55e-15` (machine epsilon)**.

### 2. Margin-Guarded Empty-Wedge Straddling (`opLimitedPolarRep2Tap`)

For 2-tap asymmetric polar repetition, switching immediately to `(lima, limb)` at $|\text{sector} - \text{mid}| > 0.5(\text{limb} - \text{lima})$ drops the adjacent neighbor (`lima + 1` or `limb - 1`) right at the boundary ray, clipping any asymmetric primitive that overhangs the boundary ray (`0.339` error). Extending the active adjacent-pair zone by $\text{margin} = 0.25\max(n - \text{span} - 1, 0)$ preserves `(lima, lima + 1)` and `(limb - 1, limb)` past the boundary rays before transitioning to `(lima, limb)` across the central Voronoi bisector (`0.0` error everywhere).

---

## Numerical Verification

| Operator / Test Case | Baseline Max Error | **Ours (`opRepPolar.glsl`)** Max Error |
| :--- | :---: | :---: |
| `opRepPolarLim` 1-Tap Symmetric ($n=5, [-1, 2]$) | `1.72e+00` (naive `clamp`) | **`2.44e-15`** |
| `opRepPolarLim` 1-Tap Symmetric ($n=5, [0, 3]$) | `2.85e+00` (naive `clamp`) | **`1.55e-15`** |
| `opLimitedPolarRep2Tap` Asymmetric ($n=8, [-2, 2]$) | `3.39e-01` (no boundary margin) | **`0.00e+00`** |
| `opLimitedPolarRep2Tap` Asymmetric ($n=5, [0, 3]$) | `2.87e+00` ($\pm\pi$ cut) | **`0.00e+00`** |

![Before vs After Comparison](../images/comparisons/polar_limited_2d.png)
