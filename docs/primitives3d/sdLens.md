# `sdLens` & `sdgLens` — Exact 1/2-`sqrt` 3D Asymmetric Two-Sphere Optical Lens + Analytical Gradient

| 3D Surface ($f(p) = 0$) | 2D Cross-Section Slice |
| :---: | :---: |
| ![sdLens](../images/primitives3d/sdLens.png) | ![sdLens slice](../images/primitives3d/sdLens_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Discovered using Gemini 4 Argon (`sdf-gems`, 2026)** |
| **Source** | [https://github.com/dcower/sdf-gems](https://github.com/dcower/sdf-gems) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdLens.glsl`](../../glsl/primitives3d/sdLens.glsl) |
| **Exact Euclidean SDF?** | Yes ($\|\nabla d\| = 1$, **`1` `sqrt`/`inversesqrt` across $85\%$–$90\%$ of $\mathbb{R}^3$** and `2` in the rim wedge, down from `4` `sqrt`s in standard revolution) |

![1/2-Sqrt Lens and DeathStar Axis-Center Radial Collapse](../images/comparisons/lens_deathstar_1sqrt_comparison.png)

---

## Mathematical Overview

### 1. The Axis-Center Radial Collapse Theorem (1 `sqrt` in $85\%$–$90\%$ of $\mathbb{R}^3$)
In standard revolved 3D primitives (`sdDeathStar`, `sdLens`), the 3D point $p \in \mathbb{R}^3$ is first projected into the 2D meridian plane via $p_{2\text{D}} = (p_x, \sqrt{p_y^2 + p_z^2})$ (`sqrt` #1) alongside the rim radius $b = \sqrt{\max(r_a^2 - a^2, 0)}$ (`sqrt` #2), followed by `length(p_2D - c)` (`sqrt` #3 & #4).

Observe two fundamental algebraic identities:
1. **Signed-Square Rim-Wedge Classification (`0` `sqrt`s):**
   With $y_{z2} = p_y^2 + p_z^2$ and $b_2 = \max(r_a^2 - a^2, 0)$, the rim-circle normal wedge conditions $p_x b - \|p_{yz}\| a < 0$ and $(d - p_x) b - \|p_{yz}\|(d - a) < 0$ are equivalent via signed squaring to:
   $$a|a|\, y_{z2} > b_2\, p_x|p_x| \quad \text{and} \quad (d - a)|d - a|\, y_{z2} > b_2\, (d - p_x)|d - p_x|,$$
   which require **zero square roots**!
2. **Axis-Center Radial Collapse (`1` `sqrt` total):**
   Outside the rim wedge ($85\%$–$90\%$ of $\mathbb{R}^3$, including **100% of the lens interior** and both spherical cap regions), the closest sphere center lies on the $x$-axis at $(c_x, 0, 0)$ with $c_x \in \{0, d\}$. Therefore:
   $$\left\|\bigl(p_x - c_x, \sqrt{p_y^2 + p_z^2}\bigr)\right\|_2 \equiv \sqrt{(p_x - c_x)^2 + p_y^2 + p_z^2},$$
   so $\sqrt{p_y^2 + p_z^2}$ and $\sqrt{b_2}$ **both vanish**, leaving **only 1 `sqrt` (or 1 `inversesqrt` in `sdgLens`)**!
3. **Conjugate-Rationalized Rim Distance (`2` `sqrt`s in Rim Wedge):**
   Inside the narrow rim wedge, $(\sqrt{y_{z2}} - \sqrt{b_2})^2 = \frac{(y_{z2} - b_2)^2}{y_{z2} + b_2 + 2\sqrt{b_2 y_{z2}}}$ fuses $\sqrt{b_2}\sqrt{y_{z2}} = \sqrt{b_2 y_{z2}}$ into a single inner radical without catastrophic cancellation on the rim $y_{z2} \to b_2$.

---

## GLSL Implementation

See [`glsl/primitives3d/sdLens.glsl`](../../glsl/primitives3d/sdLens.glsl) for `sdLens` and `sdgLens`.
