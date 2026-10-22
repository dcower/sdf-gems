# `sdCutSphere` & `sdgCutSphere` — Exact 0/1/2-SFU 3D Cut Sphere & Cancellation-Free Conjugate Rim Gradient

| 3D Surface ($f(p) = 0$) | 2D Cross-Section Slice |
| :---: | :---: |
| ![sdCutSphere](../images/primitives3d/sdCutSphere.png) | ![sdCutSphere slice](../images/primitives3d/sdCutSphere_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez — upgraded to unified 0-`sqrt` cap-Voronoi & 2-SFU rim gradient using Gemini 4 Argon (sdf-gems, 2026)** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdCutSphere.glsl`](../../glsl/primitives3d/sdCutSphere.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`, `0` SFUs in exterior + interior cap zone, `1` SFU in sphere, `2` SFUs in rim) |

---

## Mathematical Overview

![Cut Sphere & Spherical Segment Cap-Voronoi Comparison](../images/comparisons/cut_sphere_cap_voronoi_comparison.png)

1. **Unified 0-`sqrt` Exterior + Interior Cap-Voronoi Tri-Partition**:
   With $x^2 = p_x^2 + p_z^2$ and $d_y = h - p_y$, the flat cap disk at $y = h$ is the closest boundary both in the exterior cylinder ($d_y \ge 0, x^2 \le r^2 - h^2$) and inside the cap paraboloid ($d_y < 0, \sqrt{x^2 + p_y^2} \le r + d_y \iff x^2 \le (r + h)(r - h + 2 d_y)$, where $p_y^2$ cancels identically). Combining both into:
   $$x^2 \le (r + h)\big(r - h + 2\min(d_y, 0)\big)$$
   resolves the entire exterior and interior cap Voronoi zone first with **0 `sqrt`s**, and guarantees $x^2 > w^2$ in the rim branch so the sphere test reduces to the single signed-square inequality $h|h|x^2 < w^2 p_y |p_y|$.
2. **Cancellation-Free Conjugate Rim Gradient (`sdgCutSphere`, 2 SFUs instead of 3)**:
   On the circular rim of radius $w = \sqrt{w^2}$ at $y = h$, the radial displacement vector is:
   $$\mathbf{v}_{xz} = \left(1 - \frac{w}{x}\right)\mathbf{p}_{xz} = \frac{x^2 - w^2}{x^2 + \sqrt{w^2 x^2}}\mathbf{p}_{xz},$$
   which eliminates subtractive cancellation as $x \to w$ and fuses $\sqrt{w^2}$ and $\sqrt{x^2}$ into a single `sqrt(w2 * x2)`!

---

## GLSL Implementation

See [`glsl/primitives3d/sdCutSphere.glsl`](../../glsl/primitives3d/sdCutSphere.glsl).
