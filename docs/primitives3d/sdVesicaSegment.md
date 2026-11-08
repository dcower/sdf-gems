# `sdVesicaSegment` & `sdgVesicaSegment` — Exact 1-to-2 Sqrt 3D Revolved Vesica / Spindle + Analytical Gradient

| 3D Surface ($`f(p) = 0`$) | 2D Cross-Section Slice |
| :---: | :---: |
| ![sdVesicaSegment](../images/primitives3d/sdVesicaSegment.png) | ![sdVesicaSegment slice](../images/primitives3d/sdVesicaSegment_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez — upgraded/optimized using Gemini 4 Argon (`sdf-gems`, 2026)** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdVesicaSegment.glsl`](../../glsl/primitives3d/sdVesicaSegment.glsl) |
| **Exact Euclidean SDF?** | Yes ($`\Vert \nabla d\Vert = 1`$, `1` `sqrt`/`rsq` in tip cones & `2` in spindle body, down from `3` `sqrt`s everywhere) |

---

## Mathematical Overview

Inigo Quilez's baseline `sdVesicaSegment(p, a, b, w)` calls **3 square roots** on every evaluation:
1. `l = length(b - a)` (`sqrt` #1) to normalize `v = (b - a) / l` and compute $`r = \tfrac{1}{2}l`$,
2. `q.x = length(p - c - y*v)` (`sqrt` #2) to extract the radial distance from the segment axis,
3. `length(q - h.xy)` (`sqrt` #3) to compute the final distance in the $`(q_x, q_y)`$ meridian plane.

We eliminate **1 to 2 of the 3 square roots** in both `sdVesicaSegment` and its analytical gradient companion `sdgVesicaSegment`:

1. **`l = length(b - a)` Cancels Identically Out of the Entire Function**:
   Let $`b_a = b - a`$, $`p_c = p - \tfrac{1}{2}(a + b)`$, $`l_2 = \|b_a\|^2 = 4 r^2`$, and $`u = \mathrm{dot}(p_c, b_a) = l y`$.
   - The arc offset $`d = \frac{r^2 - w^2}{2 w} = \frac{l_2 - 4 w^2}{8 w}`$ depends only on $`l_2`$, not $`l = \sqrt{l_2}`$.
   - The orthogonal radial vector is $`\mathbf{u}_\perp = p_c - (u / l_2) b_a`$ with squared norm $`x_2 = \|\mathbf{u}_\perp\|^2`$.
   - Multiplying the tip-cone inequality $`r q_x < d(q_y - r)`$ by $`2 l = 4 r > 0`$ gives:
     $`\displaystyle l_2 \sqrt{x_2} < d\bigl(2 |u| - l_2\bigr) = z \iff z > 0 \;\wedge\; l_2^2 x_2 < z^2,`$
     which is completely free of both $`\sqrt{l_2}`$ and $`\sqrt{x_2}`$!
2. **Pythagorean Recombination & Factored 3D Gradient (`sdgVesicaSegment`)**:
   - **Tip-Cone Branch (`1` `sqrt` / `inversesqrt` total)**: The closest point on the spindle is the nearer endpoint $`a`$ or $`b`$, with displacement $`\mathbf{v}_t = p_c - 0.5\mathrm{sign}(u)b_a`$. A single `inversesqrt(dot(vt, vt))` yields both the exact distance and exact 3D unit gradient!
   - **Spindle-Body Branch (`2` `sqrt` / `inversesqrt`s total)**: With radial unit vector $`\mathbf{e}_r = \mathbf{u}_\perp / \sqrt{x_2}`$, the 3D vector from the revolved virtual circle center $`-d\,\mathbf{e}_r`$ to $`p_c`$ is $`\mathbf{w}_v = p_c + d\,\mathbf{e}_r`$, whose squared length is $`\|p_c\|^2 + d^2 + 2 d\sqrt{x_2}`$. A second `inversesqrt` normalizes $`\mathbf{w}_v`$ to give the exact 3D gradient and distance simultaneously!

![Vesica 1-to-2 Sqrt Reduction](../images/comparisons/vesica_comparison.png)

---

## GLSL Implementation

See [`glsl/primitives3d/sdVesicaSegment.glsl`](../../glsl/primitives3d/sdVesicaSegment.glsl) for `sdVesicaSegment` and `sdgVesicaSegment`.
