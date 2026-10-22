# `sdCappedCone` — Exact 0-Division, 0/1/2-Sqrt 3D Capped Cone

| 3D Surface ($f(p) = 0$) | 2D Cross-Section Slice ($z = 0$) |
| :---: | :---: |
| ![sdCappedCone](../images/primitives3d/sdCappedCone.png) | ![sdCappedCone slice](../images/primitives3d/sdCappedCone_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez — upgraded to 0-div, 0-clamp, 0/1/2-`sqrt` using Gemini 4 Argon (sdf-gems, 2026)** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdCappedCone.glsl`](../../glsl/primitives3d/sdCappedCone.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`, `0` `sqrt`s in flat-cap Voronoi zone, `1` coord `sqrt` on slanted wall, `0` divs) |

---

## Mathematical Overview

![Capped Cone & Gradient Comparison](../images/comparisons/capped_cone_gradient_comparison.png)

1. **0-`sqrt` Exterior & Interior Cap-Voronoi Early Exit**:
   Outside the caps ($c_y = |p_y| - h \ge 0$), the exact Cap-Voronoi cylinder is $x^2 = p_x^2 + p_z^2 \le r_c^2$. Inside the cone ($c_y < 0$), using $\sqrt{k_x^2 + 4h^2} \le |k_x| + 2h$ yields the unified 0-`sqrt`, 0-division test $h_{\text{rc}} = h r_c + \min(c_y, 0)(h + \max(s_y k_x, 0)) \ge 0$ and $x^2 h^2 \le h_{\text{rc}}^2$, returning $c_y$ without even computing `length(p.xz)`.
2. **0-Division, 0-Clamp Slanted-Wall Partition**:
   With $l = \sqrt{x^2}$, $\mathbf{q} = (l - r_2, p_y - h)$, and signed cross product $w_c = 2h q_x - k_x q_y$, a point lies inside the frustum iff $c_y \le 0 \;\land\; w_c \le 0$ (where the cap wins iff $c_y^2 m \le w_c^2$ with **0 divisions and 0 clamps**), and lies in the exterior slanted-wall normal strip iff $q_k = \mathbf{q} \cdot \mathbf{k} \in (-m, 0)$, where $d = w_c / \sqrt{m}$ requires only **1 coordinate `sqrt`**.

---

## GLSL Implementation

See [`glsl/primitives3d/sdCappedCone.glsl`](../../glsl/primitives3d/sdCappedCone.glsl).
