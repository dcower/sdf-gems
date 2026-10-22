# `sdCone` & `sdConeInf` — Exact 0/1/2-`sqrt` 3D Finite & 1-`sqrt` Infinite Cone (Slanted-Wall Cap-Voronoi Theorem)

| 3D Surface ($f(p) = 0$) | 2D Cross-Section Slice ($z = 0$) |
| :---: | :---: |
| ![sdCone](../images/primitives3d/sdCone.png) | ![sdCone slice](../images/primitives3d/sdCone_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez — upgraded to 0/1/2-`sqrt` via Slanted-Wall Cap-Voronoi Theorem using Gemini 4 Argon (`sdf-gems`, 2026)** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdCone.glsl`](../../glsl/primitives3d/sdCone.glsl) |
| **Exact Euclidean SDF?** | Yes ($\|\nabla d\| = 1$, **`0` `sqrt`s** in base-cap Voronoi, `1` in apex/mantle, `2` in rim wedge; down from `2` `sqrt`s & `3` divs in IQ) |

![0/1/2-Sqrt Cone, RoundedCylinder, and Rhombus Comparison](../images/comparisons/cone_rcyl_rhombus_comparison.png)

---

## Mathematical Overview

### 1. The Slanted-Wall Cap-Voronoi Radius Theorem (`0` `sqrt`s in Cap Voronoi)
Let $c = (\sin\alpha, \cos\alpha)$ be the unit cone aperture vector, $h$ be the cone height, $r_b = h\, c_x / c_y$ be the base radius, and $d = -p_y - h$ be the signed distance below the flat base disk $y = -h$.
- **Outside below the base disk ($d \ge 0$):** The closest point on the cone lies in the interior of the flat base disk whenever $p_x^2 + p_z^2 \le r_b^2$, with exact distance $d = -p_y - h$.
- **Inside the cone ($d < 0$):** The signed distance to the slanted mantle at $p_y = -h - d$ is:
  $$d_{\text{mantle}} = c_y \sqrt{x_2} + c_x p_y = c_y(\sqrt{x_2} - r_b) - c_x d.$$
  Thus the flat base disk dominates the slanted mantle ($d \ge d_{\text{mantle}}$) if and only if:
  $$(1 + c_x) d \ge c_y(\sqrt{x_2} - r_b) \iff \sqrt{x_2} \le r_b + \frac{1 + c_x}{c_y} d.$$

Unifying both exterior ($d \ge 0$) and interior ($d < 0$) into a single cap-Voronoi radius:
$$r_c = r_b + \frac{1 + c_x}{c_y}\min(d, 0) = \frac{h c_x + (1 + c_x)\min(d, 0)}{c_y},$$
whenever $r_c \ge 0$ and $p_x^2 + p_z^2 \le r_c^2$, the exact 3D Euclidean distance is simply **`d = -p.y - h` with `0` `sqrt`s**!

### 2. 1-`sqrt` Apex Cone & Mantle
- **Top Apex Cone ($p_y \ge 0 \;\wedge\; c_x^2 x_2 \le c_y^2 p_y^2$):** Closest point is the apex $(0,0,0)$ at distance $\sqrt{x_2 + p_y^2}$ (**1 `sqrt` total**, skipping $\sqrt{x_2}$).
- **Slanted Mantle:** Since the base-cap Voronoi region is already handled by $x_2 \le r_c^2$, the remaining interior and mantle region evaluates $c_y \sqrt{x_2} + c_x p_y$ (**1 `sqrt` total**), needing a 2nd `sqrt` only in the narrow exterior base-rim wedge.

---

## GLSL Implementation

See [`glsl/primitives3d/sdCone.glsl`](../../glsl/primitives3d/sdCone.glsl) (and [`glsl/gradients3d/sdgCone.glsl`](../../glsl/gradients3d/sdgCone.glsl) for `sdgCone` and `sdgConeInf`).
