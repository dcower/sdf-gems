# `sdCone` — Exact 1-to-2 Sqrt 3D Finite & Infinite Cone

| 3D Surface ($f(p) = 0$) | 2D Cross-Section Slice ($z = 0$) |
| :---: | :---: |
| ![sdCone](../images/primitives3d/sdCone.png) | ![sdCone slice](../images/primitives3d/sdCone_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez — upgraded/optimized using Gemini 4 Argon (sdf-gems, 2026)** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdCone.glsl`](../../glsl/primitives3d/sdCone.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`, `1` `sqrt` everywhere except base-rim wedge & `1` division, down from `2` `sqrt`s & `3` divisions) |

---

## Mathematical Overview

Inigo Quilez's finite `sdCone(p, c, h)` receives the unit aperture vector `c = vec2(sin(alpha), cos(alpha))` ($\|c\| = 1$), un-normalizes it to `q = h * vec2(c.x/c.y, -1.0)` (division #1), projects `w = vec2(length(p.xz), p.y)` onto `q` via `dot(w, q) / dot(q, q)` (division #2) and `w.x / q.x` (division #3), and evaluates **2 `sqrt`s everywhere** (`length(p.xz)` and `sqrt(d)`).

Because $\|c\| = 1$:
1. **In the Top Apex Cone ($p_y \ge 0 \;\wedge\; c_x^2 x_2 \le c_y^2 p_y^2$)**: The closest point is the apex $(0, 0, 0)$ at distance $\sqrt{x_2 + p_y^2}$, which does not even require $\sqrt{x_2}$ (**1 `sqrt` total**).
2. **In the Mantle, Base Disk, and Interior**: The signed normal distance to the cone mantle is linear (`c.y * wx + c.x * p.y` where $w_x = \sqrt{x_2}$) and the signed distance to the base cap is linear (`d = -p.y - h`), so their maximum `max(c.y * wx + c.x * p.y, d)` is exact with **only 1 `sqrt` total**. A second `sqrt` (`length(vec2(rx, d))`) is needed only inside the narrow exterior base-rim wedge ($r_x > 0 \;\wedge\; c_x r_x + c_y d \ge 0$).
3. **Infinite Cone (`sdConeInf`)**: Collapses from 2 `sqrt`s to **1 `sqrt` everywhere in $\mathbb{R}^3$** (`4.44e-16` max error).

---

## GLSL Implementation

```glsl
// c = vec2(sin(alpha), cos(alpha)) normalized cone aperture, h = cone height
float sdCone( vec3 p, vec2 c, float h )
{
    float x2 = dot(p.xz, p.xz);
    if( p.y >= 0.0 && c.x * c.x * x2 <= c.y * c.y * p.y * p.y )
        return sqrt(x2 + p.y * p.y);
    float wx = sqrt(x2), rx = wx - h * c.x / c.y, d = -p.y - h;
    return (rx > 0.0 && c.x * rx + c.y * d >= 0.0)
         ? length(vec2(rx, d))
         : max(c.y * wx + c.x * p.y, d);
}

// Infinite 3D cone (1 `sqrt` everywhere in R^3, down from 2 `sqrt`s)
float sdConeInf( vec3 p, vec2 c )
{
    float x2 = dot(p.xz, p.xz);
    return (p.y >= 0.0 && c.x * c.x * x2 <= c.y * c.y * p.y * p.y)
         ? sqrt(x2 + p.y * p.y)
         : sqrt(x2) * c.y + p.y * c.x;
}
```
