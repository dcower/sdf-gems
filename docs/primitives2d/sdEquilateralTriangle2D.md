# `sdEquilateralTriangle2D` — Exact 0/1-`sqrt` 2D Equilateral Triangle

![sdEquilateralTriangle2D](../images/primitives2d/sdEquilateralTriangle2D.png)

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez — upgraded to 0-`sqrt` interior/face-strip using Gemini 4 Argon (`sdf-gems`, 2026)** |
| **Source** | [https://iquilezles.org/articles/distfunctions2d/](https://iquilezles.org/articles/distfunctions2d/) |
| **Category** | `primitives2d` |
| **GLSL File** | [`glsl/primitives2d/sdEquilateralTriangle2D.glsl`](../../glsl/primitives2d/sdEquilateralTriangle2D.glsl) |
| **Exact Euclidean SDF?** | Yes (`\|\|grad d\|\| = 1`, `0` `sqrt`s in interior & face strips, `1` in corner wedges) |

---

## Mathematical Overview

After folding $`p \in \mathbb{R}^2`$ into the $`D_3`$ fundamental chamber via `p.x = abs(p.x) - r; p.y += r/k; if (p.x + k*p.y > 0.0) p = vec2(p.x - k*p.y, -k*p.x - p.y)/2.0;`, the chamber is bounded by the two bisector rays $`p_x = -2r + k p_y`$ and $`p_x = -k p_y`$:
```math
-2r + k p_y \le p_x \le -k p_y \iff |p_x + r| - r \le -k p_y.
```
Therefore, for $`d_x = |p_x + r| - r`$, having $`d_x > 0`$ implies $`p_y < 0`$ automatically (the $`3`$ exterior corner wedges). Whenever $`d_x \le 0`$ (the entire triangle interior and all $`3`$ exterior flat-edge normal strips), the exact signed distance is simply $`-p_y`$, eliminating `clamp`, `sign`, and the `sqrt` (`length`) across the interior and face strips (`4.44e-16` max error).

---

## GLSL Implementation

```glsl
float sdEquilateralTriangle2D( vec2 p, float r )
{
    const float k = sqrt(3.0);
    p.x = abs(p.x) - r;
    p.y = p.y + r/k;
    if( p.x + k*p.y > 0.0 ) p = vec2(p.x - k*p.y, -k*p.x - p.y)/2.0;
    float dx = abs(p.x + r) - r;
    return (dx > 0.0) ? length(vec2(dx, p.y)) : -p.y;
}
```
