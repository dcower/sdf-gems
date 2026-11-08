# `sdCappedCylinderSeg` — 1/2-Sqrt Exact 3D Arbitrary-Segment Capped Cylinder

| 3D Surface ($`f(p) = 0`$) | 2D Cross-Section Slice |
| :---: | :---: |
| ![sdCappedCylinderSeg](../images/primitives3d/sdCappedCylinderSeg.png) | ![sdCappedCylinderSeg slice](../images/primitives3d/sdCappedCylinderSeg_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez — upgraded to 1/2 `sqrt`s (down from 3) & FP16-safe using Gemini 4 Argon (sdf-gems, 2026)** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdCappedCylinderSeg.glsl`](../../glsl/primitives3d/sdCappedCylinderSeg.glsl) |
| **Exact Euclidean SDF?** | Yes (`\|\|grad d\|\| = 1`, `1` `sqrt` in cap & wall zones, `2` `sqrt`s in rim wedge; down from `3` `sqrt`s) |

---

## Mathematical Overview

Let $`b_a = b - a`$, $`l_2 = \|b_a\|_2^2`$, $`t = \frac{(p-a)\cdot b_a}{l_2}`$, $`u = (p-a) - t b_a`$, $`x_2 = \|u\|_2^2`$, and $`u_y = |t - 0.5| - 0.5`$. Because the squared axial distance past the caps is $`w_y^2 = u_y^2 l_2`$ (which requires **0 `sqrt`s**), the segment length $`L = \sqrt{l_2}`$ is never needed in the cylindrical wall zone or the exterior circular rim wedge $`\sqrt{(\sqrt{x_2} - r)^2 + u_y^2 l_2}`$! Moreover, all intermediates scale as $`O(L^2)`$ instead of $`O(L^6)`$, preventing FP16 overflow while cutting the `sqrt` count from **3 down to 1 (cap & wall) or 2 (rim)**.

---

## GLSL Implementation

```glsl
float sdCappedCylinderSeg( vec3 p, vec3 a, vec3 b, float r )
{
    vec3  ba = b - a, pa = p - a;
    float l2 = dot(ba, ba), t = dot(pa, ba) / l2;
    vec3  u  = pa - t * ba;
    float x2 = dot(u, u), uy = abs(t - 0.5) - 0.5;
    if( uy >= 0.0 ) {
        if( x2 <= r * r ) return uy * sqrt(l2);
        float wx = sqrt(x2) - r;
        return sqrt(wx * wx + uy * uy * l2);
    }
    if( x2 > r * r ) return sqrt(x2) - r;
    float wx = sqrt(x2) - r, wy2 = uy * uy * l2;
    return (wx * wx <= wy2) ? wx : uy * sqrt(l2);
}
```
