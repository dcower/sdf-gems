# `sdRoundedCylinder` — 0/1/2-Sqrt Exact 3D Rounded Cylinder

| 3D Surface ($`f(p) = 0`$) | 2D Cross-Section Slice |
| :---: | :---: |
| ![sdRoundedCylinder](../images/primitives3d/sdRoundedCylinder.png) | ![sdRoundedCylinder slice](../images/primitives3d/sdRoundedCylinder_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez — upgraded to 0/1/2 `sqrt`s via Cap-Voronoi Radius Theorem using Gemini 4 Argon (sdf-gems, 2026)** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdRoundedCylinder.glsl`](../../glsl/primitives3d/sdRoundedCylinder.glsl) |
| **Exact Euclidean SDF?** | Yes (`\|\|grad d\|\| = 1`, `0` `sqrt`s in cap zones, `1` in wall zone, `2` in rim wedge) |

---

## Mathematical Overview

Applies the **Unified Cap-Voronoi Radius Theorem** to the inner core cylinder of radius $`r = 2 r_a - r_b`$ and half-height $`h`$, evaluating with **0 `sqrt`s** in the exterior and interior flat-cap Voronoi zones, **1 `sqrt`** in the cylindrical wall zone, and **2 `sqrt`s** in the rounded corner rim.

---

## GLSL Implementation

```glsl
float sdRoundedCylinder( vec3 p, float ra, float rb, float h )
{
    float r = 2.0 * ra - rb, r2 = dot(p.xz, p.xz), wy = abs(p.y) - h, rc = r + min(wy, 0.0);
    if( rc >= 0.0 && r2 <= rc * rc ) return wy - rb;
    float wx = sqrt(r2) - r;
    return ((wy <= 0.0) ? wx : sqrt(wx * wx + wy * wy)) - rb;
}
```
