# `sdSphericalZone` — Exact 3D Hollow Spherical Zone (Two-Plane Cut Shell)

| 3D Surface ($f(p) = 0$) | 2D Cross-Section Slice |
| :---: | :---: |
| ![sdSphericalZone](../images/primitives3d/sdSphericalZone.png) | ![sdSphericalZone slice](../images/primitives3d/sdSphericalZone_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Discovered using Gemini 4 Argon (sdf-gems, 2026)** |
| **Source** | [https://github.com/dcower/sdf-gems](https://github.com/dcower/sdf-gems) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdSphericalZone.glsl`](../../glsl/primitives3d/sdSphericalZone.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`) |

---

## Mathematical Overview

Exact Euclidean SDF of a spherical shell of radius $r$ and wall thickness $t$ truncated by two arbitrary planes $y = h_1$ and $y = h_2$.

---

## GLSL Implementation

```glsl
float sdSphericalZone( vec3 p, float r, float h1, float h2, float t )
{
    float w1 = sqrt(r*r - h1*h1);
    float w2 = sqrt(r*r - h2*h2);
    vec2  q  = vec2( length(p.xz), p.y );
    float d  = (h1*q.x > w1*q.y) ? length(q - vec2(w1, h1)) :
               (w2*q.y > h2*q.x) ? length(q - vec2(w2, h2)) :
                                   abs(length(q) - r);
    return d - t;
}
```
