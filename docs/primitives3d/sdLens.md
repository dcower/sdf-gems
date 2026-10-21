# `sdLens` — Exact 3D Asymmetric Two-Sphere Optical Lens

| 3D Surface ($f(p) = 0$) | 2D Cross-Section Slice |
| :---: | :---: |
| ![sdLens](../images/primitives3d/sdLens.png) | ![sdLens slice](../images/primitives3d/sdLens_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Dillon Cower (sdf-gems, 2026)** |
| **Source** | [https://github.com/dcower/sdf-gems](https://github.com/dcower/sdf-gems) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdLens.glsl`](../../glsl/primitives3d/sdLens.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`) |

---

## Mathematical Overview

Computes the exact 3D Euclidean distance to the intersection of two spheres of radii $r_a, r_b$ separated by $d$ (the intersection dual to `sdDeathStar`).

---

## GLSL Implementation

```glsl
float sdLens( vec3 p2, float ra, float rb, float d )
{
    float a = (ra*ra - rb*rb + d*d)/(2.0*d);
    float b = sqrt(max(ra*ra - a*a, 0.0));
    vec2  p = vec2( p2.x, length(p2.yz) );
    float k = p.x*b - p.y*a;
    if( k < 0.0 && k > d*(b - p.y) )
        return length(p - vec2(a, b));
    else
        return max( length(p) - ra, length(p - vec2(d, 0.0)) - rb );
}
```
