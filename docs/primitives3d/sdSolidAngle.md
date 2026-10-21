# `sdSolidAngle` — Exact 3D Solid Angle (Spherical Sector)

| 3D Surface ($f(p) = 0$) | 2D Cross-Section Slice |
| :---: | :---: |
| ![sdSolidAngle](../images/primitives3d/sdSolidAngle.png) | ![sdSolidAngle slice](../images/primitives3d/sdSolidAngle_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdSolidAngle.glsl`](../../glsl/primitives3d/sdSolidAngle.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`) |

---

## Mathematical Overview

Intersection of a cone of half-angle $(\sin\alpha, \cos\alpha) = c$ with a sphere of radius $r_a$.

---

## GLSL Implementation

```glsl
float sdSolidAngle( vec3 p, vec2 c, float ra )
{
    vec2 q = vec2( length(p.xz), p.y );
    float l = length(q) - ra;
    float m = length(q - c*clamp(dot(q,c), 0.0, ra) );
    return max(l, m*sign(c.y*q.x - c.x*q.y));
}
```
