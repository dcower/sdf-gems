# `sdRoundCone` — Exact 3D Round Cone (Swept Sphere-Tapered Capsule)

| 3D Surface ($f(p) = 0$) | 2D Cross-Section Slice |
| :---: | :---: |
| ![sdRoundCone](../images/primitives3d/sdRoundCone.png) | ![sdRoundCone slice](../images/primitives3d/sdRoundCone_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdRoundCone.glsl`](../../glsl/primitives3d/sdRoundCone.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`) |

---

## Mathematical Overview

Computes the exact Euclidean distance to the convex hull of two spheres of radii $r_1, r_2$ separated by $h$.

---

## GLSL Implementation

```glsl
float sdRoundCone( vec3 p, float r1, float r2, float h )
{
    float b = (r1 - r2)/h;
    float a = sqrt(1.0 - b*b);
    vec2 q = vec2( length(p.xz), p.y );
    float k = dot(q, vec2(-b, a));
    if( k < 0.0 ) return length(q) - r1;
    if( k > a*h ) return length(q - vec2(0.0, h)) - r2;
    return dot(q, vec2(a, b)) - r1;
}
```
