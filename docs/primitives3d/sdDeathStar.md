# `sdDeathStar` — Exact 3D Sphere Subtracted by Sphere (Death Star)

| 3D Surface ($f(p) = 0$) | 2D Cross-Section Slice |
| :---: | :---: |
| ![sdDeathStar](../images/primitives3d/sdDeathStar.png) | ![sdDeathStar slice](../images/primitives3d/sdDeathStar_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdDeathStar.glsl`](../../glsl/primitives3d/sdDeathStar.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`) |

---

## Mathematical Overview

Exact Euclidean SDF of a sphere of radius $r_a$ with a second sphere of radius $r_b$ at $(d,0,0)$ subtracted from it.

---

## GLSL Implementation

```glsl
float sdDeathStar( vec3 p2, float ra, float rb, float d )
{
    float a = (ra*ra - rb*rb + d*d)/(2.0*d);
    float b = sqrt(max(ra*ra - a*a, 0.0));
    vec2 p = vec2( p2.x, length(p2.yz) );
    if( p.x*b - p.y*a > d*max(b - p.y, 0.0) )
        return length(p - vec2(a,b));
    else
        return max( (length(p) - ra), -(length(p - vec2(d,0.0)) - rb));
}
```
