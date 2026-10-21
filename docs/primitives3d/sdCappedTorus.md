# `sdCappedTorus` — Exact 3D Capped Torus Arc

| 3D Surface ($f(p) = 0$) | Primary 2D Cross-Section | Secondary 2D Cross-Section |
| :---: | :---: | :---: |
| ![sdCappedTorus](../images/primitives3d/sdCappedTorus.png) | ![sdCappedTorus slice 1](../images/primitives3d/sdCappedTorus_slice.png) | ![sdCappedTorus slice 2](../images/primitives3d/sdCappedTorus_slice2.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdCappedTorus.glsl`](../../glsl/primitives3d/sdCappedTorus.glsl) |
| **Exact Euclidean SDF?** | Yes |

---

## Mathematical Overview

Computes the Euclidean distance to a circular arc of major radius $r_a$, tube radius $r_b$, and half-angle $(\sin\theta, \cos\theta) = sc$.

---

## GLSL Implementation

```glsl
float sdCappedTorus( vec3 p, vec2 sc, float ra, float rb )
{
    p.x = abs(p.x);
    float k = (sc.y*p.x > sc.x*p.y) ? dot(p.xy, sc) : length(p.xy);
    return sqrt( dot(p,p) + ra*ra - 2.0*ra*k ) - rb;
}
```
