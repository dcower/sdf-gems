# `sdCutHollowSphere` — Exact 3D Cut Hollow Spherical Shell

| 3D Surface ($f(p) = 0$) | 2D Cross-Section Slice |
| :---: | :---: |
| ![sdCutHollowSphere](../images/primitives3d/sdCutHollowSphere.png) | ![sdCutHollowSphere slice](../images/primitives3d/sdCutHollowSphere_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdCutHollowSphere.glsl`](../../glsl/primitives3d/sdCutHollowSphere.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`) |

---

## Mathematical Overview

Spherical shell of radius $r$ and thickness $t$ truncated at plane $y = h$.

---

## GLSL Implementation

```glsl
float sdCutHollowSphere( vec3 p, float r, float h, float t )
{
    float w = sqrt(r*r - h*h);
    vec2 q = vec2( length(p.xz), p.y );
    return ((h*q.x < w*q.y) ? length(q - vec2(w,h)) : abs(length(q) - r)) - t;
}
```
