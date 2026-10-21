# `sdTriPrism` — 3D Triangular Prism (Plane Bound)

| 3D Surface ($f(p) = 0$) | Primary 2D Cross-Section | Secondary 2D Cross-Section |
| :---: | :---: | :---: |
| ![sdTriPrism](../images/primitives3d/sdTriPrism.png) | ![sdTriPrism slice 1](../images/primitives3d/sdTriPrism_slice.png) | ![sdTriPrism slice 2](../images/primitives3d/sdTriPrism_slice2.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdTriPrism.glsl`](../../glsl/primitives3d/sdTriPrism.glsl) |
| **Exact Euclidean SDF?** | No (Plane Bound) |

---

## Mathematical Overview

Evaluates the maximum of the bounding half-space planes of an equilateral triangular prism.

---

## GLSL Implementation

```glsl
float sdTriPrism( vec3 p, vec2 h )
{
    vec3 q = abs(p);
    return max(q.z - h.y, max(q.x*0.866025 + p.y*0.5, -p.y) - h.x*0.5);
}
```
