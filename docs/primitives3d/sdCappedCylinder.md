# `sdCappedCylinder` — Exact 3D Vertical Capped Cylinder

| 3D Surface ($f(p) = 0$) | 2D Cross-Section Slice |
| :---: | :---: |
| ![sdCappedCylinder](../images/primitives3d/sdCappedCylinder.png) | ![sdCappedCylinder slice](../images/primitives3d/sdCappedCylinder_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdCappedCylinder.glsl`](../../glsl/primitives3d/sdCappedCylinder.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`) |

---

## Mathematical Overview

Revolves a 2D box of radius $r$ and half-height $h$ around the $y$-axis.

---

## GLSL Implementation

```glsl
float sdCappedCylinder( vec3 p, float h, float r )
{
    vec2 d = abs(vec2(length(p.xz), p.y)) - vec2(r, h);
    return min(max(d.x, d.y), 0.0) + length(max(d, 0.0));
}
```
