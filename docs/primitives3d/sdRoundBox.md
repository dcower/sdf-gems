# `sdRoundBox` — Exact 3D Rounded Box

| 3D Surface ($f(p) = 0$) | Primary 2D Cross-Section | Secondary 2D Cross-Section |
| :---: | :---: | :---: |
| ![sdRoundBox](../images/primitives3d/sdRoundBox.png) | ![sdRoundBox slice 1](../images/primitives3d/sdRoundBox_slice.png) | ![sdRoundBox slice 2](../images/primitives3d/sdRoundBox_slice2.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdRoundBox.glsl`](../../glsl/primitives3d/sdRoundBox.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`) |

---

## Mathematical Overview

Offsets an inner box of half-extents $b - r$ by rounding radius $r$ via the Minkowski sum with a sphere of radius $r$.

---

## GLSL Implementation

```glsl
float sdRoundBox( vec3 p, vec3 b, float r )
{
    vec3 q = abs(p) - b + r;
    return length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0) - r;
}
```
