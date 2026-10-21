# `sdBoxFrame` — Exact 3D Wireframe Box Frame

| 3D Surface ($f(p) = 0$) | Primary 2D Cross-Section | Secondary 2D Cross-Section |
| :---: | :---: | :---: |
| ![sdBoxFrame](../images/primitives3d/sdBoxFrame.png) | ![sdBoxFrame slice 1](../images/primitives3d/sdBoxFrame_slice.png) | ![sdBoxFrame slice 2](../images/primitives3d/sdBoxFrame_slice2.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdBoxFrame.glsl`](../../glsl/primitives3d/sdBoxFrame.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`) |

---

## Mathematical Overview

Evaluates the exact Euclidean distance to the 12 edges of a rectangular box frame of half-extents $b$ and bar thickness $e$ by taking the minimum of the three axis-aligned edge struts in the positive octant.

---

## GLSL Implementation

```glsl
float sdBoxFrame( vec3 p, vec3 b, float e )
{
    p = abs(p) - b;
    vec3 q = abs(p + e) - e;
    return min(min(
        length(max(vec3(p.x, q.y, q.z), 0.0)) + min(max(p.x, max(q.y, q.z)), 0.0),
        length(max(vec3(q.x, p.y, q.z), 0.0)) + min(max(q.x, max(p.y, q.z)), 0.0)),
        length(max(vec3(q.x, q.y, p.z), 0.0)) + min(max(q.x, max(q.y, p.z)), 0.0));
}
```
