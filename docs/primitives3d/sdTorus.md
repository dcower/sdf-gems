# `sdTorus` — Exact 3D Torus

| 3D Surface ($f(p) = 0$) | Primary 2D Cross-Section | Secondary 2D Cross-Section |
| :---: | :---: | :---: |
| ![sdTorus](../images/primitives3d/sdTorus.png) | ![sdTorus slice 1](../images/primitives3d/sdTorus_slice.png) | ![sdTorus slice 2](../images/primitives3d/sdTorus_slice2.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdTorus.glsl`](../../glsl/primitives3d/sdTorus.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`) |

---

## Mathematical Overview

Reduces the 3D torus of major radius $t_x$ and minor radius $t_y$ via radial revolution $q = (\|p_{xz}\| - t_x, p_y)$ to a 2D circle of radius $t_y$.

---

## GLSL Implementation

```glsl
float sdTorus( vec3 p, vec2 t )
{
    vec2 q = vec2(length(p.xz) - t.x, p.y);
    return length(q) - t.y;
}
```
