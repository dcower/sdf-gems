# `opTwist` — 3D Axial Twist Deformation

| 3D Surface ($f(p) = 0$) | Primary 2D Cross-Section | Secondary 2D Cross-Section |
| :---: | :---: | :---: |
| ![opTwist](../images/deform/opTwist.png) | ![opTwist slice 1](../images/deform/opTwist_slice.png) | ![opTwist slice 2](../images/deform/opTwist_slice2.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `deform` |
| **GLSL File** | [`glsl/deform/opTwist.glsl`](../../glsl/deform/opTwist.glsl) |
| **Exact Euclidean SDF?** | No (Distorts metric) |

---

## Mathematical Overview

Rotates $p_{xz}$ proportionally to $p_y$.

---

## GLSL Implementation

```glsl
vec3 opTwist( in vec3 p, in float k )
{
    float c = cos(k*p.y);
    float s = sin(k*p.y);
    mat2  m = mat2(c, -s, s, c);
    return vec3(m*p.xz, p.y);
}
```
