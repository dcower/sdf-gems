# `opBend` — 3D Bend Deformation (opCheapBend)

| 3D Surface ($f(p) = 0$) | Primary 2D Cross-Section | Secondary 2D Cross-Section |
| :---: | :---: | :---: |
| ![opBend](../images/deform/opBend.png) | ![opBend slice 1](../images/deform/opBend_slice.png) | ![opBend slice 2](../images/deform/opBend_slice2.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `deform` |
| **GLSL File** | [`glsl/deform/opBend.glsl`](../../glsl/deform/opBend.glsl) |
| **Exact Euclidean SDF?** | No (Shear distortion) |

---

## Mathematical Overview

Shear-rotates $p_{xy}$ by angle $k p_x$.

---

## GLSL Implementation

```glsl
vec3 opCheapBend( in vec3 p, in float k )
{
    float c = cos(k*p.x);
    float s = sin(k*p.x);
    mat2  m = mat2(c, -s, s, c);
    return vec3(m*p.xy, p.z);
}
```
