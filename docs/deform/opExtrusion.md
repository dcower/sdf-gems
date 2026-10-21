# `opExtrusion` — Exact 2D-to-3D Linear Extrusion

| 3D Surface ($f(p) = 0$) | Primary 2D Cross-Section | Secondary 2D Cross-Section |
| :---: | :---: | :---: |
| ![opExtrusion](../images/deform/opExtrusion.png) | ![opExtrusion slice 1](../images/deform/opExtrusion_slice.png) | ![opExtrusion slice 2](../images/deform/opExtrusion_slice2.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `deform` |
| **GLSL File** | [`glsl/deform/opExtrusion.glsl`](../../glsl/deform/opExtrusion.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`) |

---

## Mathematical Overview

Extrudes any exact 2D SDF `sdf2d` along the $z$-axis to half-height $h$ with flat caps.

---

## GLSL Implementation

```glsl
float opExtrusion( in vec3 p, in float sdf2d, in float h )
{
    vec2 w = vec2( sdf2d, abs(p.z) - h );
    return min(max(w.x, w.y), 0.0) + length(max(w, 0.0));
}
```
