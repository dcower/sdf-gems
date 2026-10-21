# `opElongate` — Exact 3D Domain Elongation

| 3D Surface ($f(p) = 0$) | Primary 2D Cross-Section | Secondary 2D Cross-Section |
| :---: | :---: | :---: |
| ![opElongate](../images/deform/opElongate.png) | ![opElongate slice 1](../images/deform/opElongate_slice.png) | ![opElongate slice 2](../images/deform/opElongate_slice2.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `deform` |
| **GLSL File** | [`glsl/deform/opElongate.glsl`](../../glsl/deform/opElongate.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`) |

---

## Mathematical Overview

Splits a 3D primitive into 8 octants separated by an axis-aligned box of half-size $h$.

---

## GLSL Implementation

```glsl
vec4 opElongate( in vec3 p, in vec3 h )
{
    vec3 q = abs(p) - h;
    return vec4( max(q, 0.0), min(max(q.x, max(q.y, q.z)), 0.0) );
}
```
