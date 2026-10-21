# `opOnion` — Exact Hollow Shell / Onion Operator

| 3D Surface ($f(p) = 0$) | 2D Cross-Section Slice |
| :---: | :---: |
| ![opOnion](../images/deform/opOnion.png) | ![opOnion slice](../images/deform/opOnion_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `deform` |
| **GLSL File** | [`glsl/deform/opOnion.glsl`](../../glsl/deform/opOnion.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`) |

---

## Mathematical Overview

Turns any solid SDF into a thin-walled shell of half-thickness `thickness` via $|d| - t$.

---

## GLSL Implementation

```glsl
float opOnion( in float d, in float thickness ) { return abs(d) - thickness; }
```
