# `opRep` — 1-Tap Unlimited & Limited Cartesian Repetition

| 3D Surface ($f(p) = 0$) | 2D Cross-Section Slice |
| :---: | :---: |
| ![opRep](../images/repetition/opRep.png) | ![opRep slice](../images/repetition/opRep_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** |
| **Source** | [https://iquilezles.org/articles/sdfrepetition/](https://iquilezles.org/articles/sdfrepetition/) |
| **Category** | `repetition` |
| **GLSL File** | [`glsl/repetition/opRep.glsl`](../../glsl/repetition/opRep.glsl) |
| **Exact Euclidean SDF?** | Yes (symmetric SDFs) |

---

## Mathematical Overview

Folds space into periodic or clamped cubic cells via `round(p/s)`.

---

## GLSL Implementation

```glsl
vec3 opRep( in vec3 p, in vec3 s )
{
    return p - s*round(p/s);
}

vec3 opRepLim( in vec3 p, in float s, in vec3 l )
{
    return p - s*clamp(round(p/s), -l, l);
}
```
