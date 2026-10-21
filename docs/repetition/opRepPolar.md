# `opRepPolar` — Polar / Rotational Domain Repetition

| 3D Surface ($f(p) = 0$) | 2D Cross-Section Slice |
| :---: | :---: |
| ![opRepPolar](../images/repetition/opRepPolar.png) | ![opRepPolar slice](../images/repetition/opRepPolar_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Mercury Demogroup** |
| **Source** | [https://mercury.sexy/hg_sdf/](https://mercury.sexy/hg_sdf/) |
| **Category** | `repetition` |
| **GLSL File** | [`glsl/repetition/opRepPolar.glsl`](../../glsl/repetition/opRepPolar.glsl) |
| **Exact Euclidean SDF?** | Yes (symmetric SDFs) |

---

## Mathematical Overview

Repeats 2D space angularly into $n$ sectors around the origin.

---

## GLSL Implementation

```glsl
vec2 opRepPolar( in vec2 p, in float repetitions )
{
    float angle = SDF_TAU / repetitions;
    float a = atan(p.y, p.x) + angle*0.5;
    float r = length(p);
    a = mod(a, angle) - angle*0.5;
    return vec2(cos(a), sin(a))*r;
}
```
