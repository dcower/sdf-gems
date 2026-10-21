# `opRep2D` — 2D Neighbor-Corrected Domain Repetition (4-Tap Cartesian)

| 3D Surface ($f(p) = 0$) | 2D Cross-Section Slice |
| :---: | :---: |
| ![opRep2D](../images/repetition/opRep2D.png) | ![opRep2D slice](../images/repetition/opRep2D_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** |
| **Source** | [https://iquilezles.org/articles/sdfrepetition/](https://iquilezles.org/articles/sdfrepetition/) |
| **Category** | `repetition` |
| **GLSL File** | [`glsl/repetition/opRep2D.glsl`](../../glsl/repetition/opRep2D.glsl) |
| **Exact Euclidean SDF?** | Yes |

---

## Mathematical Overview

Evaluates the 4 closest 2D Cartesian grid cells surrounding $p$.

---

## GLSL Implementation

```glsl
void opRep2D_4Tap( in vec2 p, in float s, out vec2 c0, out vec2 c1, out vec2 c2, out vec2 c3 )
{
    vec2 q = p / s;
    vec2 id = round(q);
    vec2 o = sign(q - id);
    c0 = id + vec2(0.0, 0.0);
    c1 = id + vec2(o.x, 0.0);
    c2 = id + vec2(0.0, o.y);
    c3 = id + vec2(o.x, o.y);
}
```
