# `opUnion` — Hard Boolean Union

![opUnion](../images/booleans/opUnion.png)

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `booleans` |
| **GLSL File** | [`glsl/booleans/opUnion.glsl`](../../glsl/booleans/opUnion.glsl) |
| **Exact Euclidean SDF?** | Yes (exterior) |

---

## Mathematical Overview

Computes $`\min(d_1, d_2)`$, which is exact in the exterior of both shapes.

---

## GLSL Implementation

```glsl
float opUnion( float d1, float d2 ) { return min(d1, d2); }
```
