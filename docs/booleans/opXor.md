# `opXor` — Exact Boolean Symmetric Difference (XOR)

![opXor](../images/booleans/opXor.png)

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** |
| **Source** | [https://iquilezles.org/articles/sdfxor/](https://iquilezles.org/articles/sdfxor/) |
| **Category** | `booleans` |
| **GLSL File** | [`glsl/booleans/opXor.glsl`](../../glsl/booleans/opXor.glsl) |
| **Exact Euclidean SDF?** | Yes (up to medial axis) |

---

## Mathematical Overview

Combines $\min(d_1, d_2)$ and $-\max(d_1, d_2)$ to preserve exact signed distances across symmetric difference boundaries.

---

## GLSL Implementation

```glsl
float opXor( float d1, float d2 ) { return max(min(d1, d2), -max(d1, d2)); }
```
