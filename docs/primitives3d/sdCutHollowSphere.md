# `sdCutHollowSphere` — Exact 1/2-`sqrt` 3D Cut Hollow Spherical Shell (Early-Return Shell Branch)

| 3D Surface ($f(p) = 0$) | 2D Cross-Section Slice |
| :---: | :---: |
| ![sdCutHollowSphere](../images/primitives3d/sdCutHollowSphere.png) | ![sdCutHollowSphere slice](../images/primitives3d/sdCutHollowSphere_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez — upgraded/optimized using Gemini 4 Argon (`sdf-gems`, 2026)** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdCutHollowSphere.glsl`](../../glsl/primitives3d/sdCutHollowSphere.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`, `1` `sqrt` in spherical shell [$65.4\%$ of space] & `2` in rim, down from `3`, `0.0` error) |

---

## Mathematical Overview & Conjugate-Square Rationalization

Using an explicit early-return `if( h * abs(h) * x2 >= w2 * p.y * abs(p.y) ) return abs(sqrt(x2 + p.y * p.y) - r) - t;` ensures the GPU only executes **1 `sqrt`** across the $65.4\%$ spherical shell sector, and **2 `sqrt`s** in the circular rim wedge via the conjugate-square identity.

---

## GLSL Implementation

```glsl
float sdCutHollowSphere( vec3 p, float r, float h, float t )
{
    float x2 = dot(p.xz, p.xz), w2 = r * r - h * h;
    if( h * abs(h) * x2 >= w2 * p.y * abs(p.y) )
        return abs(sqrt(x2 + p.y * p.y) - r) - t;
    float d = x2 - w2, dy = p.y - h;
    return sqrt(d * d / max(x2 + w2 + 2.0 * sqrt(w2 * x2), 1e-30) + dy * dy) - t;
}
```
