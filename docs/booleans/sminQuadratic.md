# `sminQuadratic` — Quadratic Polynomial Smooth Minimum

![sminQuadratic](../images/booleans/sminQuadratic.png)

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** |
| **Source** | [https://iquilezles.org/articles/smin/](https://iquilezles.org/articles/smin/) |
| **Category** | `booleans` |
| **GLSL File** | [`glsl/booleans/sminQuadratic.glsl`](../../glsl/booleans/sminQuadratic.glsl) |
| **Exact Euclidean SDF?** | Conservative (`||grad|| <= 1`) |

---

## Mathematical Overview

$C^1$ quadratic Closest-Distance (CD) smooth minimum.

---

## GLSL Implementation

```glsl
float sminQuadratic( float a, float b, float k )
{
    k *= 4.0;
    float h = max( k - abs(a - b), 0.0 )/k;
    return min(a, b) - h*h*k*(1.0/4.0);
}
```
