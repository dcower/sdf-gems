# `sminCubic` — Cubic Polynomial Smooth Minimum

![sminCubic](../images/booleans/sminCubic.png)

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** |
| **Source** | [https://iquilezles.org/articles/smin/](https://iquilezles.org/articles/smin/) |
| **Category** | `booleans` |
| **GLSL File** | [`glsl/booleans/sminCubic.glsl`](../../glsl/booleans/sminCubic.glsl) |
| **Exact Euclidean SDF?** | Conservative (`||grad|| <= 1`) |

---

## Mathematical Overview

$C^2$ cubic Closest-Distance (CD) smooth minimum.

---

## GLSL Implementation

```glsl
float sminCubic( float a, float b, float k )
{
    k *= 6.0;
    float h = max( k - abs(a - b), 0.0 )/k;
    return min(a, b) - h*h*h*k*(1.0/6.0);
}
```
