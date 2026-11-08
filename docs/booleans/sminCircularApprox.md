# `sminCircularApprox` — [3/1] Padé Rational Circular CD Smooth Minimum

![sminCircularApprox](../images/booleans/sminCircularApprox.png)

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)** |
| **Source** | [https://iquilezles.org/articles/smin/](https://iquilezles.org/articles/smin/) |
| **Category** | `booleans` |
| **GLSL File** | [`glsl/booleans/sminCircularApprox.glsl`](../../glsl/booleans/sminCircularApprox.glsl) |
| **Exact Euclidean SDF?** | Conservative (`9.80e-4` max error — `11.0x` more accurate than IQ's quartic) |

---

## Mathematical Overview

Matches the exact endpoint values and derivatives of $`1 - \sqrt{1 - h(2-h)}`$ at $`h=0`$ and $`h=1`$ using a `[3/1]` rational Padé approximant with $`\sqrt{2}`$ coefficients.

---

## GLSL Implementation

```glsl
float sminCircularApprox( float a, float b, float k )
{
    k *= 1.0/(1.0 - sqrt(0.5));
    float h = max( k - abs(a - b), 0.0 )/k;
    return min(a, b) - k * h * h * (1.0 + 0.41421356 * h) / (2.0 * (1.0 + 1.41421356 * h));
}
```

---

## Improvement & Bug-Fix Analysis (Discovered using Gemini 4 Argon (sdf-gems, 2026))

- **Before (IQ Quartic Polynomial)**: Max error `1.08e-2`.
- **After (`sdf-gems` `[3/1]` Padé Rational)**: Max error **`9.80e-4` (`11.0x` more accurate)** with 1 fewer multiply.

### Code Diff vs. Baseline (Inigo Quilez)
```diff
-    float x = max( k - abs(a - b), 0.0 )/k;
-    return min(a, b) - k*0.5*x*x*(x*(x*(1.0 - sqrt(0.5)) - (2.0 - sqrt(2.0))) + 1.0);
+    float h = max( k - abs(a - b), 0.0 )/k;
+    return min(a, b) - k * h * h * (1.0 + 0.41421356 * h) / (2.0 * (1.0 + 1.41421356 * h));
```
