# `sdCappedCylinderSeg` — Exact 3D Arbitrary-Segment Capped Cylinder (FP16 Overflow Fixed)

| 3D Surface ($f(p) = 0$) | 2D Cross-Section Slice |
| :---: | :---: |
| ![sdCappedCylinderSeg](../images/primitives3d/sdCappedCylinderSeg.png) | ![sdCappedCylinderSeg slice](../images/primitives3d/sdCappedCylinderSeg_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdCappedCylinderSeg.glsl`](../../glsl/primitives3d/sdCappedCylinderSeg.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`) |

---

## Mathematical Overview

Normalizes the segment axis $n = (b-a)/\|b-a\|$ first so all intermediate quantities scale linearly as $O(\|b-a\|^1)$ instead of $O(\|b-a\|^6)$.

---

## GLSL Implementation

```glsl
float sdCappedCylinderSeg( vec3 p, vec3 a, vec3 b, float r )
{
    vec3  ba = b - a, pa = p - a;
    float L  = length(ba);
    vec3  n  = ba / L;
    float y  = dot(pa, n);
    float x  = length(pa - n * y) - r;
    vec2  d  = vec2(x, abs(y - 0.5 * L) - 0.5 * L);
    return min(max(d.x, d.y), 0.0) + length(max(d, 0.0));
}
```

---

## Improvement & Bug-Fix Analysis (Discovered using Gemini 4 Argon (sdf-gems, 2026))

- **Bug in Baseline**: Scaling intermediates by `baba = |b-a|^2` causes `x2 = x*x` to scale as `O(|b-a|^6)`, overflowing FP16 (`mediump`, max 65504) whenever `|b - a| > 11.5`.
- **Fix**: Normalizing `ba` first reduces intermediate scaling to `O(|b-a|^1)` and shortens the code by 4 lines.

### Code Diff vs. Baseline (Inigo Quilez)
```diff
-    float baba = dot(ba,ba);
-    float paba = dot(pa,ba);
-    float x = length(pa*baba - ba*paba) - r*baba;
-    float y = abs(paba - baba*0.5) - baba*0.5;
-    float x2 = x*x;
-    float y2 = y*y*baba;
-    float d = (max(x,y) < 0.0) ? -min(x2,y2) : (((x > 0.0) ? x2 : 0.0) + ((y > 0.0) ? y2 : 0.0));
-    return sign(d)*sqrt(abs(d))/baba;
+    float L  = length(ba);
+    vec3  n  = ba / L;
+    float y  = dot(pa, n);
+    float x  = length(pa - n * y) - r;
+    vec2  d  = vec2(x, abs(y - 0.5 * L) - 0.5 * L);
+    return min(max(d.x, d.y), 0.0) + length(max(d, 0.0));
```
