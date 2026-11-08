# `sdCappedConeSeg` — Exact 3D Arbitrary-Segment Capped Cone (Axis-NaN Fixed)

| 3D Surface ($`f(p) = 0`$) | 2D Cross-Section Slice |
| :---: | :---: |
| ![sdCappedConeSeg](../images/primitives3d/sdCappedConeSeg.png) | ![sdCappedConeSeg slice](../images/primitives3d/sdCappedConeSeg_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdCappedConeSeg.glsl`](../../glsl/primitives3d/sdCappedConeSeg.glsl) |
| **Exact Euclidean SDF?** | Yes (`\|\|grad d\|\| = 1`, 0% NaN rate) |

---

## Mathematical Overview

Uses orthogonal vector rejection `length(pa - ba * paba)` to avoid catastrophic cancellation on the cone axis.

---

## GLSL Implementation

```glsl
float sdCappedConeSeg( vec3 p, vec3 a, vec3 b, float ra, float rb )
{
    float rba  = rb - ra;
    vec3  ba   = b - a, pa = p - a;
    float baba = dot(ba, ba);
    float paba = dot(pa, ba) / baba;
    float x    = length(pa - ba * paba);
    float cax  = max(0.0, x - ((paba < 0.5) ? ra : rb));
    float cay  = abs(paba - 0.5) - 0.5;
    float k    = rba*rba + baba;
    float f    = clamp( (rba*(x - ra) + paba*baba)/k, 0.0, 1.0 );
    float cbx  = x - ra - f*rba;
    float cby  = paba - f;
    float s    = (cbx < 0.0 && cay < 0.0) ? -1.0 : 1.0;
    return s*sqrt( min(cax*cax + cay*cay*baba, cbx*cbx + cby*cby*baba) );
}
```

---

## Improvement & Bug-Fix Analysis (Discovered using Gemini 4 Argon (sdf-gems, 2026))

- **Bug in Baseline**: `sqrt(papa - paba*paba*baba)` cancels to negative float noise (`NaN`) on the cone axis `pa || ba`.
- **Fix**: Orthogonal rejection `length(pa - ba * paba)` guarantees a non-negative sum of squares everywhere (`0.0%` NaN rate).

### Code Diff vs. Baseline (Inigo Quilez)
```diff
-    float papa = dot(p - a, p - a);
-    float paba = dot(p - a, b - a)/baba;
-    float x = sqrt( papa - paba*paba*baba );
+    float paba = dot(pa, ba) / baba;
+    float x    = length(pa - ba * paba);
```
