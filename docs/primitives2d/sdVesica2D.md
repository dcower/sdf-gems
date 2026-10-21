# `sdVesica2D` — Exact 1-to-2 Sqrt 2D Vesica Piscis (Signed-Square Cone Reduction)

| 2D SDF Field (Isolines & Zero-Set) |
| :---: |
| ![sdVesica2D](../images/primitives2d/sdVesica2D.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez — upgraded/optimized using Gemini 4 Argon (sdf-gems, 2026)** |
| **Source** | [https://iquilezles.org/articles/distfunctions2d/](https://iquilezles.org/articles/distfunctions2d/) |
| **Category** | `primitives2d` |
| **GLSL File** | [`glsl/primitives2d/sdVesica2D.glsl`](../../glsl/primitives2d/sdVesica2D.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`, `1` `sqrt` in main body & `1` `sqrt` everywhere for width-parameterized `sdVesicaWidth2D`, `0.0` error) |

---

## Mathematical Overview

Inigo Quilez's 2D `sdVesica(p, r, d)` computes `b = sqrt(r*r - d*d)` upfront (`sqrt` #1) solely to test the tip-cone condition `(p.y - b)*d > p.x*b`, even though the main body branch `length(p - vec2(-d, 0.0)) - r` (`sqrt` #2) never uses `b`!

Rewriting the tip-cone condition in the first quadrant $p = |p| \ge 0$ (with $0 \le d < r$) as:
$$d \, p_y > b(p_x + d) \iff d^2 p_y^2 > (r^2 - d^2)(p_x + d)^2$$
eliminates `sqrt(r*r - d*d)` from the branch condition and from the entire main-body region (`0.0` bitwise error vs. baseline). Furthermore, when parameterized by physical half-length $r$ and half-width $w$ (`sdVesicaWidth2D(p, r, w)`), $d = \frac{r^2 - w^2}{2 w}$ is rational and `sdVesicaWidth2D` executes **only 1 `sqrt` everywhere in $\mathbb{R}^2$**.

![Vesica 1-to-2 Sqrt Reduction](../images/comparisons/vesica_comparison.png)

---

## GLSL Implementation

```glsl
float sdVesica2D( vec2 p, float r, float d )
{
    p = abs(p);
    float b2 = r * r - d * d;
    float x  = p.x + d;
    return (d * d * p.y * p.y > b2 * x * x)
         ? length(vec2(p.x, p.y - sqrt(b2)))
         : length(vec2(x, p.y)) - r;
}

// Physical half-length (r) & half-width (w) parameterization: 1 `sqrt` everywhere in R^2!
float sdVesicaWidth2D( vec2 p, float r, float w )
{
    p = abs(p);
    float d = 0.5 * (r * r - w * w) / w;
    return (r * p.x < d * (p.y - r))
         ? length(vec2(p.x, p.y - r))
         : length(vec2(p.x + d, p.y)) - (d + w);
}
```

### Code Diff vs. Baseline (Inigo Quilez)
```diff
  float sdVesica( vec2 p, float r, float d )
  {
      p = abs(p);
-     float b = sqrt(r*r-d*d);
-     return ((p.y-b)*d > p.x*b) ? length(p-vec2(0.0,b))
-                                : length(p-vec2(-d,0.0))-r;
+     float b2 = r * r - d * d;
+     float x  = p.x + d;
+     return (d * d * p.y * p.y > b2 * x * x)
+          ? length(vec2(p.x, p.y - sqrt(b2)))
+          : length(vec2(x, p.y)) - r;
  }
```
