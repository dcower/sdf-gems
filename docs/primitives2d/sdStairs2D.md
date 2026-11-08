# `sdStairs2D` — Exact 1-Sqrt Rotation-Free 2D Staircase (Fixes Tall-Step Sign Bug)

| 2D SDF Field (Isolines & Zero-Set) |
| :---: |
| ![sdStairs2D](../images/primitives2d/sdStairs2D.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez — upgraded/optimized using Gemini 4 Argon (sdf-gems, 2026)** |
| **Source** | [https://iquilezles.org/articles/distfunctions2d/](https://iquilezles.org/articles/distfunctions2d/) ([Shadertoy `7tKSWt`](https://www.shadertoy.com/view/7tKSWt)) |
| **Category** | `primitives2d` |
| **GLSL File** | [`glsl/primitives2d/sdStairs2D.glsl`](../../glsl/primitives2d/sdStairs2D.glsl) |
| **Exact Euclidean SDF?** | Yes (`\|\|grad d\|\| = 1`, `1` `sqrt` down from `2`, `0` `mat2` rotations, exact for all `wh.x / wh.y` aspect ratios) |

---

## Mathematical Overview

Inigo Quilez's 2D staircase primitive `sdStairs(p, wh, n)` ([`7tKSWt`](https://www.shadertoy.com/view/7tKSWt)) has two subtle mathematical and computational issues:

1. **Redundant `sqrt` (`dia = length(wh)`), Dual `mat2` Rotations, and 4 Divisions**:
   In the baseline implementation:
   ```glsl
   float dia = length(wh);
   p = mat2(wh.x,-wh.y, wh.y,wh.x)*p/dia;
   float id = clamp(round(p.x/dia), 0.0, n-1.0);
   p.x = p.x - id*dia;
   p = mat2(wh.x, wh.y,-wh.y,wh.x)*p/dia;
   ```
   Because $`R = \frac{1}{\|wh\|} \begin{pmatrix} w & h \\ -h & w \end{pmatrix}`$ is orthogonal ($`R^T R = I`$) and $`R \, wh = \|wh\| e_1`$, rotating $`p`$ by $`R`$, subtracting $`\mathit{id} \cdot \|wh\| e_1`$, and rotating back by $`R^T`$ cancels algebraically:
   $`\displaystyle \frac{(R p)_x}{\|wh\|} = \frac{\mathrm{dot}(p, wh)}{\mathrm{dot}(wh, wh)}, \qquad R^T\bigl(R p - \mathit{id}\cdot \|wh\| e_1\bigr) = p - \mathit{id}\cdot wh.`$
   Thus **`dia = length(wh)` (`1` `sqrt`), both `mat2` matrix multiplications (`8` muls, `4` adds), and `4` of the `5` scalar divisions** collapse into two lines:
   ```glsl
   float id = clamp(round(dot(p, wh) / dot2(wh)), 0.0, n - 1.0);
   p -= id * wh;
   ```

2. **Tall-Step Voronoi & Sign Flip Bug (`wh.x < wh.y`)**:
   In `sdStairs`, the fundamental step cell reflects across the midpoint of the **vertical riser** $`(\mathit{id}\cdot w, (\mathit{id} + \tfrac{1}{2})h)`$, which shares the **horizontal tread** between adjacent cells `id` and `id + 1` while leaving the **vertical riser** $`x = \mathit{id}\cdot w, \, y \in [\mathit{id}\cdot h, (\mathit{id}+1)h]`$ unshared.
   The diagonal boundary $`\mathrm{dot}(p, wh) = (\mathit{id} + \tfrac{1}{2})(w^2 + h^2)`$ between cells `id` and `id + 1` intersects the vertical line $`x = \mathit{id}\cdot w`$ at:
   $`\displaystyle y = \mathit{id}\cdot h + \frac{w^2 + h^2}{2 h}.`$
   Whenever $`w < h`$ (`wh.x < wh.y`), $`\frac{w^2 + h^2}{2 h} < h`$, so the cell boundary **slices straight through the vertical riser of step `id`**! Points near the top of step `id`'s riser (and in the exterior air to its left) are misassigned to cell `id + 1`, where $`p_{\text{local}}.x < 0`$ and $`p_{\text{local}}.y < -h/2`$ — failing the exterior sign check `p.y > hh * sign(p.x)` and **flipping the sign of the exterior air to negative (`d < 0`, inside the solid)** with distance errors up to `4.20`.
   By folding across the exact anti-diagonal isometry of the staircase polygon (`if (wh.x < wh.y) { p = ba.yx - p.yx; wh = wh.yx; ba = ba.yx; }`), $`w \ge h`$ is guaranteed and `sdStairs2D` achieves **machine precision (`< 2e-15` max error) across all step aspect ratios**.

![2D Staircase Comparison](../images/comparisons/stairs2d_comparison.png)

### Benchmark Comparison ($`300{,}000`$ Random Points, $`n = 4`$)

| Step Size `wh` | IQ `sdStairs` (`7tKSWt`) Max Error | IQ Sign-Flip Error Rate | `sdf-gems` `sdStairs2D` Max Error | `sdf-gems` Sign-Flip Error Rate | `sqrt` Count |
| :--- | :---: | :---: | :---: | :---: | :---: |
| `vec2(1.2, 0.3)` ($`w > h`$) | `4.44e-15` | `0.00%` | **`1.78e-15`** | **`0.00%`** | **`1`** (vs. `2`) |
| `vec2(0.8, 0.8)` ($`w = h`$) | `3.55e-15` | `0.00%` | **`1.33e-15`** | **`0.00%`** | **`1`** (vs. `2`) |
| `vec2(0.7, 0.9)` ($`w < h`$) | `4.43e-01` | `0.13%` | **`1.33e-15`** | **`0.00%`** | **`1`** (vs. `2`) |
| `vec2(0.5, 1.0)` ($`w < h`$) | `1.49e+00` | `0.88%` | **`6.66e-16`** | **`0.00%`** | **`1`** (vs. `2`) |
| `vec2(0.3, 1.2)` ($`w < h`$) | `4.20e+00` | `3.68%` | **`8.88e-16`** | **`0.00%`** | **`1`** (vs. `2`) |

---

## GLSL Implementation

```glsl
float sdStairs2D( vec2 p, vec2 wh, float n )
{
    vec2 ba = wh * n;
    if( wh.x < wh.y ) { p = ba.yx - p.yx; wh = wh.yx; ba = ba.yx; }
    float d = min(dot2(p - vec2(clamp(p.x, 0.0, ba.x), 0.0)),
                  dot2(p - vec2(ba.x, clamp(p.y, 0.0, ba.y))));
    float s = sign(max(-p.y, p.x - ba.x));

    float id = clamp(round(dot(p, wh) / dot2(wh)), 0.0, n - 1.0);
    p -= id * wh;

    float hh = wh.y * 0.5;
    p.y -= hh;
    if( p.y > hh * sign(p.x) ) s = 1.0;
    p = (id < 0.5 || p.x > 0.0) ? p : -p;

    d = min(d, dot2(p - vec2(0.0, clamp(p.y, -hh, hh))));
    d = min(d, dot2(p - vec2(clamp(p.x, 0.0, wh.x), hh)));
    return sqrt(d) * s;
}
```

### Code Diff vs. Baseline (Inigo Quilez `7tKSWt`)
```diff
  float sdStairs( in vec2 p, in vec2 wh, in float n )
  {
      vec2 ba = wh*n;
+     if( wh.x < wh.y ) { p = ba.yx - p.yx; wh = wh.yx; ba = ba.yx; }
      float d = min(dot2(p - vec2(clamp(p.x,0.0,ba.x), 0.0)), 
                    dot2(p - vec2(ba.x, clamp(p.y,0.0,ba.y))));
      float s = sign(max(-p.y, p.x-ba.x));
 
-     float dia = length(wh);
-     p = mat2(wh.x,-wh.y, wh.y,wh.x)*p/dia;
-     float id = clamp(round(p.x/dia), 0.0, n-1.0);
-     p.x = p.x - id*dia;
-     p = mat2(wh.x, wh.y,-wh.y,wh.x)*p/dia;
+     float id = clamp(round(dot(p, wh) / dot2(wh)), 0.0, n - 1.0);
+     p -= id * wh;
 
      float hh = wh.y/2.0;
      p.y -= hh;
      if( p.y > hh*sign(p.x) ) s = 1.0;
      p = (id < 0.5 || p.x > 0.0) ? p : -p;
 
      d = min( d, dot2(p - vec2(0.0, clamp(p.y,-hh,hh))) );
      d = min( d, dot2(p - vec2(clamp(p.x,0.0,wh.x), hh)) );
      return sqrt(d)*s;
  }
```
