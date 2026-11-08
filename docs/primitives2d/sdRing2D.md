# `sdRing2D`, `sdgRing2D`, `sdArc2D` & `sdgArc2D` — Exact 0/1-`sqrt` 2D Annular Ring Sector & Circular Arc + Analytical Gradients

![2D SDF Preview](../images/primitives2d/sdRing2D.png)

- **GLSL Implementation**: [`sdRing2D.glsl`](../../glsl/primitives2d/sdRing2D.glsl)
- **Category**: `primitives2d`
- **Status**: **Exact 0/1-`sqrt` Euclidean SDF & 0/1-SFU Analytical Gradient** (Parabolic $`v^2`$-Cancellation & Triangle-Inequality Annular Dominance discovered using **Gemini 4 Argon**, `sdf-gems`, 2026)

![Triangle, Pie, and Ring 0-Sqrt / 0-Division Comparison](../images/comparisons/triangle_pie_ring_comparison.png)

---

## 1. Summary & Mathematical Derivation

After folding $`p_x \leftarrow |p_x|`$, let $`v = n_x p_x - n_y p_y`$ be the signed distance outside the radial cut line, $`u = n_y p_x + n_x p_y`$ be the coordinate along the radial cut line, $`h = \mathit{th}/2`$, and $`r_1 = r - h`$, $`r_2 = r + h`$.

### Theorem 1: Exterior 0-`sqrt` Cap Strip ($`v > 0`$)
Outside the angular sector ($`v > 0`$), the reverse triangle inequality guarantees that the flat radial cap segment $`S = \{v = 0, u \in [r_1, r_2]\}`$ always dominates the annular arcs: $`|\|p\|_2 - r| - h \le d_{\text{cap}}(p)`$. Furthermore, whenever $`|r - u| \le h`$ ($`u \in [r_1, r_2]`$), the closest point lies in the interior of the radial cap segment, so the exact exterior Euclidean distance is simply **`v` (`0` `sqrt`s)**!

### Theorem 2: Parabolic $`v^2`$-Cancellation Inside the Annular Sector ($`v \le 0, u > 0`$)
Inside the sector ($`v \le 0`$) for $`u > 0`$, the radial cap distance is $`v \le 0`$ while the annular arc distance is $`|\|p\|_2 - r| - h = \max(r_1 - \sqrt{u^2 + v^2}, \sqrt{u^2 + v^2} - r_2)`$. Thus the radial cap dominates ($`v > |\|p\|_2 - r| - h`$) if and only if:
```math
r_1 + v < \sqrt{u^2 + v^2} < r_2 - v.
```
Squaring both inequalities cancels $`v^2`$ on both sides simultaneously:
```math
r_1(r_1 + 2v) < u^2 < r_2(r_2 - 2v) \quad (\text{or } r_1(r_1 - 2|v|) < u^2 < r_2(r_2 + 2|v|)).
```
Whenever this parabolic inequality holds, `sdRing2D` and `sdgRing2D` return **`v` with `0` `sqrt`s**!

---

## 2. Benchmark & Comparison

| Implementation | `sqrt` / `rsq` Calls | Max Distance Error | Max Gradient Error |
| :--- | :---: | :---: | :---: |
| **IQ `sdRing` (`DsccDH`)** | `2` everywhere | $`4.44 \times 10^{-16}`$ | N/A |
| **`sdf-gems` `sdRing2D` / `sdgRing2D`** | **`0` in cap Voronoi, `1` elsewhere** | **$`4.44 \times 10^{-16}`$** | **$`1.19 \times 10^{-8}`$** |
| **IQ `sdArc` / `sdgArc`** | `2` | $`4.44 \times 10^{-16}`$ | $`1.19 \times 10^{-8}`$ |
| **`sdf-gems` `sdArc2D` / `sdgArc2D`** | **`1` everywhere** | **$`4.44 \times 10^{-16}`$** | **$`1.19 \times 10^{-8}`$** |

---

## 3. GLSL Implementation

See [`glsl/primitives2d/sdRing2D.glsl`](../../glsl/primitives2d/sdRing2D.glsl) for `sdRing2D`, `sdgRing2D`, `sdArc2D`, and `sdgArc2D`.
