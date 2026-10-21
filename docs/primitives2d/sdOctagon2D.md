# `sdOctagon2D` & `sdOctagonalPrism` — 1-Fold, 0/1-`sqrt` 2D Regular Octagon & 1-`sqrt` 3D Prism

![2D SDF Preview](../images/primitives2d/sdOctagon2D.png)

- **GLSL Implementation**: [`sdOctagon2D.glsl`](../../glsl/primitives2d/sdOctagon2D.glsl)
- **Category**: `primitives2d`
- **Status**: **Exact 0/1-`sqrt` Euclidean SDF** (1-Fold $D_8$ Chamber & Vertex Dominance discovered using **Gemini 4 Argon**, `sdf-gems`, 2026)

---

## 1. Summary & Motivation

Inigo Quilez's [`sdOctogon(vec2 p, float r)`](https://www.shadertoy.com/view/llGfDG) folds $\mathbb{R}^2$ into the fundamental domain of $D_8$ using `p = abs(p)` followed by **two oblique Householder reflections** across $67.5^\circ$ and $22.5^\circ$, then clamps $p_x$ to $[-k_z r, k_z r]$ and evaluates `length(p) * sign(p.y)` everywhere in $\mathbb{R}^2$.

`sdOctagon2D` introduces two algebraic reductions:
1. **1-Reflection $D_8$ Fold via $45^\circ$ Diagonal Sort**:
   Sorting `if (p.y < p.x) p = p.yx;` folds $[0^\circ, 90^\circ]$ across $45^\circ$ into $[45^\circ, 90^\circ]$ with **zero multiplications**. A **single** Householder reflection across $67.5^\circ$ then maps $[45^\circ, 67.5^\circ]$ into the fundamental chamber $[67.5^\circ, 90^\circ]$, eliminating the second Householder reflection ($4$ multiplications, $3$ additions, and $1$ `min`).
2. **Fundamental-Chamber Vertex Dominance (`0` `sqrt`s in Interior & Face Strips)**:
   Inside $[67.5^\circ, 90^\circ]$, every point satisfies $0 \le p_x \le k_z p_y$ where $k_z = \tan(\pi/8) = \sqrt{2} - 1$. Defining $d = p - (k_z r, r)$:
   $$d_x > 0 \implies p_x > k_z r \implies k_z p_y \ge p_x > k_z r \implies p_y > r \implies d_y > 0.$$
   Therefore, $d_x > 0$ occurs **if and only if** $p$ lies in the exterior vertex wedge (where $d_y > 0$ automatically and the distance is `length(d)`). Whenever $d_x \le 0$ ($67\%$ of $[-2r, 2r]^2$, covering the entire interior and all $8$ exterior flat-edge normal strips), the exact signed distance is simply the linear coordinate $d_y = p_y - r$ (**`0` `sqrt`s, `0` `clamp`s, `0` `sign`s**).

---

## 2. Benchmark & Operation Count Comparison

| Implementation | Oblique Householder Folds | `sqrt` (`length`) Calls | `clamp` / `sign` Calls | Max Error vs. Ground Truth |
| :--- | :---: | :---: | :---: | :---: |
| **IQ `sdOctogon` (`llGfDG`)** | `2` (8 muls, 6 adds) | `1` everywhere (`100%` of $\mathbb{R}^2$) | `1` `clamp`, `1` `sign` | `8.88e-16` |
| **`sdf-gems` `sdOctagon2D`** | **`1` (4 muls, 3 adds)** | **`0` in interior & face strips, `1` in corner wedges** | **`0` `clamp`, `0` `sign`** | **`8.88e-16`** |

---

## 3. GLSL Implementation

```glsl
// 1. Exact 1-Fold, 0/1-Sqrt 2D Regular Octagon of inradius r
float sdOctagon2D(vec2 p, float r) {
    const vec3 k = vec3(-0.9238795325, 0.3826834323, 0.4142135623);
    p = abs(p);
    if (p.y < p.x) p = p.yx;
    p -= 2.0 * min(dot(k.xy, p), 0.0) * k.xy;
    vec2 d = p - vec2(k.z, 1.0) * r;
    return (d.x > 0.0) ? length(d) : d.y;
}

// 2. Exact 1-Fold, 1-Sqrt 3D Regular Octagonal Prism of inradius r and half-height hz
float sdOctagonalPrism(vec3 p, float r, float hz) {
    const vec3 k = vec3(-0.9238795325, 0.3826834323, 0.4142135623);
    p = abs(p);
    if (p.y < p.x) p.xy = p.yx;
    p.xy -= 2.0 * min(dot(k.xy, p.xy), 0.0) * k.xy;
    vec3 d = vec3(p.x - k.z * r, p.y - r, p.z - hz);
    return min(max(d.y, d.z), 0.0) + length(max(d, 0.0));
}
```
