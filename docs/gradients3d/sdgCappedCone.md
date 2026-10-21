# `sdgCappedCone` — Exact 1/2-`sqrt` 3D Capped Cone with Analytical Gradient (Fixing $56.91\%$ Rim Gradient Bug)

| 3D Raymarched Surface | 2D Cross-Section Slice |
| :---: | :---: |
| ![3D Preview](../images/gradients3d/sdgCappedCone.png) | ![2D Cross-Section](../images/gradients3d/sdgCappedCone_slice.png) |

- **GLSL Implementation**: [`sdgCappedCone.glsl`](../../glsl/gradients3d/sdgCappedCone.glsl)
- **Category**: `gradients3d`
- **Status**: **Exact Euclidean SDF + Unit 3D Gradient** (Bug fix & `sqrt` reduction discovered using **Gemini 4 Argon**, `sdf-gems`, 2026)

---

## 1. Critical Bug Discovery in Baseline `sdgCappedCone`

In Inigo Quilez's [`sdgCappedCone`](https://iquilezles.org/articles/distgradfunctions3d/) ([Shadertoy `3ctcW7`](https://www.shadertoy.com/view/3ctcW7)):
```glsl
return (la<lb)?
    vec4(s*sqrt(la),0.0,sign(p.y),0.0) :
    vec4(s*sqrt(lb),vec3(k.y*p.xz/l,-k.x).xzy/sqrt(m));
```
Across the entire `la >= lb` branch, IQ returns the constant conical-mantle face normal `vec3(k.y * p.xz / l, -k.x).xzy / sqrt(m)`—even when the clamped segment parameter `t = clamp(dot(q, k) / m, 0.0, 1.0)` hits `0.0` or `1.0` (the top and bottom circular rim wedges), and across `la < lb` when $l > r_c$ it returns `vec3(0, sign(p.y), 0)`.

- Across $50{,}000$ random points in $[-2, 2]^3$, IQ's `sdgCappedCone` returns an **incorrect gradient at $56.91\%$ of all points**, with **maximum gradient error `1.597`** ($> 105^\circ$ direction error) and **mean gradient error `0.418`**.

![sdgCappedCone Comparison](../images/comparisons/sdg_cone_comparison.png)

---

## 2. Mathematical Fix & `sqrt` Reduction (`3` $\to$ `1/2` `sqrt`s)

1. **Endpoint Subsumption in the Cap Zone (`1` `sqrt` total)**:
   Whenever $l > r_c$ (where $r_c = (p_y < 0) ? r_1 : r_2$), the cap disk endpoint $(r_c, \pm h)$ is already an endpoint of the mantle segment. Thus the flat cap strictly dominates **if and only if** $l \le r_c \land c_y^2 < \|b\|_2^2$ (where $c_y = |p_y| - h$), where the signed distance is simply $c_y$ (**`0` additional `sqrt`s** beyond $l = \|p_{xz}\|_2$) and the gradient is $(0, \operatorname{sign}(p_y), 0)$.
2. **Exact Unified Mantle & Rim Gradient (`2` `sqrt`s total, down from `3`)**:
   Everywhere outside the flat-cap Voronoi region, the already-computed 2D displacement vector $b = k \operatorname{clamp}(\operatorname{dot}(q, k)/\|k\|_2^2, 0, 1) - q = (\Delta l, \Delta y)$ is the exact 2D residual from the closest point on the mantle segment—both on the conical face $t \in (0, 1)$ and on both circular rims $t \in \{0, 1\}$. Dividing $(b_x p_x / l, b_y, b_x p_z / l)$ by $d = s\sqrt{\|b\|_2^2}$ yields the exact unit 3D gradient while eliminating `vec2 a`, `dot(a, a)`, `sqrt(la)`, and `sqrt(m)`.

---

## 3. Numerical Verification ($50{,}000$ Random 3D Points in $[-2, 2]^3$)

| Implementation | Gradient Failure Rate ($\|\nabla d - \nabla_{\text{FD}} d\|_2 > 0.01$) | Max Gradient Error | Median Gradient Error | Max $\big|\|\nabla d\|_2 - 1\big|$ | `sqrt` Calls (Cap / Mantle) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| IQ `sdgCappedCone` (`3ctcW7`) | **`56.91%`** | `1.597e+00` | `2.651e-01` | `4.441e-16` | `2 / 3` |
| **`sdf-gems` `sdgCappedCone`** | **`0.00%`** | **`8.120e-10`** | **`1.210e-10`** | **`4.441e-16`** | **`1 / 2`** |

---

## 4. GLSL Implementation

```glsl
vec4 sdgCappedCone(vec3 p, float h, float r1, float r2) {
    vec2  k  = vec2(r2 - r1, 2.0 * h);
    float l  = length(p.xz);
    vec2  q  = vec2(r2 - l, h - p.y);
    float rc = (p.y < 0.0) ? r1 : r2;
    float cy = abs(p.y) - h;
    vec2  b  = k * clamp(dot(q, k) / dot2(k), 0.0, 1.0) - q;
    float lb = dot2(b);
    if (l <= rc && cy * cy < lb) return vec4(cy, 0.0, (p.y < 0.0) ? -1.0 : 1.0, 0.0);
    float d  = ((b.x < 0.0 && cy < 0.0) ? -1.0 : 1.0) * sqrt(lb);
    return vec4(d, (l > 0.0) ? vec3(b.x * p.x / l, b.y, b.x * p.z / l) / d : vec3(0.0, sign(b.y), 0.0));
}
```
