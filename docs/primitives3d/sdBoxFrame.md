# `sdBoxFrame` — Exact 1-Sqrt 3D Wireframe Box Frame (Closest-Beam Reduction)

| 3D Surface ($`f(p) = 0`$) | Primary 2D Cross-Section | Secondary 2D Cross-Section |
| :---: | :---: | :---: |
| ![sdBoxFrame](../images/primitives3d/sdBoxFrame.png) | ![sdBoxFrame slice 1](../images/primitives3d/sdBoxFrame_slice.png) | ![sdBoxFrame slice 2](../images/primitives3d/sdBoxFrame_slice2.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez — upgraded/optimized using Gemini 4 Argon (sdf-gems, 2026)** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdBoxFrame.glsl`](../../glsl/primitives3d/sdBoxFrame.glsl) |
| **Exact Euclidean SDF?** | Yes (`\|\|grad d\|\| = 1`, `1` `sqrt` down from `3`; exterior identical to the 3-beam form, interior corrected in the corner cubes, see below) |

---

## Mathematical Overview & Closest-Beam Reduction Theorem

Inigo Quilez's baseline `sdBoxFrame` evaluates **three separate 3D box SDFs** (one for each axis-aligned beam in the positive octant) and takes their minimum, requiring **3 `length()` (`sqrt`) calls** and **23 scalar `min`/`max` calls**.

We prove that **all three beams never need to be evaluated**: the closest beam at any point $`p_0 \in \mathbb{R}^3`$ is **identically the single beam along the coordinate axis with the smallest shifted coordinate $`s = \min(p_x, p_y, p_z)`$** ($`p = |p_0| - b`$).

1. **2-Step Tournament Minimum Partition**:
   Let $`a = \min(p_x, p_y)`$. Then with **only 4 scalar `min`/`max` comparisons and zero arithmetic**, the multiset $`\{p_x, p_y, p_z\}`$ partitions bit-exactly into:
   - its global minimum $`s = \min(a, p_z) = \min(p_x, p_y, p_z)`$, and
   - its remaining two coordinates $`u = (\max(a, p_z), \, \max(p_x, p_y)) \in \mathbb{R}^2`$.
2. **Closest-Beam Dominance Theorem**:
   For $`e > 0`$, the transverse bar distance $`q(t) = |t + e| - e`$ satisfies $`q(t) \ge t`$ everywhere, with strict inequality $`q(t) > t \iff t < -e`$, on which $`q(t) = -t - 2e`$ is **strictly decreasing** in $`t`$ while $`t`$ is **strictly increasing** in $`t`$. Because the 3D box SDF is symmetric and coordinate-wise non-decreasing, replacing $`q(p_i)`$ with $`p_i`$ achieves the maximum reduction when $`p_i`$ is minimal:
   $`\displaystyle \min_{i \in \{x,y,z\}} \mathrm{sdBox}(p_i, q(p_j), q(p_k)) \equiv \mathrm{sdBox}\big(\min(a, p_z), \; q(u_x), \; q(u_y)\big).`$
3. **Interior Longitudinal Elimination**:
   Because $`q(t) = |t + e| - e \ge -e`$ for all $`t \in \mathbb{R}`$ and $`q(t) = t`$ for $`t \ge -e`$, the minimum coordinate $`s = \min(a, p_z)`$ satisfies $`s \le \max(q(u_x), q(u_y))`$ **identically across all of $`\mathbb{R}^3`$**, eliminating $`s`$ from the interior term $`\min(\max(q_x, q_y), 0.0)`$. **Correction:** that interior term is exact only outside the 8 corner cubes; see *Interior Correction* below for the exact interior now used.
4. **Companion Cylindrical-Edge Wireframe (`sdBoxEdges`)**:
   Replacing the square strut cross-section with a circular cross-section of radius $`r`$ yields the exact 3-line, 1-`sqrt` Euclidean SDF to all 12 edges of a 3D box: $`\sqrt{\|p\|^2 - \min(\min(p_x, p_y, p_z), 0)^2} - r`$.

## Interior Correction: Corner Cubes (Claude Fable 5.2, sdf-gems 2026)

The interior term $`\min(\max(q_x, q_y), 0)`$ of the closest-beam form (and of the 3-beam baseline it reproduces) treats the inner planes $`w_j = -2e`$ as surface everywhere. Inside the 8 corner cubes ($`w_i \in [-2e, 0]`$ on all three axes) those planes lie in the interior of the solid, so the depth was under-estimated there (a conservative bound, not the exact distance) and the implied closest point $`p - d\,\nabla d`$ landed inside the solid.

A point of the hollow needs **two** coordinates below $`b - 2e`$, so with $`\mathbf{t} = \max(\mathbf{w} + 2e, 0)`$ the exact interior depth is

```math
D = \min\Big(-\max_i w_i,\; \sqrt{t_{(1)}^2 + t_{(2)}^2}\Big), \qquad d = -D,
```

where $`t_{(1)} \le t_{(2)}`$ are the two smallest components of $`\mathbf{t}`$ (an axis with $`b_i \le 2e`$ has no hollow side and is excluded). Outside the corner cubes $`t_{(1)} = 0`$ and this reduces to the previous formula; the exterior branch is unchanged.

| Check (literal GLSL, FP32) | Before | After |
| :--- | :---: | :---: |
| Max interior error vs FP64 closed form (24k points, 3 shapes) | `1.4e-2` to `2.8e-2` | `1.1e-7` |
| Implied closest point off the surface, interior | up to `9.1e-2` (1.4% to 3.1% of interior points above `1e-2`) | within the `4.5e-3` resolution of the reference cloud |
| Self-projection $`\vert d(p - d\,\nabla d)\vert`$ (90k points, 4 shapes incl. $`b_i \le 2e`$) | fails in corner cubes | `1.5e-7` |
| 1-Lipschitz ratio / sign mismatches vs inside predicate | — | `1.00000` / `0` |

The independent reference is a surface point cloud of the union of the 12 struts built from the inside predicate only (spacing `0.003`), which the result matches to its resolution (`3.5e-3`) on both sides. Cost: the interior now uses 1 `sqrt` (previously 0).

---

## GLSL Implementation

```glsl
// 1. Exact 1-Sqrt Square-Strut Box Frame (struts occupy |p_i| in [b_i - 2e, b_i] on two axes)
//    Outside the solid: closest-beam reduction, length(max(vec3(min w, q), 0)).
//    Inside the solid: the nearest non-solid point is either across an outer face (-max w) or in the hollow,
//    which needs TWO coordinates below b - 2e: sqrt of the two smallest t_i^2, t = max(w + 2e, 0).
//    (The former interior term min(max(q.x, q.y), 0) treated the inner planes w_j = -2e as surface inside the
//    8 corner cubes, where they are interior to the solid: it under-estimated the depth there by up to
//    (sqrt(2) - 1) * 2e. Fixed with Claude Fable 5.2, sdf-gems 2026.)
float sdBoxFrame( vec3 p, vec3 b, float e )
{
    p = abs(p) - b;
    float a = min(p.x, p.y);
    vec2  q = abs(vec2(max(a, p.z), max(p.x, p.y)) + e) - e;
    float l = length(max(vec3(min(a, p.z), q), 0.0));
    if( l > 0.0 ) return l;
    // An axis with b_i <= 2e has no hollow side: its struts span the whole box.
    vec3  t  = vec3((b.x > 2.0 * e) ? max(p.x + 2.0 * e, 0.0) : 1e15,
                    (b.y > 2.0 * e) ? max(p.y + 2.0 * e, 0.0) : 1e15,
                    (b.z > 2.0 * e) ? max(p.z + 2.0 * e, 0.0) : 1e15);
    float lo = min(t.x, t.y), hi = max(t.x, t.y);
    float h  = length(vec2(min(lo, t.z), max(lo, min(hi, t.z))));   // two smallest of t
    return -min(-max(p.x, max(p.y, p.z)), h);
}

// 2. Exact 3-Line 1-Sqrt Cylindrical-Edge Wireframe Cage
float sdBoxEdges( vec3 p, vec3 b, float r )
{
    p = abs(p) - b;
    float s = min(min(p.x, min(p.y, p.z)), 0.0);
    return sqrt(dot(p, p) - s * s) - r;
}
```

---

## Improvement & Optimization Analysis (Discovered using Gemini 4 Argon (sdf-gems, 2026))

![Box Frame 1-Sqrt Reduction](../images/comparisons/box_frame_comparison.png)

| Metric | Baseline `sdBoxFrame` (3 Beams) | `sdf-gems` `sdBoxFrame` (Closest-Beam Reduction) |
| :--- | :---: | :---: |
| **`sqrt` (`length`) Calls** | `3` | **`1` (`3.0x` fewer `sqrt`s)** |
| **Scalar `min` / `max` Ops** | `23` | **`9` (`2.56x` fewer)** |
| **Scalar `abs` Ops** | `6` | **`5`** |
| **Total Scalar Operations** | `56` | **`29` (`1.93x` fewer total ops)** |
| **Max Error vs. 3-Beam Baseline** | — | **`0.0` (`8.88e-16` in FP64 across $`10^6`$ points) for the closest-beam reduction; the shipped interior additionally differs from the baseline inside the corner cubes, where the baseline is only a bound** |

### Code Diff vs. Baseline (Inigo Quilez)
```diff
 float sdBoxFrame( vec3 p, vec3 b, float e )
 {
     p = abs(p) - b;
-    vec3 q = abs(p + e) - e;
-    return min(min(
-        length(max(vec3(p.x, q.y, q.z), 0.0)) + min(max(p.x, max(q.y, q.z)), 0.0),
-        length(max(vec3(q.x, p.y, q.z), 0.0)) + min(max(q.x, max(p.y, q.z)), 0.0)),
-        length(max(vec3(q.x, q.y, p.z), 0.0)) + min(max(q.x, max(q.y, p.z)), 0.0));
+    float a = min(p.x, p.y);
+    vec2  q = abs(vec2(max(a, p.z), max(p.x, p.y)) + e) - e;
+    return length(max(vec3(min(a, p.z), q), 0.0)) + min(max(q.x, q.y), 0.0);
 }
```
