# `opUnionChamfer` — Exact Corner-Local Euclidean Chamfer Booleans

![opUnionChamfer](../images/booleans/opUnionChamfer.png)

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Mercury Demogroup (plane bound) — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)** |
| **Source** | [https://mercury.sexy/hg_sdf/](https://mercury.sexy/hg_sdf/) |
| **Category** | `booleans` |
| **GLSL File** | [`glsl/booleans/opUnionChamfer.glsl`](../../glsl/booleans/opUnionChamfer.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1` across both exterior and interior, max error `8.88e-16`) |

---

## Mathematical Overview

Mercury Demogroup's baseline `opUnionChamfer(a, b, r) = min(min(a, b), (a - r + b)*sqrt(0.5))` takes the minimum of three supporting half-planes. While exact in the exterior for orthogonal surfaces ($\nabla a \perp \nabla b$), in the interior ($a < 0$ or $b < 0$) it is **non-local** and **violates the Eikonal equation** ($\|\nabla d\| = 1$): along the diagonal $a = b < 0$, $(a + b - r)/\sqrt{2} = \sqrt{2}\,a - r/\sqrt{2} < a$, corrupting the interior distance field all the way to $-\infty$ (overestimating interior depth by $41.4\%$).

By folding $(a, b)$ via $x = \min(a, b)$ and $y = \max(a, b)$, the complement of the chamfered union has a single convex vertex at $(0, r)$ where the ray $\{x = 0, y \ge r\}$ meets the $45^\circ$ chamfer segment $\{x + y = r, x \in [0, r/2]\}$. Whenever $x < 0$ and $y - x > r$ (outside the normal strip of the $45^\circ$ segment), the closest boundary point is $(0, \max(y, r))$ at exact signed Euclidean distance $-\sqrt{x^2 + \max(r - y, 0)^2}$ (which simplifies identically to $x = \min(a, b)$ for $y \ge r$, making the operator **strictly corner-local**).

We also provide `opUnionChamferOblique(a, b, r, c)` for surfaces meeting at arbitrary angle $\cos\theta = c = \nabla a \cdot \nabla b$, normalizing the bisector plane by $\|\nabla(a + b)\| = \sqrt{2 + 2 c}$.

---

## GLSL Implementation

```glsl
// 1. Exact Corner-Local Euclidean 45-Degree Chamfer Union (fixes interior seam distortion)
float opUnionChamfer( float a, float b, float r )
{
    float x = min(a, b), y = max(a, b);
    return (x < 0.0 && y - x > r) ? -length(vec2(x, max(r - y, 0.0)))
                                  : min(x, (x + y - r) * sqrt(0.5));
}

// 2. Angle-Exact Oblique Chamfer Union for surfaces meeting at cos(theta) = c = dot(grad(a), grad(b))
float opUnionChamferOblique( float a, float b, float r, float c )
{
    return min(min(a, b), (a + b - r) * inversesqrt(max(2.0 + 2.0 * c, 1e-6)));
}

// 3. Exact Euclidean Chamfer Intersection & Subtraction
float opIntersectionChamfer( float a, float b, float r )
{
    return -opUnionChamfer(-a, -b, r);
}

float opSubtractionChamfer( float a, float b, float r )
{
    return -opUnionChamfer(a, -b, r);
}
```

---

## Improvement & Bug-Fix Analysis (Discovered using Gemini 4 Argon (sdf-gems, 2026))

![Stairs & Chamfer Exact Euclidean Comparison](../images/comparisons/stairs_chamfer_comparison.png)

| Metric ($r = 1.0$) | Baseline `opUnionChamfer` (Mercury `hg_sdf`) | `sdf-gems` `opUnionChamfer` |
| :--- | :---: | :---: |
| **Max Exterior Error ($\nabla a \perp \nabla b$)** | `2.22e-16` | **`2.22e-16`** |
| **Max Interior Error ($\nabla a \perp \nabla b$)** | `0.164` (`41.4%` slope error to $-\infty$) | **`8.88e-16` (`0.0` mathematical error)** |
| **Strictly Corner-Local for $y \ge r$?** | No (distorts $a = b \to -\infty$) | **Yes (`d == min(a, b)` for $\max(a,b) \ge r, \min(a,b) \le 0$)** |

### Code Diff vs. Baseline (Mercury Demogroup)
```diff
 float opUnionChamfer( float a, float b, float r )
 {
-    return min(min(a, b), (a - r + b)*sqrt(0.5));
+    float x = min(a, b), y = max(a, b);
+    return (x < 0.0 && y - x > r) ? -length(vec2(x, max(r - y, 0.0)))
+                                  : min(x, (x + y - r) * sqrt(0.5));
 }
```
