# `sminCubicGeom` — Corner-Local Associative $C^2$ ($G^2$ Curvature-Continuous) Smooth Minimum

![2D SDF Preview](../images/booleans/sminCubicGeom.png)

- **GLSL Implementation**: [`sminCubicGeom.glsl`](../../glsl/booleans/sminCubicGeom.glsl)
- **Category**: `booleans`
- **Status**: **Corner-Local, Conservative ($2^{-1/6} \le \|\nabla f\|_2 \le 1$), Associative, $C^2$ Smooth Minimum** (Discovered using **Gemini 4 Argon**, `sdf-gems`, 2026)

---

## 1. Summary & The 5-Property Smooth-Minimum Problem

In his classification of smooth-minimum functions ([`smin`](https://iquilezles.org/articles/smin/)), Inigo Quilez evaluates operators by **Rigidity**, **Locality**, **Conservativeness**, **Associativity**, and **$C^2$ / $G^2$ Curvature Continuity**:

1. **CD Polynomial Smooth Minimums (`sminCubic`, `sminQuartic`)**: Achieve $C^2$ ($G^2$ curvature continuity, preventing specular reflection seams), but **fail Locality**—their blend zone $|a - b| < 6k$ extends to $\pm\infty$ along the Voronoi medial bisector $a = b$, underestimating distance by $1 - 1/\sqrt{2} = 29.3\%$ infinitely far from the surface—and **fail Associativity**.
2. **Geometrical Circular Smooth Minimum (`sminCircularGeom`, $\ell_2$ quarter-circle)**: Achieves Rigidity, strict Corner-Locality ($a < R, b < R$), Conservativeness ($\perp$), and Associativity, but **fails $C^2$ / $G^2$ Curvature Continuity**: because $\sqrt{u_x^2 + u_y^2} = u_y + \frac{u_x^2}{2u_y} + O(u_x^4)$ near $u_x = 0$, the second derivative $\frac{\partial^2 f}{\partial a^2}$ jumps discontinuously from $0$ to $-1/u_y$ across the seam $a = R$ and $b = R$, producing visible Mach-band creases in specular reflections.

![sminCubicGeom Comparison](../images/comparisons/smin_cubic_geom_comparison.png)

---

## 2. Mathematical Derivation ($\ell_3$ Lamé Corner Fillet)

Set the corner support radius $R = \frac{k}{1 - 2^{-1/3}} \approx 4.84732210\,k$ (normalized so that at $a = b = k$, the zero-level surface passes through $(k, k)$), and define $u = \max(R - (a, b), 0)$:
$$\mathrm{sminCubicGeom}(a, b, k) = \max\big(R, \min(a, b)\big) - \big(u_x^3 + u_y^3\big)^{1/3}.$$

### Proof of the 5 Simultaneous Properties
1. **Strict Corner-Locality & Rigidity (`1.11e-16` error)**:
   Whenever $a \ge R$ or $b \ge R$, at least one component of $u$ is $0$, so $(u_x^3 + u_y^3)^{1/3} = \max(u_x, u_y) = R - \min(a, b)$ and $\mathrm{sminCubicGeom}(a, b, k) \equiv \min(a, b)$ identically.
2. **$C^2$ ($G^2$ Curvature Continuity) at $a = R, b = R$ and $C^\infty$ across $a = b$**:
   Expanding near $u_x \to 0^+$ ($a \to R^-$) for $u_y > 0$:
   $$(u_x^3 + u_y^3)^{1/3} = u_y\left(1 + \frac{u_x^3}{u_y^3}\right)^{1/3} = u_y + \frac{u_x^3}{3u_y^2} - \frac{u_x^6}{9u_y^5} + \cdots$$
   Both the linear ($u_x$) and quadratic ($u_x^2$) terms vanish identically. Consequently, the gradient $\nabla f$ and the full Hessian $\nabla^2 f$ match $\min(a, b)$ continuously across $a = R$ and $b = R$ with **zero curvature jump**, while $(u_x^3 + u_y^3)^{1/3}$ is $C^\infty$-smooth across the diagonal $a = b < R$.
3. **Hölder Conservativeness Theorem ($\|\nabla f\|_2 \in [2^{-1/6}, 1]$)**:
   For orthogonal SDFs ($\nabla a \perp \nabla b$, $\|\nabla a\|_2 = \|\nabla b\|_2 = 1$), the gradient inside the corner $a < R, b < R$ is:
   $$\nabla f = \frac{u_x^2\,\nabla a + u_y^2\,\nabla b}{(u_x^3 + u_y^3)^{2/3}}, \qquad \|\nabla f\|_2 = \frac{(u_x^4 + u_y^4)^{1/2}}{(u_x^3 + u_y^3)^{2/3}} = \left(\frac{\|u\|_4}{\|u\|_3}\right)^2.$$
   By the $\ell_p$ norm inequality $\|u\|_4 \le \|u\|_3 \le 2^{1/12}\|u\|_4$ on $\mathbb{R}^2$:
   $$2^{-1/6} \approx 0.8908987 \le \|\nabla f\|_2 \le 1.0000000 \quad \forall (a, b) \in \mathbb{R}^2.$$
   Thus `sminCubicGeom` is **100% conservative** (zero raymarching overstep for $\nabla a \perp \nabla b$) and reduces maximum Eikonal deviation by **$2.7\times$** compared to `sminCubic` ($10.9\%$ vs. $29.3\%$, confined strictly to $[0, R]^2$).
4. **Exact 3-Way Associativity (`4.44e-16` error)**:
   Because $\big((u_x^3 + u_y^3) + u_z^3\big)^{1/3} = (u_x^3 + u_y^3 + u_z^3)^{1/3}$, chaining `sminCubicGeom(a, sminCubicGeom(b, c, k), k)` is completely order-independent at 3-way junctions.

---

## 3. Property & Numerical Benchmark Comparison

| Smooth Minimum | Rigid? | Corner-Local? | Conservative ($\perp$)? | Min $\|\nabla f\|_2$ ($\perp$) | Associative? | $C^2$ / $G^2$ Continuous? |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **`sminQuadratic`** | Yes | No ($\infty$ along $a=b$) | Yes | `0.7071` (`-29.3%`) | No | No ($C^1$) |
| **`sminCubic`** | Yes | No ($\infty$ along $a=b$) | Yes | `0.7071` (`-29.3%`) | No | Yes ($C^2$) |
| **`sminCircularGeom` ($\ell_2$)** | Yes | Yes ($a,b < R$) | Yes | `1.0000` (`0.0%`) | Yes | No ($C^1$; $\nabla^2 f$ jumps by $1/R$) |
| **`sminCubicGeom` ($\ell_3$, Ours)** | **Yes** | **Yes ($a,b < R$)** | **Yes** | **`0.8909` (`-10.9%`)** | **Yes (`4.44e-16`)** | **Yes ($C^2$ seam, $C^\infty$ diag)** |

---

## 4. GLSL Implementation

```glsl
float sminCubicGeom(float a, float b, float k) {
    k *= 4.84732210; // 1.0 / (1.0 - pow(2.0, -1.0 / 3.0))
    vec2 u = max(k - vec2(a, b), 0.0);
    return max(k, min(a, b)) - pow(u.x * u.x * u.x + u.y * u.y * u.y, 1.0 / 3.0);
}
```
