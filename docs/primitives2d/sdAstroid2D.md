# `sdAstroid2D`, `sdgAstroid2D`, `sdAstroidRevolved3D`, `sdgAstroidRevolved3D` & `sdAstroid3D` — Exact 0-`pow`, 0-Trig 2D & 3D Astroid

![sdAstroid2D Preview](../images/primitives2d/sdAstroid2D.png)

## Overview

The **Astroid** $|x|^{2/3} + |y|^{2/3} = r^{2/3}$ (the 4-cusped hypocycloid, $\ell_{2/3}$ unit ball, and evolute of the ellipse) and its 3D counterparts — the **Revolved Astroid** $(x^2 + z^2)^{1/3} + |y|^{2/3} = r^{2/3}$ and the **Octahedral 3D Astroid** $|x|^{2/3} + |y|^{2/3} + |z|^{2/3} = r^{2/3}$ — are classic non-convex algebraic curves and surfaces with sharp cusps.

Standard shader implementations evaluate the algebraic level set $(|x|^{2/3} + |y|^{2/3})^{3/2} - r$, which is **not** a Euclidean SDF: it underestimates the distance at the origin by $50\%$ (`0.5 * r`) and its gradient norm blows up to $\infty$ along the coordinate axes $x = 0$ and $y = 0$. We derive the first **exact, 0-`pow`, 0-trig Euclidean SDF and analytical gradient** (`sdAstroid2D`, `sdgAstroid2D`, `sdAstroidRevolved3D`, `sdgAstroidRevolved3D`, `sdAstroid3D`).

---

## Mathematical Formulation

![Astroid Comparison](../images/comparisons/astroid_comparison.png)

### 1. Theorem: Pythagorean Co-Normal Projection Identity

In the folded first octant $p_x \ge p_y \ge 0$, the closest point on the astroid lies in $\theta \in [0, \pi/4]$, parameterized algebraically by $s = \sin\theta \in [0, 1/\sqrt{2}]$, $c = \sqrt{1 - s^2}$, with $\mathbf{q}(s) = (r c^3, r s^3)$. Because $\mathbf{q}'(\theta) = 3 r c s (-c, s)$, the unit tangent and unit outward normal are simply:
$$\hat{\mathbf{t}}(s) = (-c, s), \qquad \hat{\mathbf{n}}(s) = (s, c).$$
Projecting $\mathbf{p} - \mathbf{q}(s)$ onto $\{\hat{\mathbf{n}}(s), \hat{\mathbf{t}}(s)\}$ and using $c^2 + s^2 = 1$ and $c^4 - s^4 = 1 - 2s^2$ yields:
$$d_n(s) = \langle \mathbf{p} - \mathbf{q}(s), \hat{\mathbf{n}}(s) \rangle = p_x s + p_y c - r s c, \qquad f(s) = \langle \mathbf{p} - \mathbf{q}(s), \hat{\mathbf{t}}(s) \rangle = p_y s - p_x c + r(1 - 2s^2),$$
with Parseval's identity $\|\mathbf{p} - \mathbf{q}(s)\|_2^2 \equiv d_n(s)^2 + f(s)^2$. At the root $f(s^*) = 0$, the exact signed Euclidean distance is $d_n(s^*)$ and the exact outward unit gradient is $(s^*, c^*)$ with **0 extra `sqrt`s**.

### 2. Theorem: Secant-Envelope Cusp Early Exit & Monotone Root Bracket

By concavity of $u \mapsto \sqrt{1 - u}$ on $u = s^2 \in [0, 1/2]$, $\sqrt{1 - s^2} \ge 1 - (2 - \sqrt{2})s^2$ with equality at both endpoints $s = 0$ and $s = 1/\sqrt{2}$. Hence $f(s)$ is globally bounded above by the **Secant-Envelope Quadratic**:
$$f(s) \le Q(s) = -A s^2 + p_y s + B, \qquad A = 2r - (2 - \sqrt{2})p_x, \qquad B = r - p_x.$$
- **0-Step Cusp Early Exit ($38.4\%$ of `[-2, 2]^2`)**: Whenever $\Delta = p_y^2 + 4 A B \le 0$ or $A \le 0$, $f(s) \le Q(s) \le 0$ on all of $[0, 1/\sqrt{2}]$, proving that the cusp $(r, 0)$ is the unique global minimum (`0` Halley steps, `1` `inversesqrt`).
- **Monotone Root Bracket**: When $\Delta > 0$ and $A > 0$, the upper root $s_0 = \min\big(\frac{p_y + \sqrt{\Delta}}{2A}, \frac{1}{\sqrt{2}}\big)$ satisfies $f(s_0) \le Q(s_0) = 0$, landing strictly on the monotone side $s_0 \ge s^*$ where $f'(s) < 0$.

---

## Numerical Verification

| Function | `pow` / Trig Calls | Max Distance Error | Mean Distance Error | Max Smooth Gradient Error |
| :--- | :---: | :---: | :---: | :---: |
| **Algebraic $(|x|^{2/3}+|y|^{2/3})^{3/2}-r$** | `2` `pow` | `5.00e-1` ($50\%$ at origin) | `1.42e-1` | $\infty$ (blows up on axes) |
| **Ours `sdAstroid2D` / `sdgAstroid2D`** | **`0` `pow`, `0` trig** | **`8.88e-16`** | **`2.91e-17`** | **`2.07e-8`** |
| **Ours `sdAstroidRevolved3D` / `sdg...`** | **`0` `pow`, `0` trig** | **`8.88e-16`** | **`3.14e-17`** | **`2.51e-9`** |
| **Ours `sdAstroid3D` (Octahedral)** | **`0` `pow`, `0` trig** | **`4.31e-10`** | **`1.35e-12`** | — |

---

## GLSL Implementation

See [`glsl/primitives2d/sdAstroid2D.glsl`](../../glsl/primitives2d/sdAstroid2D.glsl).
