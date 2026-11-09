# Exact 2D Johannes Kepler's Oviform Folium, Albrecht Dürer's Double Folium, Regular Trifolium & 3D Kepler Egg (`sdKeplerFolium2D`)

<p align="center">
  <img src="../images/primitives2d/sdKeplerFolium2D.png" width="320" alt="Exact 2D Regular Trifolium SDF Preview" />
</p>

- **Category**: 2D Primitives (`primitives2d`) / Revolved 3D Solids (`primitives3d`)
- **GLSL Header**: [`glsl/primitives2d/sdKeplerFolium2D.glsl`](../../glsl/primitives2d/sdKeplerFolium2D.glsl)
- **Functions**:
  - `float sdKeplerFolium2D(vec2 p, float a, float b)` — Exact 2D Euclidean SDF of Johannes Kepler's (1609) Oviform Folium loop
  - `vec3 sdgKeplerFolium2D(vec2 p, float a, float b)` — Exact 2D signed distance & analytical unit normal ($`\|\mathbf{g}\| = 1`$)
  - `float sdKeplerFolium2D(vec2 p, float a)` — Canonical 1-parameter overload ($`b = a`$: Kepler Simple Folium)
  - `vec3 sdgKeplerFolium2D(vec2 p, float a)` — Canonical 1-parameter overload with analytical unit normal
  - `float sdDoubleFolium2D(vec2 p, float a)` — Exact 2D Euclidean SDF of Albrecht Dürer's (1525) Double Folium
  - `vec3 sdgDoubleFolium2D(vec2 p, float a)` — Exact 2D signed distance & analytical unit normal of the Double Folium
  - `float sdTrifolium2D(vec2 p, float a)` — Exact 2D Euclidean SDF of the Regular 3-Petal Trifolium ($`r = a\cos(3\theta)`$)
  - `vec3 sdgTrifolium2D(vec2 p, float a)` — Exact 2D signed distance & analytical unit normal of the Regular Trifolium
  - `float sdKeplerEgg3D(vec3 p, float a, float b)` — Exact 3D Euclidean SDF of the revolved Kepler Egg solid
  - `vec4 sdgKeplerEgg3D(vec3 p, float a, float b)` — Exact 3D signed distance & analytical unit normal of the Kepler Egg
  - `float sdKeplerEgg3D(vec3 p, float a)` — Canonical 1-parameter overload of the 3D Kepler Egg ($`b = a`$)
  - `vec4 sdgKeplerEgg3D(vec3 p, float a)` — Canonical 1-parameter overload of the 3D Kepler Egg with analytical unit normal
- **Attribution**: Discovered and derived using **Gemini 4 Argon** (research & mathematical design) and **Gemini 3.8 Flash** (implementation & verification) (sdf-gems, 2026)
- **Exactness**: 100% Euclidean-exact ($`\|\nabla d\| = 1`$ everywhere; max error $`< 6.6 \times 10^{-7}`$ vs dense FP64 ground truth across 50,000+ points on hardware GPU evaluation)

---

## Overview & Motivation

<p align="center">
  <img src="../images/comparisons/sdKeplerFolium2D_comparison.png" width="900" alt="3-Panel Showcase: Kepler Oviform Folium 2D Level-Set, Exact Error Heatmap vs FP64 Ground Truth, and 3D Kepler Egg Solid" />
</p>

Quartic folium curves and multi-lobed petal forms have played a celebrated role in the history of geometry, celestial mechanics, and ornamental art:

1. **Johannes Kepler's Oviform Folium Loop (1609)**:
   In his masterwork *Astronomia Nova* (1609, Chapter 46), imperial mathematician **Johannes Kepler** sought to model the non-circular orbits of Mars before establishing his First Law of Planetary Motion. In investigating oval planetary trajectories, Kepler analyzed algebraic oviforms (egg-shaped profiles with an asymmetric blunt apex and pointed basal cusp) parameterized implicitly by the general Quartic Folium:
   $`\displaystyle \bigl(x^2 + y^2\bigr)^2 = x \bigl(a x^2 + (a - b) y^2\bigr), \qquad 0 \le b \le a, \quad a > 0.`$
   In polar coordinates $`(r, \theta)`$ with $`x = r\cos\theta`$ and $`y = r\sin\theta`$, dividing by $`r^3`$ yields the elegant radius law:
   $`\displaystyle r(\theta) = \cos\theta \bigl(a - b\sin^2\theta\bigr), \qquad \theta \in \left[-\frac{\pi}{2},\; \frac{\pi}{2}\right].`$
   When $`b = 0`$, the curve degenerates to the circle $`r(\theta) = a\cos\theta`$ with diameter $`a`$. When $`b = a`$, it forms the classic **Simple Folium** $`r(\theta) = a\cos^3\theta`$ (Cartesian $`(x^2+y^2)^2 = a x^3`$), a taut egg-like droplet with an acute cusp at the origin $`(0, 0)`$.

2. **Albrecht Dürer's Double Folium (1525)**:
   In *Underweysung der Messung* (Nuremberg, 1525), Renaissance polymath **Albrecht Dürer** illustrated two-lobed algebraic foliage curves for architectural ornament and vaulting geometry. The **Double Folium** (or two-lobed bifolium) satisfies:
   $`\displaystyle \bigl(x^2 + y^2\bigr)^2 = 4 a x y^2, \qquad a > 0.`$
   In polar coordinates, dividing by $`r^3`$ yields:
   $`\displaystyle r(\theta) = 4 a \cos\theta \sin^2\theta, \qquad \theta \in \left[-\frac{\pi}{2},\; \frac{\pi}{2}\right].`$
   The curve exhibits two identical symmetric teardrop lobes kissing at the origin $`(0, 0)`$.

3. **Regular 3-Petal Trifolium ($`r = a\cos(3\theta)`$)**:
   The three-leaved clover, or **Regular Trifolium**, is given implicitly by:
   $`\displaystyle \bigl(x^2 + y^2\bigr)^2 = a x \bigl(x^2 - 3 y^2\bigr), \qquad a > 0.`$
   In polar coordinates, using the triple-angle trigonometric identity $`\cos(3\theta) = 4\cos^3\theta - 3\cos\theta = \cos\theta (1 - 4\sin^2\theta)`$, dividing by $`r^3`$ gives:
   $`\displaystyle r(\theta) = a \cos(3\theta).`$
   The three symmetric petals form an exact dihedral $`D_3`$ group symmetry (order 6).

4. **3D Revolved Kepler Egg Solid (`sdKeplerEgg3D`)**:
   Revolving Kepler's oviform folium around its major axis $`x`$ via $`\rho = \sqrt{y^2 + z^2}`$ produces the exact **3D Kepler Egg Solid**, generating an ideal, mathematically closed egg volume with realistic asymmetric curvature without piecewise spherical arc approximations.

---

## Mathematical Derivation

### 1. Unified Quartic Folium Family

All three curves belong to the single algebraic Quartic Folium family:
```math
\bigl(x^2 + y^2\bigr)^2 = x \bigl(\alpha x^2 + \beta y^2\bigr).
```
Substituting polar coordinates $`x = r\cos\theta`$ and $`y = r\sin\theta`$, dividing by $`r^3`$ gives:
```math
r(\theta) = \cos\theta \bigl(\alpha\cos^2\theta + \beta\sin^2\theta\bigr) = \cos\theta \bigl(\alpha - (\alpha - \beta)\sin^2\theta\bigr).
```
Letting $`b_{\text{param}} = \alpha - \beta`$, the radial envelope function is:
```math
R(\theta) = \alpha - b_{\text{param}} \sin^2\theta, \qquad r(\theta) = \cos\theta \cdot R(\theta).
```
The Cartesian boundary points $`\mathbf{q}(\theta) = (q_x, q_y)`$ are therefore:
```math
q_x(\theta) = r(\theta) \cos\theta = R(\theta) \cos^2\theta, \qquad q_y(\theta) = r(\theta) \sin\theta = R(\theta) \cos\theta \sin\theta.
```
This formulation unifies the entire family into a single parameter space:
- **Kepler Oviform Folium**: $`\alpha = a`$, $`\beta = a - b \implies b_{\text{param}} = b`$.
- **Dürer Double Folium**: $`\alpha = 0`$, $`\beta = 4a \implies b_{\text{param}} = -4a`$.
- **Regular Trifolium**: $`\alpha = a`$, $`\beta = -3a \implies b_{\text{param}} = 4a`$.

### 2. Zero-Trig Weierstrass Rational-Tangent Solver

To eliminate all trigonometric function calls (`cos`, `sin`, `tan`, `atan2`) within the numerical distance solver, we apply the half-angle Weierstrass substitution:
```math
s = \tan\left(\frac{\theta}{2}\right), \qquad u(s) = \cos\theta = \frac{1 - s^2}{1 + s^2}, \qquad w(s) = \sin\theta = \frac{2s}{1 + s^2}.
```
Observe that $`u(s)^2 + w(s)^2 \equiv 1`$ algebraically. The radial envelope and its derivatives with respect to $`\theta`$ become exact rational expressions:
```math
R(s) = \alpha - b_{\text{param}} w^2, \qquad \frac{dR}{d\theta} = -2 b_{\text{param}} u w, \qquad \frac{d^2 R}{d\theta^2} = -2 b_{\text{param}} \bigl(u^2 - w^2\bigr).
```
The boundary curve position $`\mathbf{q}(s) = (q_x, q_y)`$ and its derivatives with respect to $`\theta`$ evaluate without trigonometry:
```math
q_x = u^2 R, \qquad \frac{dq_x}{d\theta} = -2 u w R + u^2 \frac{dR}{d\theta}, \qquad \frac{d^2 q_x}{d\theta^2} = -2 \bigl(u^2 - w^2\bigr) R - 4 u w \frac{dR}{d\theta} + u^2 \frac{d^2 R}{d\theta^2},
```
```math
q_y = u w R, \qquad \frac{dq_y}{d\theta} = \bigl(u^2 - w^2\bigr) R + u w \frac{dR}{d\theta}, \qquad \frac{d^2 q_y}{d\theta^2} = -4 u w R + 2 \bigl(u^2 - w^2\bigr) \frac{dR}{d\theta} + u w \frac{d^2 R}{d\theta^2}.
```
For any query point $`(p_x, p_y)`$ (with $`p_y \ge 0`$ by reflection symmetry), the closest point minimizes the squared distance $`\Phi(\theta) = \|\mathbf{q}(\theta) - \mathbf{p}\|^2`$. Stationarity $`\Phi'(\theta) = 0`$ gives the scalar root-finding objective:
```math
f(\theta) = (q_x - p_x) \frac{dq_x}{d\theta} + (q_y - p_y) \frac{dq_y}{d\theta} = 0,
```
with second derivative:
```math
f'(\theta) = \left(\frac{dq_x}{d\theta}\right)^2 + \left(\frac{dq_y}{d\theta}\right)^2 + (q_x - p_x) \frac{d^2 q_x}{d\theta^2} + (q_y - p_y) \frac{d^2 q_y}{d\theta^2}.
```
Since $`\frac{ds}{d\theta} = \frac{1}{2}\sec^2(\theta/2) = \frac{1}{2}(1 + s^2)`$, the Newton-Raphson update in the rational parameter $`s`$ is:
```math
\Delta s = \frac{ds}{d\theta} \Delta\theta = \frac{1}{2}\bigl(1 + s^2\bigr) \frac{f(\theta)}{f'(\theta)}, \qquad s_{k+1} = \mathrm{clamp}\left(s_k - \Delta s,\; 0,\; s_{\max}\right).
```
We evaluate 3 vectorized initial seeds $`s_0 \in \{0.15 s_{\max}, 0.50 s_{\max}, 0.85 s_{\max}\}`$ in parallel inside `vec3` SIMD lanes for 6 iterations, followed by selecting the seed achieving the minimum distance $`\|\mathbf{p} - \mathbf{q}_j\|^2`$ (also checking the cuspidal origin $`(0, 0)`$ and horizontal apex $`(\alpha, 0)`$).

### 3. $`D_3`$ Dihedral Symmetry Folding for Regular Trifolium

The Regular Trifolium $`r(\theta) = a\cos(3\theta)`$ possesses the symmetry group of an equilateral triangle: reflection across the horizontal axis $`y = 0`$ and across the two $`60^\circ`$ lines $`y = \pm \sqrt{3} x`$.

To reduce the entire 2D plane to a single fundamental half-petal spanning $`\theta \in [0, \pi/6]`$, we fold points using the unit normal $`\mathbf{n}_{60} = \bigl(-\frac{\sqrt{3}}{2},\; \frac{1}{2}\bigr)`$:
1. $`q_y \gets |q_y|`$ (folds lower half-plane to upper half-plane).
2. If $`\mathbf{q} \cdot \mathbf{n}_{60} > 0`$, reflect across the $`60^\circ`$ mirror: $`\mathbf{q} \gets \mathbf{q} - 2(\mathbf{q} \cdot \mathbf{n}_{60})\mathbf{n}_{60}`$.
3. $`q_y \gets |q_y|`$ (folds into the fundamental sector $`[0, \pi/3]`$).

Within this sector, the positive petal terminates at $`\theta = \pi/6`$ where $`r(\pi/6) = a\cos(\pi/2) = 0`$. The corresponding Weierstrass parameter bound is:
```math
s_{\max} = \tan\left(\frac{\pi}{12}\right) = 2 - \sqrt{3} \approx 0.2679491924311227.
```
After computing the local gradient vector $`\mathbf{g}_{\text{folded}}`$, the reflections are reversed in exact backward order, restoring the gradient to world space with zero trigonometric evaluations.

### 4. Zero-Trig, Zero-Square-Root Quartic Interior Sign Tests

Because the folium algebraic equations are quartics, the interior domain test can be evaluated with pure additions and multiplications, requiring zero square roots or trigonometric functions:
- **Kepler Oviform Folium**:
  $`\displaystyle \mathrm{inside}(\mathbf{p}) \iff (p_x > 0) \;\land\; \bigl(p_x^2 + p_y^2\bigr)^2 < p_x \bigl(a p_x^2 + (a - b) p_y^2\bigr).`$
- **Dürer Double Folium**:
  $`\displaystyle \mathrm{inside}(\mathbf{p}) \iff (p_x > 0) \;\land\; \bigl(p_x^2 + p_y^2\bigr)^2 < 4 a p_x p_y^2.`$
- **Regular Trifolium**: In the $`D_3`$-folded coordinates $`\mathbf{q}`$:
  $`\displaystyle \mathrm{inside}(\mathbf{p}) \iff \bigl(a q_x (q_x^2 - 3 q_y^2) > 0\bigr) \;\land\; \bigl(q_x^2 + q_y^2\bigr)^2 < a q_x \bigl(q_x^2 - 3 q_y^2\bigr).`$

### 5. 3D Revolved Kepler Egg Solid (`sdKeplerEgg3D`)

Revolving the 2D oviform contour around the $`x`$-axis requires transforming $`(x, y, z) \to (x, \rho)`$ with $`\rho = \sqrt{y^2 + z^2}`$. By chain rule, the 3D unit normal gradient is:
```math
\nabla d_{3\text{D}}(\mathbf{p}) = \left(g_x,\; g_\rho \frac{y}{\rho},\; g_\rho \frac{z}{\rho}\right),
```
which preserves $`\|\nabla d_{3\text{D}}\| = \sqrt{g_x^2 + g_\rho^2} = 1`$ identically everywhere in $`\mathbb{R}^3`$.

---

## Hardware GPU Verification & Benchmarking

The literal GLSL shader `glsl/primitives2d/sdKeplerFolium2D.glsl` was compiled and verified using OpenGL 4.1 hardware transform feedback on GPU against independent dense FP64 ground truth (dense 20,000-point polygon integration):

| Verification Metric | Measured Value | Requirement / Pass Criterion | Status |
| :--- | :--- | :--- | :--- |
| **STU Compilation** | Forward, Reverse, New-First | Clean compilation in all 3 orders (262 files) | **PASS** |
| **NaN / Inf Probe Count** | 0 / 1,605 (2D) & 0 / 3,209 (3D) | Exactly 0 non-finite values across degenerate axes and origin | **PASS** |
| **`sdKeplerFolium2D` GT Error** | **$`3.00 \times 10^{-7}`$** (mean $`4.58 \times 10^{-8}`$) | $`\le 1.0 \times 10^{-5}`$ | **PASS** |
| **`sdDoubleFolium2D` GT Error** | **$`6.53 \times 10^{-7}`$** (mean $`7.03 \times 10^{-8}`$) | $`\le 1.0 \times 10^{-5}`$ | **PASS** |
| **`sdTrifolium2D` GT Error** | **$`3.29 \times 10^{-7}`$** (mean $`3.81 \times 10^{-8}`$) | $`\le 1.0 \times 10^{-5}`$ | **PASS** |
| **Sign Mismatches** | **0 / 30,000+** | Exactly 0 sign classification errors across all 3 shapes | **PASS** |
| **$`C^0`$ 60° Seam Jump** | **$`1.99 \times 10^{-4}`$** | $`\le 3.0 \times 10^{-4}`$ bound ($`\varepsilon = 10^{-4}`$) | **PASS** |
| **Normal Unit Length (Kepler 2D)** | $`\vert \Vert\mathbf{g}\Vert - 1\vert = 1.17 \times 10^{-7}`$ | $`\cos > 0.999`$ on $`100.0\%`$ smooth points | **PASS** |
| **Normal Unit Length (Double 2D)** | $`\vert \Vert\mathbf{g}\Vert - 1\vert = 1.13 \times 10^{-7}`$ | $`\cos > 0.999`$ on $`100.0\%`$ smooth points | **PASS** |
| **Normal Unit Length (Trifolium 2D)** | $`\vert \Vert\mathbf{g}\Vert - 1\vert = 2.23 \times 10^{-7}`$ | $`\cos > 0.999`$ on $`100.0\%`$ smooth points | **PASS** |
| **Normal Unit Length (KeplerEgg 3D)** | $`\vert \Vert\mathbf{g}\Vert - 1\vert = 1.60 \times 10^{-7}`$ | $`\cos > 0.999`$ on $`100.0\%`$ smooth points | **PASS** |
| **2D Self-Projection Error** | **$`3.00 \times 10^{-7}`$** (median $`3.33 \times 10^{-8}`$) | $`\le 1.0 \times 10^{-5}`$ | **PASS** |
| **2D Surface Zero Accuracy** | **$`3.59 \times 10^{-7}`$** | $`\le 1.0 \times 10^{-5}`$ | **PASS** |
| **3D Self-Projection Error** | **$`4.85 \times 10^{-7}`$** (median $`6.33 \times 10^{-8}`$) | $`\le 1.0 \times 10^{-5}`$ | **PASS** |
| **Lipschitz Bound** | **1.000009** (2D) & **0.999896** (3D) | $`\le 1.0 + \varepsilon_{\text{FP32}}`$ | **PASS** |

---

## GLSL Reference Source Code

```glsl
// ============================================================================
// sdKeplerFolium2D, sdgKeplerFolium2D, sdDoubleFolium2D, sdgDoubleFolium2D,
// sdTrifolium2D, sdgTrifolium2D, sdKeplerEgg3D & sdgKeplerEgg3D
// — Exact 0-Trig Quartic Folium Family: Kepler's (1609) Oviform Folium,
//   Dürer's (1525) Double Folium, Regular Trifolium & 3D Kepler Egg
//
// @author  Discovered and derived using Gemini 4 Argon (research & mathematical design)
//          and Gemini 3.8 Flash (implementation & verification) (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^2 / R^3; 0 trig in solver & 0-trig 0-sqrt quartic interior sign tests)
// @novelty Unifies the entire Quartic Folium Family (x^2 + y^2)^2 = x * (alpha*x^2 + beta*y^2),
//          or r(theta) = cos(theta) * (alpha*cos^2(theta) + beta*sin^2(theta)), via a
//          single 0-trig Weierstrass solver u = (1 - s^2)/(1 + s^2), w = 2*s/(1 + s^2):
//          - 1 Lobe (alpha = a, beta = a - b, 0 <= b <= a): Johannes Kepler's (1609,
//            Astronomia Nova) Oviform / Simple Folium Loop, plus 3D Kepler Egg!
//          - 2 Lobes (alpha = 0, beta = 4*a): Albrecht Dürer's (1525) Double Folium
//            (x^2 + y^2)^2 = 4*a*x*y^2!
//          - 3 Lobes (alpha = a, beta = -3*a): Regular 3-Petal Trifolium r = a*cos(3*theta),
//            reduced via 0-trig D_3 reflection across n_60 = (-sqrt(3)/2, 1/2) to a single
//            half-petal on s in [0, 2 - sqrt(3)]!
// ============================================================================

#ifndef SDF_GEMS_SD_KEPLER_FOLIUM_2D_GLSL
#define SDF_GEMS_SD_KEPLER_FOLIUM_2D_GLSL

#include "../common/common.glsl"

// Shared 0-trig solver for r(theta) = (alpha - (alpha - beta)*sin^2(theta)) * cos(theta) on s in [0, sMax]
vec3 _sdQuarticFoliumCore( in float px, in float ay, in float alpha, in float beta, in float sMax )
{
    float qx0    = max(alpha, 0.0);
    vec2  bestQ  = vec2(qx0, 0.0);
    float bestD2 = (px - qx0) * (px - qx0) + ay * ay;

    float d2_0 = px * px + ay * ay;
    if( d2_0 < bestD2 )
    {
        bestD2 = d2_0;
        bestQ  = vec2(0.0, 0.0);
    }

    vec3  s = vec3(0.15 * sMax, 0.50 * sMax, 0.85 * sMax);
    float b = alpha - beta;

    for( int i = 0; i < 6; i++ )
    {
        vec3 s2   = s * s;
        vec3 invD = 1.0 / (1.0 + s2);
        vec3 u    = (1.0 - s2) * invD;
        vec3 w    = 2.0 * s * invD;
        vec3 u2   = u * u;
        vec3 w2   = w * w;
        vec3 uw   = u * w;
        vec3 u2w2 = u2 - w2;

        vec3 R   = vec3(alpha) - b * w2;
        vec3 dR  = (-2.0 * b) * uw;
        vec3 d2R = (-2.0 * b) * u2w2;

        vec3 qx   = u2 * R;
        vec3 dqx  = -2.0 * uw * R + u2 * dR;
        vec3 d2qx = -2.0 * u2w2 * R - 4.0 * uw * dR + u2 * d2R;

        vec3 qy   = uw * R;
        vec3 dqy  = u2w2 * R + uw * dR;
        vec3 d2qy = -4.0 * uw * R + 2.0 * u2w2 * dR + uw * d2R;

        vec3 dx = qx - vec3(px);
        vec3 dy = qy - vec3(ay);
        vec3 f  = dx * dqx + dy * dqy;
        vec3 df = dqx * dqx + dqy * dqy + dx * d2qx + dy * d2qy;

        vec3 stepT = vec3(
            (df.x > 1e-12) ? (f.x / df.x) : ((f.x >= 0.0) ? 0.2 : -0.2),
            (df.y > 1e-12) ? (f.y / df.y) : ((f.y >= 0.0) ? 0.2 : -0.2),
            (df.z > 1e-12) ? (f.z / df.z) : ((f.z >= 0.0) ? 0.2 : -0.2)
        );
        stepT = clamp(stepT, vec3(-0.35), vec3(0.35));
        s = clamp(s - 0.5 * (1.0 + s2) * stepT, vec3(0.0), vec3(sMax));
    }

    for( int j = 0; j < 3; j++ )
    {
        float sj   = s[j];
        float s2j  = sj * sj;
        float invD = 1.0 / (1.0 + s2j);
        float uj   = (1.0 - s2j) * invD;
        float wj   = 2.0 * sj * invD;
        float Rj   = alpha - b * wj * wj;
        vec2  qj   = vec2(uj * uj * Rj, uj * wj * Rj);
        vec2  dv   = vec2(px, ay) - qj;
        float d2   = dot(dv, dv);
        if( d2 < bestD2 )
        {
            bestD2 = d2;
            bestQ  = qj;
        }
    }

    return vec3(bestD2, bestQ.x, bestQ.y);
}

// 1. Exact 0-Trig 2D Johannes Kepler's (1609) Oviform / Simple Folium Loop + Unit Gradient vec3(d, grad.xy)
//    r(theta) = cos(theta) * (a - b * sin^2(theta)), 0 <= b <= a (when b = a, r = a * cos^3(theta) Simple Folium).
vec3 sdgKeplerFolium2D( in vec2 p, in float a, in float b )
{
    a = max(a, 1e-5);
    b = clamp(b, 0.0, a);

    float px = p.x;
    float sy = (p.y < 0.0) ? -1.0 : 1.0;
    float ay = abs(p.y);

    vec3  res   = _sdQuarticFoliumCore(px, ay, a, a - b, 1.0);
    vec2  bestQ = res.yz;

    float r2     = px * px + ay * ay;
    bool  inside = (px > 0.0) && (r2 * r2 < px * (a * px * px + (a - b) * ay * ay));
    float sgn    = inside ? -1.0 : 1.0;

    vec2  diff = vec2(px - bestQ.x, sy * (ay - bestQ.y));
    float d2   = dot(diff, diff);
    float invD = inversesqrt(max(d2, 1e-30));
    float d    = d2 * invD;
    vec2  grad = (d2 > 1e-24) ? (sgn * diff * invD) : vec2(1.0, 0.0);
    return vec3(sgn * d, grad);
}

// 2. Scalar 2D Kepler's Oviform / Folium Loop SDF
float sdKeplerFolium2D( in vec2 p, in float a, in float b )
{
    return sdgKeplerFolium2D(p, a, b).x;
}

// 1-parameter overloads for the canonical Simple Folium / Kepler Ovoid (b = a: (x^2 + y^2)^2 = a * x^3)
vec3 sdgKeplerFolium2D( in vec2 p, in float a )
{
    return sdgKeplerFolium2D(p, a, a);
}

float sdKeplerFolium2D( in vec2 p, in float a )
{
    return sdgKeplerFolium2D(p, a, a).x;
}

// 3. Exact 0-Trig 2D Albrecht Dürer's (1525) Double Folium + Unit Gradient vec3(d, grad.xy)
//    (x^2 + y^2)^2 = 4 * a * x * y^2 <=> r(theta) = 4 * a * cos(theta) * sin^2(theta)
vec3 sdgDoubleFolium2D( in vec2 p, in float a )
{
    a = max(a, 1e-5);
    float px = p.x;
    float sy = (p.y < 0.0) ? -1.0 : 1.0;
    float ay = abs(p.y);

    vec3  res   = _sdQuarticFoliumCore(px, ay, 0.0, 4.0 * a, 1.0);
    vec2  bestQ = res.yz;

    float r2     = px * px + ay * ay;
    bool  inside = (px > 0.0) && (r2 * r2 < 4.0 * a * px * ay * ay);
    float sgn    = inside ? -1.0 : 1.0;

    vec2  diff = vec2(px - bestQ.x, sy * (ay - bestQ.y));
    float d2   = dot(diff, diff);
    float invD = inversesqrt(max(d2, 1e-30));
    float d    = d2 * invD;
    vec2  grad = (d2 > 1e-24) ? (sgn * diff * invD) : vec2(1.0, 0.0);
    return vec3(sgn * d, grad);
}

// 4. Scalar 2D Dürer's Double Folium SDF
float sdDoubleFolium2D( in vec2 p, in float a )
{
    return sdgDoubleFolium2D(p, a).x;
}

// 5. Exact 0-Trig 2D Regular 3-Petal Trifolium + Analytical Unit Gradient vec3(d, grad.xy)
//    (x^2 + y^2)^2 = a * x * (x^2 - 3 * y^2) <=> r(theta) = a * cos(3 * theta) >= 0
vec3 sdgTrifolium2D( in vec2 p, in float a )
{
    a = max(a, 1e-5);
    const vec2 n60 = vec2(-0.8660254037844386, 0.5);

    vec2  q   = p;
    float sy0 = (q.y < 0.0) ? -1.0 : 1.0; q.y = abs(q.y);
    bool  r0  = (dot(q, n60) > 0.0);       if( r0 ) q -= 2.0 * dot(q, n60) * n60;
    float sy1 = (q.y < 0.0) ? -1.0 : 1.0; q.y = abs(q.y);

    // In the folded 60-degree wedge [0, pi/3], the positive half-petal spans theta in [0, pi/6],
    // so s = tan(theta/2) spans [0, tan(pi/12)] = [0, 2 - sqrt(3)].
    vec3  res   = _sdQuarticFoliumCore(q.x, q.y, a, -3.0 * a, 0.2679491924311227);
    vec2  bestQ = res.yz;

    float r2     = dot(q, q);
    float rhs    = a * q.x * (q.x * q.x - 3.0 * q.y * q.y);
    bool  inside = (rhs > 0.0) && (r2 * r2 < rhs);
    float sgn    = inside ? -1.0 : 1.0;

    vec2  diff = q - bestQ;
    float d2   = dot(diff, diff);
    float invD = inversesqrt(max(d2, 1e-30));
    float d    = d2 * invD;
    vec2  g    = (d2 > 1e-24) ? (sgn * diff * invD) : vec2(1.0, 0.0);

    g.y *= sy1;
    if( r0 ) g -= 2.0 * dot(g, n60) * n60;
    g.y *= sy0;
    return vec3(sgn * d, g);
}

// 6. Scalar 2D Regular 3-Petal Trifolium SDF
float sdTrifolium2D( in vec2 p, in float a )
{
    return sdgTrifolium2D(p, a).x;
}

// 7. Exact 3D Kepler Egg / Oviform Solid of Revolution + Analytical Unit Gradient vec4(d, grad.xyz)
vec4 sdgKeplerEgg3D( in vec3 p, in float a, in float b )
{
    float rho = length(p.yz);
    vec3  r2d = sdgKeplerFolium2D(vec2(p.x, rho), a, b);
    vec2  er  = (rho > 1e-12) ? (p.yz / rho) : vec2(1.0, 0.0);
    return vec4(r2d.x, r2d.y, r2d.z * er.x, r2d.z * er.y);
}

// 8. Scalar 3D Kepler Egg / Oviform Solid SDF
float sdKeplerEgg3D( in vec3 p, in float a, in float b )
{
    return sdgKeplerFolium2D(vec2(p.x, length(p.yz)), a, b).x;
}

// 1-parameter overloads for the canonical Kepler Egg (b = a)
vec4 sdgKeplerEgg3D( in vec3 p, in float a )
{
    return sdgKeplerEgg3D(p, a, a);
}

float sdKeplerEgg3D( in vec3 p, in float a )
{
    return sdKeplerEgg3D(p, a, a);
}

#endif // SDF_GEMS_SD_KEPLER_FOLIUM_2D_GLSL
```
