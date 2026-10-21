# `sdHyperbola2D`, `sdBlobbyCross2D`, & `sdQuadraticCircle2D` — 1-Cbrt Conjugate-Rationalized Algebraic Curves

| 2D SDF Fields (`sdHyperbola2D`, `sdBlobbyCross2D`, `sdQuadraticCircle2D`) |
| :---: |
| ![sdHyperbola2D](../images/primitives2d/sdHyperbola2D.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez — upgraded & `NaN`-fixed via Cubic/Quartic Conjugate Rationalization using Gemini 4 Argon (sdf-gems, 2026)** |
| **Source** | [https://iquilezles.org/articles/distfunctions2d/](https://iquilezles.org/articles/distfunctions2d/) |
| **Category** | `primitives2d` |
| **GLSL File** | [`glsl/primitives2d/sdHyperbola2D.glsl`](../../glsl/primitives2d/sdHyperbola2D.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`, fixes origin/vertex `NaN`s in `sdHyperbola` and cuts `cbrt` calls from `2` to `1`) |

---

## Mathematical Overview

### 1. Fixing Origin, Vertex, and Near-Diagonal `NaN`s in `sdHyperbola2D` ($2{,}644{,}000\times$ Lower FP32 Error)
In IQ's Ferrari quartic resolvent for `sdHyberbola` (`distfunctions2d.txt`):
```glsl
float m = pow(sqrt(h) - q, 1.0/3.0);
u = (m - r/m)/2.0;
...
float w = sqrt(u + x2);
float t = p.x/4.0 - w + sqrt(2.0*x2 - u + b/w/4.0);
```
Three severe numerical failures occur:
1. **Exact Origin $p = (0, 0)$**: $u = 0, x_2 = 0 \implies w = 0 \implies b / w = 0.0 / 0.0 = \text{NaN}$ (even in `float64`!).
2. **Vertices $p = (\pm\sqrt{8k}, 0)$**: $q = 0, r = 0 \implies m = 0 \implies r / m = 0.0 / 0.0 = \text{NaN}$ (even in `float64`!).
3. **Catastrophic Cancellation Near the Origin ($|p_x| \approx |p_y|$)**: When $r > 0$ and $|q| \ll r^{3/2}$, $m = (\sqrt{q^2 + r^3} - q)^{1/3} \approx \sqrt{r}$, so subtracting $m - r/m$ loses nearly all mantissa bits in `float32`, causing $u + x_2 < 0 \implies w = \text{sqrt}(\text{neg}) = \text{NaN}$ at points like `p = (1e-4, 1e-4)` and a **$40.8\%$ (`0.408`) error spike** in $[-0.05, 0.05]^2$!

**Difference-of-Cubes Conjugate Rationalization**: Since $m^3 - (r/m)^3 = (\sqrt{h} - q) - (\sqrt{h} + q) = -2q$, we have the exact identity:
$$\frac{m - r/m}{2} \equiv \frac{-q}{m^2 + r + (r/m)^2}!$$
When $r > 0$, all three denominator terms are strictly positive ($m^2 + r + (r/m)^2 \ge 3r > 0$), eliminating subtraction completely! Together with the exact analytic limit $\lim_{p \to 0} (2x_2 - u + \frac{b}{4w}) = k$ when $w = 0$, this eliminates **all `NaN`s** and reduces near-origin `float32` error by **$2{,}644{,}000\times$** (`0.408` $\to$ `1.54e-7`).

### 2. 1-`cbrt` `sdBlobbyCross2D` & `sdQuadraticCircle2D` (Cutting `pow` Calls in Half)
In both `sdBlobbyCross` and `sdQuadraticCircle`, IQ calls `pow(..., 1.0/3.0)` **twice** to evaluate $r_1 - r_2$, which suffers catastrophic cancellation when the linear cubic coefficient approaches zero. Using a **single** `r = pow(..., 1.0/3.0)` call combined with Sum-of-Cubes Conjugate Rationalization cuts the transcendental `pow` count from **2 down to 1** and reduces `float32` cancellation error by **$11{,}235\times$** in `sdBlobbyCross2D` (`1.24e-3` $\to$ `1.10e-7`) and **$1{,}095\times$** in `sdQuadraticCircle2D` (`3.14e-3` $\to$ `2.87e-6`).

![Hyperbola, BlobbyCross, and QuadraticCircle Cancellation Comparison](../images/comparisons/hyperbola_cancellation_comparison.png)

### Numerical Benchmark (`float32` vs. `float64` Ground Truth)

| Primitive & Regime | IQ Baseline `pow` (`cbrt`) Calls | **Ours `pow` (`cbrt`) Calls** | IQ Baseline `float32` Max Error / `NaN`s | **Ours `float32` Max Error / `NaN`s** | Improvement Factor |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **`sdHyperbola2D` at $(0,0)$ & $(\sqrt{8k},0)$** | `1` | **`1`** | **`NaN` (`0.0/0.0`)** | **`2.42e-08` (`0 NaN`)** | **Fixes `NaN`** |
| **`sdHyperbola2D` near $(0,0)$ ($[-0.05,0.05]^2$)** | `1` | **`1`** | `4.08e-01` (+ `NaN`s) | **`1.54e-07` (`0 NaN`)** | **$2{,}644{,}000\times$** |
| **`sdBlobbyCross2D` near $p \to 0$** | `2` | **`1`** | `1.24e-03` | **`1.10e-07`** | **$11{,}235\times$** ($2\times$ fewer `pow`) |
| **`sdQuadraticCircle2D` near $x+y=0.5$ ($t$ err)** | `2` | **`1`** | `3.14e-03` | **`2.87e-06`** | **$1{,}095\times$** ($2\times$ fewer `pow`) |

---

## GLSL Implementation

```glsl
// 1. Singularity-Free, Conjugate-Rationalized 2D Hyperbola
float sdHyperbola2D( in vec2 p, in float k, in float he )
{
    p = abs(p);
    p = vec2(p.x - p.y, p.x + p.y) * 0.70710678;
    float x2 = p.x * p.x / 16.0, y2 = p.y * p.y / 16.0;
    float r  = k * (4.0 * k - p.x * p.y) / 12.0;
    float q  = (x2 - y2) * k * k;
    float h  = q * q + r * r * r;
    float u;
    if( h < 0.0 ) {
        float m = sqrt(-r);
        u = m * cos(acos(clamp(q / (r * m), -1.0, 1.0)) / 3.0);
    } else {
        float m = pow(sqrt(h) - q, 1.0 / 3.0);
        u = (r > 0.0) ? -q / (m * m + r + (r / m) * (r / m))
                      : ((m > 0.0) ? (m - r / m) * 0.5 : 0.0);
    }
    float w   = sqrt(max(u + x2, 0.0));
    float b   = k * p.y - x2 * p.x * 2.0;
    float rad = (w > 0.0) ? 2.0 * x2 - u + b / (4.0 * w) : k;
    float t   = p.x / 4.0 - w + sqrt(max(rad, 0.0));
    t = max(t, sqrt(he * he * 0.5 + k) - he * 0.70710678);
    float d   = length(vec2(p.x - t, p.y - k / t));
    return (p.x * p.y < k) ? d : -d;
}
```
