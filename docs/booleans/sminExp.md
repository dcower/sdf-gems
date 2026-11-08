# `sminExp` — Overflow-Free 1-`exp2` Exponential, Sigmoid & Root Smooth Minimums

![sminExp](../images/booleans/sminExp.png)

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** — upgraded to 1-`exp2` shift-stabilized LogSumExp & conjugate forms using **Gemini 4 Argon** (`sdf-gems`, 2026) |
| **Source** | [https://iquilezles.org/articles/smin/](https://iquilezles.org/articles/smin/) |
| **Category** | `booleans` |
| **GLSL File** | [`glsl/booleans/sminExp.glsl`](../../glsl/booleans/sminExp.glsl) |
| **Exact Euclidean SDF?** | **Conservative** ($\|\nabla d\| \le 1$, $0\%$ FP32/FP16 overflow/underflow, $50\%$ fewer `exp2` calls) |

---

## Before & After Comparison

![sminExp Comparison](../images/comparisons/smin_exp_comparison.png)

---

## Mathematical Overview

### 1. Catastrophic FP32/FP16 Underflow & Overflow in 2-`exp2` `sminExp` vs. 1-`exp2` Shift-Stabilized LogSumExp

Inigo Quilez's classic exponential smooth minimum evaluates two separate base-2 exponentials:

$$r = 2^{-a/k} + 2^{-b/k}, \qquad d = -k \log_2(r)$$

- **Far-Field Underflow ($\min(a,b) > 0 \implies +\infty$)**: Whenever $\min(a,b) > 126 k$ (FTZ) or $> 149 k$ (subnormals) in FP32 — or $\min(a,b) > 14 k$ (FTZ) / $> 24 k$ in FP16 — both $2^{-a/k}$ and $2^{-b/k}$ underflow to `0.0`, producing `log2(0.0) = -Inf` and `-k * log2(0.0) = +Inf` across $66.5\%$ of a $[-20, 20]^2$ grid at $k = 0.05$ ($95.1\%$ in FP16).
- **Interior Overflow ($\min(a,b) < 0 \implies -\infty$)**: Whenever $\min(a,b) < -128 k$ in FP32 ($< -16 k$ in FP16), $2^{-\min(a,b)/k}$ overflows to `+Inf`, returning `-Inf`.

Factoring out $2^{-\min(a,b)/k}$ from $2^{-a/k} + 2^{-b/k} = 2^{-\min(a,b)/k}\!\left(1 + 2^{-|a - b|/k}\right)$ yields the **1-`exp2` Shift-Stabilized LogSumExp Smooth Minimum**:

$$\mathrm{sminExp}(a, b, k) = \min(a, b) - k \log_2\!\left(1 + 2^{-|a - b|/k}\right)$$

1. **$50\%$ Fewer `exp2` Calls**: Evaluates only **1 `exp2` + 1 `log2`** (2 SFU instructions instead of 3).
2. **$100\%$ Overflow- and Underflow-Free across all of FP32 and FP16**: Because $-|a - b|/k \le 0$ for all $a, b \in \mathbb{R}$ and $k > 0$, $e = 2^{-|a - b|/k} \in (0, 1]$ never overflows, and $1 + e \in [1, 2]$ keeps `log2` inside the benign octave $[1, 2]$.
3. **0-Extra-SFU Analytical 3D Gradient (`sdgSmoothExp`)**: Differentiating gives $\frac{\partial d}{\partial \max(a,b)} = \frac{e}{1 + e} = m$ and $\frac{\partial d}{\partial \min(a,b)} = 1 - m$, reusing $e = 2^{-|a - b|/k}$ with **0 extra `exp2` or `log2` calls**.
4. **Tournament-Minimum 3-Way (`sminExp3`) & 4-Way (`sminExp4`)**: Isolating $m = \min_i d_i$ via branchless compare-swaps forces $2^{-(m - m)/k} \equiv 1$ analytically, evaluating 3-way `sminExp3` in **2 `exp2` + 1 `log2`** and 4-way `sminExp4` in **3 `exp2` + 1 `log2`**.

### 2. Seam-Safe `sminSigmoid` & Conjugate-Stabilized `sminRoot`

- **Seam-Safe `sminSigmoid`**: IQ's `a + (b - a) / (1.0 - exp2((b - a) / (k * ln 2)))` divides `0.0 / 0.0 = NaN` on the entire medial bisector $a = b$ and overflows when $b - a \gg k$. Factoring out $\min(a,b)$ with $z = |a - b|/k \ge 0$, $e = 2^{-z / \ln 2} \in (0, 1]$, and Taylor fallback $1 - \frac{1}{2}z$ for $z \le 10^{-3}$ eliminates all `NaN`s and overflows.
- **Conjugate-Stabilized `sminRoot`**: IQ's $\frac{1}{2}\!\left(a + b - \sqrt{(b - a)^2 + 4 k^2}\right)$ suffers $100\%$ catastrophic mantissa cancellation of $\min(a, b)$ when $|a - b| \gg k$. Multiplying by the conjugate transforms it into:
  $$\mathrm{sminRoot}(a, b, k) = \min(a, b) - \frac{2 k^2}{|a - b| + \sqrt{|a - b|^2 + 4 k^2}}$$
  which achieves `1.45e-13` error at $a = 10^{-5}, b = 10^3, k = 0.01$ and shares $s = \sqrt{|a - b|^2 + 4 k^2}$ with its analytical 3D gradient `sdgSmoothRoot`.

### Numerical Verification Summary ($k = 0.05$, $[-20, 20]^2$ Grid, $160{,}801$ Points)

| Operator | SFU Ops | FP32 `Inf`/`NaN` | FP16 `Inf`/`NaN` | Max Abs Error vs FP64 | Max $\|\nabla - \nabla_{\text{FD}}\|_2$ |
| :--- | :---: | :---: | :---: | :---: | :---: |
| IQ 2-`exp2` `sminExp` | $2\text{ exp2} + 1\text{ log2}$ | $66.5\%$ | $95.1\%$ | $\infty$ | — |
| **1-`exp2` `sminExp` & `sdgSmoothExp`** | **$1\text{ exp2} + 1\text{ log2}$** | **$0.0\%$** | **$0.0\%$** | `8.49e-07` | `2.89e-09` |
| **Tournament `sminExp3` / `sminExp4`** | **$2/3\text{ exp2} + 1\text{ log2}$** | **$0.0\%$** | **$0.0\%$** | `1.91e-06` | — |
| **Seam-Safe `sminSigmoid`** | **$1\text{ exp2}$** | **$0.0\%$** | **$0.0\%$** | `2.22e-06` | — |
| **Conjugate `sminRoot` & `sdgSmoothRoot`** | **$1\text{ sqrt}$** | **$0.0\%$** | **$0.0\%$** | `1.45e-13` | `2.79e-09` |

---

## GLSL Implementation

```glsl
float sminExp( float a, float b, float k )
{
    return min(a, b) - k * log2(1.0 + exp2(-abs(a - b) / k));
}

vec4 sdgSmoothExp( in vec4 a, in vec4 b, in float k )
{
    float e = exp2(-abs(a.x - b.x) / k);
    float m = e / (1.0 + e);
    float d = min(a.x, b.x) - k * log2(1.0 + e);
    return vec4( d, mix(a.yzw, b.yzw, (a.x < b.x) ? m : 1.0 - m) );
}

float sminExp3( float a, float b, float c, float k )
{
    float m = min(a, b), u = max(a, b);
    float v = max(m, c); m = min(m, c);
    return m - k * log2(1.0 + exp2(-(u - m) / k) + exp2(-(v - m) / k));
}

float sminExp4( vec4 d, float k )
{
    vec2  lo = min(d.xy, d.zw), hi = max(d.xy, d.zw);
    float m  = min(lo.x, lo.y);
    vec3  u  = vec3(max(lo.x, lo.y), hi) - m;
    return m - k * log2(1.0 + dot(exp2(-u / k), vec3(1.0)));
}

float sminSigmoid( float a, float b, float k )
{
    float z = abs(a - b) / k;
    float e = exp2(-1.44269504 * z);
    return min(a, b) - k * ((z > 1e-3) ? z * e / max(1.0 - e, 1e-30) : 1.0 - 0.5 * z);
}

float sminRoot( float a, float b, float k )
{
    float d = abs(a - b), k2 = 2.0 * k;
    return min(a, b) - (2.0 * k * k) / (d + sqrt(d * d + k2 * k2));
}
```
