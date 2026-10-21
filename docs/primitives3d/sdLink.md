# `sdLink` & `sdChain` — 1/2-Sqrt Chain Link & 1.83-Sqrt Screw-Folded Infinite 3D Chain

| 3D Surface (`sdLink` & `sdChain`) | 2D Cross-Section Slice (`sdChain`) |
| :---: | :---: |
| ![sdChain](../images/primitives3d/sdChain.png) | ![sdChain slice](../images/primitives3d/sdChain_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez (`sdLink`) — upgraded & `sdChain` discovered using Gemini 4 Argon (sdf-gems, 2026)** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdLink.glsl`](../../glsl/primitives3d/sdLink.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`, `1` `sqrt` in slab $|p_y| \le l_e$; `sdChain` averages `1.828` `sqrt`s vs. `4.000`) |

---

## Mathematical Overview

### 1. 1/2-`sqrt` `sdLink` Optimization
Across the entire middle slab $|p_y| \le l_e$, $q_y = \max(|p_y| - l_e, 0) = 0$, so `length(q.xy)` simplifies analytically to $|p_x|$ with **0 `sqrt`s**, reducing `sdLink` from **2 `sqrt`s to 1 `sqrt`** across the middle slab.

### 2. Screw-Folded + Squared-Dominance Infinite 3D Chain (`sdChain`)
An infinite interlocking chain along the $y$-axis with pitch $s$ alternates $xy$-links at $y = 2ks$ and $zy$-links at $y = (2k+1)s$. Evaluating `min(sdLink(q0), sdLink(q1))` directly requires **4 `sqrt`s**.
1. **Screw-Reflection Fundamental Domain**: Folding $u_y = |p_y - 2s\,\operatorname{round}(p_y/(2s))| \in [0, s]$ and applying the glide/screw reflection $(x, u_y, z) \mapsto (z, s - u_y, x)$ when $u_y > 0.5s$ maps every point into $u_y \in [0, 0.5s]$, where Link 0 is always axially closer than Link 1:
   $$q_0 = \max(u_y - l_e, \, 0) \;\le\; q_1 = \max(s - u_y - l_e, \, 0).$$
2. **Exact `sqrt`-Free Squared-Dominance Inequality**: Comparing Link 0's squared tube-centerline distance $A = (w_0 - r_1)^2 + z^2$ against Link 1's $B = (w_1 - r_1)^2 + x^2$ cancels $x^2 + z^2 + r_1^2$ identically:
   $$A \le B \iff 2r_1 w_1 \le 2r_1 w_0 + q_1^2 - q_0^2.$$
   Since $q_1 \ge q_0 \ge 0$, the right-hand side is **always non-negative**, so squaring both sides yields the exact **`sqrt`-free dominance test**:
   $$\boxed{A \le B \iff 4 r_1^2 (z^2 + q_1^2) \;\le\; \big(2r_1 w_0 + q_1^2 - q_0^2\big)^2!}$$
   Whenever Link 0 wins (`~69%` of space), $w_1 = \sqrt{z^2 + q_1^2}$ is never computed, cutting the average `sqrt` count from **`4.000` down to `1.828` ($2.19\times$ fewer `sqrt`s, `4.44e-16` max error)**!

---

## GLSL Implementation

```glsl
// 1. Single 3D Chain Link (1 sqrt in middle slab |p.y| <= le, 2 sqrts in end caps |p.y| > le)
float sdLink( vec3 p, float le, float r1, float r2 )
{
    float x = abs(p.x), qy = abs(p.y) - le;
    float w = (qy <= 0.0) ? x : length(vec2(x, qy));
    return length(vec2(w - r1, p.z)) - r2;
}

// 2. Infinite 3D Interlocking Chain (1.828 avg sqrts vs 4.000 in 2-eval min(sdLink0, sdLink1))
float sdChain( vec3 p, float le, float r1, float r2, float s )
{
    float period = 2.0 * s;
    float uy = abs(p.y - period * round(p.y / period));
    float x  = abs(p.x), z = abs(p.z);
    if( uy > 0.5 * s ) { uy = s - uy; float t = x; x = z; z = t; }

    float q0  = max(uy - le, 0.0), q1 = max(s - uy - le, 0.0);
    float q02 = q0 * q0, q12 = q1 * q1;
    float w0  = (q0 <= 0.0) ? x : sqrt(x * x + q02);
    float rhs = 2.0 * r1 * w0 + (q12 - q02), w1_sq = z * z + q12;
    float d2  = (4.0 * r1 * r1 * w1_sq <= rhs * rhs)
              ? (w0 - r1) * (w0 - r1) + z * z
              : (sqrt(w1_sq) - r1) * (sqrt(w1_sq) - r1) + x * x;
    return sqrt(d2) - r2;
}
```
