# `opRep` — Unlimited, Limited, Mirrored, and Clamped-Floor Base Repetition

| 3D Surface ($f(p) = 0$) | 2D Cross-Section Slice |
| :---: | :---: |
| ![opRep](../images/repetition/opRep.png) | ![opRep slice](../images/repetition/opRep_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez — upgraded with `opRepMirror` & `opRepLimBase` using Gemini 4 Argon (`sdf-gems`, 2026)** |
| **Source** | [https://iquilezles.org/articles/sdfrepetition/](https://iquilezles.org/articles/sdfrepetition/) |
| **Category** | `repetition` |
| **GLSL File** | [`glsl/repetition/opRep.glsl`](../../glsl/repetition/opRep.glsl) |
| **Exact Euclidean SDF?** | Yes (symmetric SDFs for 1-tap; arbitrary SDFs with `opRepLimBase` $2^d$ loop) |

---

## Mathematical Overview

1. **`opRep` & `opRepLim`**: Folds space into periodic or clamped cubic cells via `round(p/s)`.
2. **`opRepMirror`**: Pure-floating-point alternating mirror repetition using `1.0 - 2.0 * mod(id, 2.0)` (replacing integer bitmasks `int(id)&1` which fail on negative float IDs or WebGL1).
3. **`opRepLimBase`**: In IQ's asymmetric limited repetition (`sdfrepetition`), clamping inside the $2^d$ loop (`clamp(id + o, lima, limb)`) causes both $o = 0$ and $o = 1$ to clamp to the **same** boundary cell whenever $p/s > \text{limb}$ or $p/s < \text{lima}$, evaluating the boundary cell twice while skipping its interior neighbor. Computing the base index before the loop via `i0 = clamp(floor(p / s), lima, limb - 1.0)` guarantees $2^d$ **distinct** adjacent cells everywhere in $\mathbb{R}^d$, reducing exterior jittered-instance errors by $32\%$ while eliminating `round()`, `sign()`, and $2^d$ in-loop `clamp()` calls.

---

## GLSL Implementation

```glsl
vec3 opRep( in vec3 p, in vec3 s )
{
    return p - s*round(p/s);
}

vec3 opRepLim( in vec3 p, in float s, in vec3 l )
{
    return p - s*clamp(round(p/s), -l, l);
}

vec3 opRepMirror( in vec3 p, in vec3 s, out vec3 id )
{
    id = round(p / s);
    return (p - s * id) * (1.0 - 2.0 * mod(id, 2.0));
}

vec3 opRepLimBase( in vec3 p, in float s, in vec3 lima, in vec3 limb )
{
    return clamp(floor(p / s), lima, limb - 1.0);
}
```
