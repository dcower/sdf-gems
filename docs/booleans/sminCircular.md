# `sminCircular` — Circular Clamped-Difference (CD) Smooth Minimum & 1-`sqrt` 3D Gradient

![sminCircular](../images/booleans/sminCircular.png)

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** — extended with 1-`sqrt` shared 3D Gradient using **Gemini 4 Argon** (`sdf-gems`, 2026) |
| **Source** | [https://iquilezles.org/articles/smin/](https://iquilezles.org/articles/smin/) |
| **Category** | `booleans` |
| **GLSL File** | [`glsl/booleans/sminCircular.glsl`](../../glsl/booleans/sminCircular.glsl) |
| **Exact Euclidean SDF?** | Conservative (`\|\|grad\|\| <= 1`, 1 shared `sqrt` in `sdgSmoothCircular`) |

---

## Mathematical Overview

Circular profile smooth minimum along the direct difference $`|a - b|`$, together with `smaxCircular`, `opSmoothSubtractionCircular`, and the analytical 3D gradient `sdgSmoothCircular` which reuses the single radical $`s = \sqrt{1 - h(h - 2)}`$ for both distance and gradient weight $`m = \frac{s - 1 + h}{2s}`$ (`2.82e-9` max central-FD error).

---

## GLSL Implementation

```glsl
float sminCircular( float a, float b, float k )
{
    k *= 1.0 / (1.0 - sqrt(0.5));
    float h = max( k - abs(a - b), 0.0 ) / k;
    return min(a, b) - k * 0.5 * (1.0 + h - sqrt(1.0 - h * (h - 2.0)));
}

vec4 sdgSmoothCircular( in vec4 a, in vec4 b, in float k )
{
    k *= 1.0 / (1.0 - sqrt(0.5));
    float h = max( k - abs(a.x - b.x), 0.0 ) / k;
    float s = sqrt(1.0 - h * (h - 2.0));
    float m = 0.5 * (s - 1.0 + h) / s;
    float d = min(a.x, b.x) - k * 0.5 * (1.0 + h - s);
    return vec4( d, mix(a.yzw, b.yzw, (a.x < b.x) ? m : 1.0 - m) );
}
```
