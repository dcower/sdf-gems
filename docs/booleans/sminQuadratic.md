# `sminQuadratic` — $C^1$ Quadratic Polynomial Smooth Minimum, Material Blend & 3D Gradient

![sminQuadratic](../images/booleans/sminQuadratic.png)

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** — extended with Material Blend & 3D Gradient using **Gemini 4 Argon** (`sdf-gems`, 2026) |
| **Source** | [https://iquilezles.org/articles/smin/](https://iquilezles.org/articles/smin/) |
| **Category** | `booleans` |
| **GLSL File** | [`glsl/booleans/sminQuadratic.glsl`](../../glsl/booleans/sminQuadratic.glsl) |
| **Exact Euclidean SDF?** | Conservative (`||grad|| <= 1`, $C^1$ continuous) |

---

## Mathematical Overview

$C^1$ quadratic Closest-Distance (CD) polynomial smooth minimum (`sminQuadratic`), smooth maximum (`smaxQuadratic`), smooth subtraction (`opSmoothSubtractionQuadratic`), material blend weight (`sminQuadraticBlend`), and analytical 3D gradient (`sdgSmoothQuadratic`, `2.80e-9` max central-FD error).

---

## GLSL Implementation

```glsl
float sminQuadratic( float a, float b, float k )
{
    k *= 4.0;
    float h = max( k - abs(a - b), 0.0 ) / k;
    return min(a, b) - h * h * k * (1.0 / 4.0);
}

vec2 sminQuadraticBlend( float a, float b, float k )
{
    float h = max( 4.0 * k - abs(a - b), 0.0 ) / (4.0 * k);
    float w = h * h;
    return vec2( min(a, b) - w * k, (a < b) ? 0.5 * w : 1.0 - 0.5 * w );
}

vec4 sdgSmoothQuadratic( in vec4 a, in vec4 b, in float k )
{
    k *= 4.0;
    float h = max( k - abs(a.x - b.x), 0.0 ) / (2.0 * k);
    return vec4( min(a.x, b.x) - h * h * k, mix(a.yzw, b.yzw, (a.x < b.x) ? h : 1.0 - h) );
}
```
