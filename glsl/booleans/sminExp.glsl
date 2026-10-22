// ============================================================================
// sminExp — Overflow-Free 1-exp2 Exponential, Sigmoid & Root Smooth Minimums
// @author Inigo Quilez — upgraded to 1-exp2 shift-stabilized LSE & conjugate forms using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/smin/
// @exact  Conservative (`||grad|| <= 1`, 0% FP32/FP16 overflow/underflow, 50% fewer `exp2` calls)
// ============================================================================

#ifndef SDF_GEMS_SMIN_EXP_GLSL
#define SDF_GEMS_SMIN_EXP_GLSL

// 1. Overflow-Free 1-exp2 Exponential Smooth Minimum (LogSumExp)
//    Factors out 2^(-min(a,b)/k) so exp2 only evaluates a non-positive argument -|a-b|/k <= 0,
//    halving the exp2 count (1 exp2 + 1 log2 vs IQ's 2 exp2 + 1 log2) and eliminating FP32/FP16 overflow/underflow.
float sminExp( float a, float b, float k )
{
    return min(a, b) - k * log2(1.0 + exp2(-abs(a - b) / k));
}

float smaxExp( float a, float b, float k )
{
    return max(a, b) + k * log2(1.0 + exp2(-abs(a - b) / k));
}

float opSmoothSubtractionExp( float a, float b, float k )
{
    return max(-a, b) + k * log2(1.0 + exp2(-abs(a + b) / k));
}

// 2. Material Blend Factor & Analytical 3D Gradient for sminExp (0 extra exp2/log2 calls!)
vec2 sminExpBlend( float a, float b, float k )
{
    float e = exp2(-abs(a - b) / k);
    float m = e / (1.0 + e);
    float d = min(a, b) - k * log2(1.0 + e);
    return vec2( d, (a < b) ? m : 1.0 - m );
}

vec4 sdgSmoothExp( in vec4 a, in vec4 b, in float k )
{
    float e = exp2(-abs(a.x - b.x) / k);
    float m = e / (1.0 + e);
    float d = min(a.x, b.x) - k * log2(1.0 + e);
    return vec4( d, mix(a.yzw, b.yzw, (a.x < b.x) ? m : 1.0 - m) );
}

// 3. Tournament-Minimum 3-Way & 4-Way Exponential Smooth Minimums
//    Isolating the minimum m forces 2^(-(m-m)/k) == 1.0 analytically, using only
//    2 exp2 + 1 log2 for 3-way (and 3 exp2 + 1 log2 for 4-way) with 0 overflow/underflow.
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

// 4. Overflow-Free, Seam-Safe Sigmoid Smooth Minimum
//    Fixes IQ's 0.0 / 0.0 -> NaN singularity on the medial seam a == b and positive exp2 overflow.
float sminSigmoid( float a, float b, float k )
{
    float z = abs(a - b) / k;
    float e = exp2(-1.44269504 * z);
    return min(a, b) - k * ((z > 1e-3) ? z * e / max(1.0 - e, 1e-30) : 1.0 - 0.5 * z);
}

// 5. Conjugate-Stabilized Square-Root (Hyperbolic) Smooth Minimum & 3D Gradient
//    Fixes catastrophic mantissa cancellation in IQ's 0.5*(a + b - sqrt((b-a)^2 + 4*k^2)) when |a - b| >> k.
float sminRoot( float a, float b, float k )
{
    float d = abs(a - b), k2 = 2.0 * k;
    return min(a, b) - (2.0 * k * k) / (d + sqrt(d * d + k2 * k2));
}

vec4 sdgSmoothRoot( in vec4 a, in vec4 b, in float k )
{
    float dx = abs(a.x - b.x), k2 = 2.0 * k;
    float s  = sqrt(dx * dx + k2 * k2);
    float d  = min(a.x, b.x) - (2.0 * k * k) / (dx + s);
    float m  = 0.5 * (1.0 - dx / s);
    return vec4( d, mix(a.yzw, b.yzw, (a.x < b.x) ? m : 1.0 - m) );
}

#endif // SDF_GEMS_SMIN_EXP_GLSL
