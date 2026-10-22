// ============================================================================
// sminCircular — Circular Clamped-Difference (CD) Smooth Minimum, Blend & 3D Gradient
// @author Inigo Quilez — extended with 1-sqrt shared 3D Gradient using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/smin/
// @exact  Conservative (`||grad|| <= 1`, exact circular profile for orthogonal planes, 1 shared `sqrt`)
// ============================================================================

#ifndef SDF_GEMS_SMIN_CIRCULAR_GLSL
#define SDF_GEMS_SMIN_CIRCULAR_GLSL

// 1. Circular Clamped-Difference (CD) Smooth Minimum
float sminCircular( float a, float b, float k )
{
    k *= 1.0 / (1.0 - sqrt(0.5));
    float h = max( k - abs(a - b), 0.0 ) / k;
    return min(a, b) - k * 0.5 * (1.0 + h - sqrt(1.0 - h * (h - 2.0)));
}

float smaxCircular( float a, float b, float k )
{
    return -sminCircular(-a, -b, k);
}

float opSmoothSubtractionCircular( float a, float b, float k )
{
    return -sminCircular(a, -b, k);
}

// 2. Circular Smooth Minimum with Analytical 3D Gradient (reuses the single sqrt s = sqrt(1 - h*(h - 2)))
vec4 sdgSmoothCircular( in vec4 a, in vec4 b, in float k )
{
    k *= 1.0 / (1.0 - sqrt(0.5));
    float h = max( k - abs(a.x - b.x), 0.0 ) / k;
    float s = sqrt(1.0 - h * (h - 2.0));
    float m = 0.5 * (s - 1.0 + h) / s;
    float d = min(a.x, b.x) - k * 0.5 * (1.0 + h - s);
    return vec4( d, mix(a.yzw, b.yzw, (a.x < b.x) ? m : 1.0 - m) );
}

#endif // SDF_GEMS_SMIN_CIRCULAR_GLSL
