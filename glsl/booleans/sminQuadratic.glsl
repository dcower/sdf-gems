// ============================================================================
// sminQuadratic — C1 Quadratic Polynomial Smooth Minimum, Blend & 3D Gradient
// @author Inigo Quilez — extended with Material Blend & 3D Gradient using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/smin/
// @exact  Conservative (`||grad|| <= 1`, C1 continuous)
// ============================================================================

#ifndef SDF_GEMS_SMIN_QUADRATIC_GLSL
#define SDF_GEMS_SMIN_QUADRATIC_GLSL

// 1. C1 Quadratic Polynomial Smooth Minimum
float sminQuadratic( float a, float b, float k )
{
    k *= 4.0;
    float h = max( k - abs(a - b), 0.0 ) / k;
    return min(a, b) - h * h * k * (1.0 / 4.0);
}

float smaxQuadratic( float a, float b, float k )
{
    return -sminQuadratic(-a, -b, k);
}

float opSmoothSubtractionQuadratic( float a, float b, float k )
{
    return -sminQuadratic(a, -b, k);
}

// 2. Quadratic Smooth Minimum with Material Mix Factor: vec2(d, blend_weight_in_0_1)
vec2 sminQuadraticBlend( float a, float b, float k )
{
    float h = max( 4.0 * k - abs(a - b), 0.0 ) / (4.0 * k);
    float w = h * h;
    return vec2( min(a, b) - w * k, (a < b) ? 0.5 * w : 1.0 - 0.5 * w );
}

// 3. Quadratic Smooth Minimum with Analytical 3D Gradient: vec4(d, gx, gy, gz)
vec4 sdgSmoothQuadratic( in vec4 a, in vec4 b, in float k )
{
    k *= 4.0;
    float h = max( k - abs(a.x - b.x), 0.0 ) / (2.0 * k);
    return vec4( min(a.x, b.x) - h * h * k, mix(a.yzw, b.yzw, (a.x < b.x) ? h : 1.0 - h) );
}

#endif // SDF_GEMS_SMIN_QUADRATIC_GLSL
