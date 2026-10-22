// ============================================================================
// sminCubic — C2 Cubic & C2-Symmetric Quartic Polynomial Smooth Minimums & 3D Gradients
// @author Inigo Quilez — extended with Quartic & 3D Gradients using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/smin/
// @exact  Conservative (`||grad|| <= 1`, C2 at blend boundaries; sminQuartic is C^inf across a = b)
// ============================================================================

#ifndef SDF_GEMS_SMIN_CUBIC_GLSL
#define SDF_GEMS_SMIN_CUBIC_GLSL

// 1. C2 Cubic Polynomial Smooth Minimum
float sminCubic( float a, float b, float k )
{
    k *= 6.0;
    float h = max( k - abs(a - b), 0.0 ) / k;
    return min(a, b) - h * h * h * k * (1.0 / 6.0);
}

float smaxCubic( float a, float b, float k )
{
    return -sminCubic(-a, -b, k);
}

float opSmoothSubtractionCubic( float a, float b, float k )
{
    return -sminCubic(a, -b, k);
}

// 2. Cubic Smooth Minimum with Material Mix Factor: vec2(d, blend_weight_in_0_1)
vec2 sminCubicBlend( float a, float b, float k )
{
    float h = max( 6.0 * k - abs(a - b), 0.0 ) / (6.0 * k);
    float w = h * h * h;
    return vec2( min(a, b) - w * k, (a < b) ? 0.5 * w : 1.0 - 0.5 * w );
}

// 3. Cubic Smooth Minimum with Analytical 3D Gradient: vec4(d, gx, gy, gz)
vec4 sdgSmoothCubic( in vec4 a, in vec4 b, in float k )
{
    k *= 6.0;
    float h = max( k - abs(a.x - b.x), 0.0 ) / k;
    float m = 0.5 * h * h;
    return vec4( min(a.x, b.x) - h * h * h * k * (1.0 / 6.0),
                 mix(a.yzw, b.yzw, (a.x < b.x) ? m : 1.0 - m) );
}

// 4. C2-Symmetric Quartic Polynomial Smooth Minimum (eliminates the |a - b|^3 third-derivative cusp at a = b)
float sminQuartic( float a, float b, float k )
{
    k *= 16.0 / 3.0;
    float h = max( k - abs(a - b), 0.0 ) / k;
    return min(a, b) - h * h * h * (4.0 - h) * k * (1.0 / 16.0);
}

float smaxQuartic( float a, float b, float k )
{
    return -sminQuartic(-a, -b, k);
}

vec4 sdgSmoothQuartic( in vec4 a, in vec4 b, in float k )
{
    k *= 16.0 / 3.0;
    float h = max( k - abs(a.x - b.x), 0.0 ) / k;
    float m = 0.25 * h * h * (3.0 - h);
    float d = min(a.x, b.x) - h * h * h * (4.0 - h) * k * (1.0 / 16.0);
    return vec4( d, mix(a.yzw, b.yzw, (a.x < b.x) ? m : 1.0 - m) );
}

#endif // SDF_GEMS_SMIN_CUBIC_GLSL
