// ============================================================================
// sminQuintic — C3/C4 Quintic Polynomial Smooth Minimum & 3D Gradient
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/smin/
// @exact  Conservative (`||grad|| <= 1`, C3 at blend boundaries, C4 at medial seam a = b)
// ============================================================================

#ifndef SDF_GEMS_SMIN_QUINTIC_GLSL
#define SDF_GEMS_SMIN_QUINTIC_GLSL

// 1. C3/C4 Quintic Polynomial Smooth Minimum (cancels both |a-b|^1 and |a-b|^3 at a = b)
float sminQuintic( float a, float b, float k )
{
    k *= 20.0 / 3.0;
    float h = max( k - abs(a - b), 0.0 ) / k;
    return min(a, b) - h * h * h * h * (5.0 - 2.0 * h) * k * (1.0 / 20.0);
}

// 2. C3/C4 Quintic Polynomial Smooth Maximum (Intersection) & Subtraction
float smaxQuintic( float a, float b, float k )
{
    return -sminQuintic(-a, -b, k);
}

float opSmoothSubtractionQuintic( float a, float b, float k )
{
    return -sminQuintic(a, -b, k);
}

// 3. C3/C4 Quintic Polynomial Smooth Minimum with Analytical 3D Gradient: vec4(d, gx, gy, gz)
vec4 sdgSmoothQuintic( in vec4 a, in vec4 b, in float k )
{
    k *= 20.0 / 3.0;
    float h = max( k - abs(a.x - b.x), 0.0 ) / k;
    float m = 0.5 * h * h * h * (2.0 - h);
    float d = min(a.x, b.x) - h * h * h * h * (5.0 - 2.0 * h) * k * (1.0 / 20.0);
    return vec4( d, mix(a.yzw, b.yzw, (a.x < b.x) ? m : 1.0 - m) );
}

#endif // SDF_GEMS_SMIN_QUINTIC_GLSL
