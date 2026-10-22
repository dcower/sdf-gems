// ============================================================================
// sdEquilateralTriangle2D & sdgEquilateralTriangle2D — Exact 0/1-Sqrt 2D Equilateral Triangle & Gradient
// @author Inigo Quilez — upgraded to 0-sqrt interior/face-strip & 0/1-SFU gradient using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1`, `0` sqrts in interior & face strips, `1` in corner wedges)
// ============================================================================

#ifndef SDF_GEMS_SD_EQUILATERAL_TRIANGLE2D_GLSL
#define SDF_GEMS_SD_EQUILATERAL_TRIANGLE2D_GLSL

float sdEquilateralTriangle2D( vec2 p, float r )
{
    const float k = sqrt(3.0);
    p.x = abs(p.x) - r;
    p.y = p.y + r / k;
    if( p.x + k * p.y > 0.0 ) p = vec2(p.x - k * p.y, -k * p.x - p.y) / 2.0;
    float dx = abs(p.x + r) - r;
    return (dx > 0.0) ? length(vec2(dx, p.y)) : -p.y;
}

vec3 sdgEquilateralTriangle2D( vec2 p, float r )
{
    const float k = 1.73205081;
    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    p.x = abs(p.x) - r;
    p.y = p.y + r / k;
    bool rf = (p.x + k * p.y > 0.0);
    if( rf ) p = vec2(p.x - k * p.y, -k * p.x - p.y) * 0.5;

    float ux = p.x + r;
    float dx = abs(ux) - r;
    float d;
    vec2  g;
    if( dx <= 0.0 )
    {
        d = -p.y;
        g = vec2(0.0, -1.0);
    }
    else
    {
        vec2  v    = vec2((ux < 0.0) ? -dx : dx, p.y);
        float v2   = dot(v, v);
        float invV = inversesqrt(max(v2, 1e-30));
        d = v2 * invV;
        g = v * invV;
    }
    if( rf ) g = vec2(g.x - k * g.y, -k * g.x - g.y) * 0.5;
    g.x *= sx;
    return vec3(d, g);
}

#endif // SDF_GEMS_SD_EQUILATERAL_TRIANGLE2D_GLSL
