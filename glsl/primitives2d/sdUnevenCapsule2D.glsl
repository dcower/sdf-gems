// ============================================================================
// sdUnevenCapsule2D — Exact 1-Sqrt 2D Uneven Capsule (Signed-Square Tangent-Cone Reduction)
// @author Inigo Quilez — upgraded/optimized using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1`, 1 `sqrt` everywhere down from 2 `sqrt`s, `0.0` error)
// ============================================================================

float sdUnevenCapsule2D( vec2 p, float r1, float r2, float h )
{
    p.x = abs(p.x);
    float b = (r1 - r2)/h, a2 = 1.0 - b*b;
    float y = p.y, z = y - h, s = b*abs(b)*p.x*p.x;
    return (a2*y*abs(y) < s) ? length(p) - r1 :
           (a2*z*abs(z) > s) ? length(vec2(p.x, z)) - r2 :
                               sqrt(a2)*p.x + b*y - r1;
}
