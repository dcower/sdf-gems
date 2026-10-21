// ============================================================================
// sdSphericalZone — Exact 3D Hollow Spherical Zone (Two-Plane Cut Shell)
// @author Discovered using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

float sdSphericalZone( vec3 p, float r, float h1, float h2, float t )
{
    float w1 = sqrt(r*r - h1*h1);
    float w2 = sqrt(r*r - h2*h2);
    vec2  q  = vec2( length(p.xz), p.y );
    float d  = (h1*q.x > w1*q.y) ? length(q - vec2(w1, h1)) :
               (w2*q.y > h2*q.x) ? length(q - vec2(w2, h2)) :
                                   abs(length(q) - r);
    return d - t;
}
