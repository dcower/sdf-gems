// ============================================================================
// sdCappedCone — Exact 1-to-2 Sqrt 3D Capped Cone via Endpoint Subsumption
// @author Inigo Quilez — upgraded/optimized using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, 1 `sqrt` in flat-cap Voronoi zone & 2 in mantle, down from 2 everywhere)
// @requires dot2 (from common.glsl)
// ============================================================================

float sdCappedCone( vec3 p, float h, float r1, float r2 )
{
    vec2  q  = vec2(length(p.xz), p.y);
    vec2  k1 = vec2(r2, h);
    vec2  k2 = vec2(r2 - r1, 2.0 * h);
    float rc = (q.y < 0.0) ? r1 : r2;
    float cy = abs(q.y) - h;
    vec2  cb = q - k1 + k2 * clamp(dot(k1 - q, k2) / dot2(k2), 0.0, 1.0);
    float b2 = dot2(cb);
    if( q.x <= rc && cy * cy <= b2 ) return cy;
    return ((cb.x < 0.0 && cy < 0.0) ? -1.0 : 1.0) * sqrt(b2);
}
