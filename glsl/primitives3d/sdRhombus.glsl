// ============================================================================
// sdRhombus — Exact 3D Extruded Rounded Rhombus (x=0 Sign Bug Fixed)
// @author Inigo Quilez — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

float sdRhombus( vec3 p, float la, float lb, float h, float ra )
{
    p = abs(p);
    float f = clamp( (la*p.x - lb*p.z + lb*lb) / (la*la + lb*lb), 0.0, 1.0 );
    vec2 w = p.xz - vec2(la, lb) * vec2(f, 1.0 - f);
    vec2 q = vec2( length(w)*sign(w.x + w.y) - ra, p.y - h );
    return min(max(q.x, q.y), 0.0) + length(max(q, 0.0));
}
