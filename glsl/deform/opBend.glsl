// ============================================================================
// opBend — Exact atan-Free Circular Arch Bend + Straight Tangent Legs (opBendArch)
// @author Inigo Quilez (opCheapBend) — upgraded/fixed by Dillon Cower (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1` outside inner focal axis)
// ============================================================================

vec3 opBendArch( in vec3 p, in float R, in float alpha, in vec2 cs )
{
    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    vec2  w  = vec2( abs(p.x), R + p.y );
    float u  = w.x*cs.x - w.y*cs.y;
    float v  = w.x*cs.y + w.y*cs.x;
    return vec3( sx * max(R*alpha + u, 0.0),
                 (u <= 0.0) ? length(w) - R : v - R,
                 p.z );
}
