// ============================================================================
// sdCross3D — Exact 1-Sqrt 3D Orthogonal Cross / Jack via Weyl-Chamber Folding
// @author Discovered using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1` in both exterior and interior; 1 `sqrt` down from 3)
// ============================================================================

// b.x = arm half-length, b.y = arm half-thickness (b.x >= b.y)
float sdCross3D( vec3 p, vec2 b )
{
    p = abs(p);
    if( p.x < p.y ) p.xy = p.yx;
    if( p.y < p.z ) p.yz = p.zy;
    if( p.x < p.y ) p.xy = p.yx;

    vec3 w = p - b.xyy;
    bool e = max(w.x, w.y) > 0.0;
    float l = length(e ? max(w, 0.0) : vec3(max(b.y - p.x, 0.0), -w.y, 0.0));
    return e ? l : -min(-w.x, l);
}

// Infinite 3D orthogonal cross (union of 3 axis-aligned square rods of half-width b)
float sdInfCross3D( vec3 p, float b )
{
    p = abs(p);
    if( p.x < p.y ) p.xy = p.yx;
    if( p.y < p.z ) p.yz = p.zy;
    if( p.x < p.y ) p.xy = p.yx;

    return (p.y > b) ?  length(vec2(p.y - b, max(p.z - b, 0.0)))
                     : -length(vec2(max(b - p.x, 0.0), b - p.y));
}
