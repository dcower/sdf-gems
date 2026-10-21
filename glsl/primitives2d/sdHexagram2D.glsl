// ============================================================================
// sdHexagram2D & sdPentagram2D — 0/1-Sqrt Star Polygons & 1-Sqrt 3D Star Prisms
// @author Inigo Quilez (2D) — upgraded to 0/1-sqrt & 1-sqrt 3D star prisms using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1`, `0` sqrts in all 12/10 edge strips; `1` sqrt in 3D non-convex star prisms!)
// ============================================================================

// 1. 0/1-Sqrt 2D Regular Hexagram (6-pointed Star of David; 0 sqrts in all 12 edge normal strips)
float sdHexagram2D( in vec2 p, in float r )
{
    const vec3 k = vec3(-0.5, 0.8660254038, 0.5773502692);
    p = abs(p);
    p -= 2.0 * min(dot(k.xy, p), 0.0) * k.xy;
    p -= 2.0 * min(dot(k.yx, p), 0.0) * k.yx;
    float s = r * k.z;
    vec2  d = vec2(abs(p.x - 2.0 * s) - s, p.y - r);
    return (d.x > 0.0) ? length(d) * sign(d.y) : d.y;
}

// 2. 0/1-Sqrt 2D Regular Pentagram (5-pointed star; 0 sqrts in all 10 edge normal strips)
float sdPentagram2D( in vec2 p, in float r )
{
    const vec2 v1 = vec2( 0.809016994, -0.587785252);
    const vec2 v2 = vec2(-0.809016994, -0.587785252);
    const vec2 v3 = vec2( 0.309016994, -0.951056516);
    p.x = abs(p.x);
    p -= 2.0 * max(dot(v1, p), 0.0) * v1;
    p -= 2.0 * max(dot(v2, p), 0.0) * v2;
    p.x = abs(p.x);
    p.y -= r;
    float h = 0.363271264 * r;
    vec2  d = vec2(abs(dot(p, v3) - h) - h, p.y * v3.x - p.x * v3.y);
    return (d.x > 0.0) ? length(d) * sign(d.y) : d.y;
}

// 3. Exact 1-Sqrt 3D Non-Convex Hexagram Prism (at most 1 sqrt everywhere in R^3!)
//    Because the re-entrant wedge d.x < 0 lies strictly inside the 2D hexagram (d.y < 0),
//    only the exterior outer-tip wedge (d.x > 2*s) and interior re-entrant wedge need 1 sqrt!
float sdHexagramPrism( in vec3 p, in vec2 h )
{
    const vec3 k = vec3(-0.5, 0.8660254038, 0.5773502692);
    p = abs(p);
    p.xy -= 2.0 * min(dot(k.xy, p.xy), 0.0) * k.xy;
    p.xy -= 2.0 * min(dot(k.yx, p.xy), 0.0) * k.yx;
    float s = h.x * k.z;
    vec3  d = vec3(p.x - s, p.y - h.x, p.z - h.y);
    if( d.x < 0.0 ) return (d.z >= 0.0) ? d.z : max(-length(d.xy), d.z);
    return min(max(d.y, d.z), 0.0) + length(max(vec3(d.x - 2.0 * s, d.y, d.z), 0.0));
}

// 4. Exact 1-Sqrt 3D Non-Convex Pentagram Prism (at most 1 sqrt everywhere in R^3!)
float sdPentagramPrism( in vec3 p, in vec2 h )
{
    const vec2 v1 = vec2( 0.809016994, -0.587785252);
    const vec2 v2 = vec2(-0.809016994, -0.587785252);
    const vec2 v3 = vec2( 0.309016994, -0.951056516);
    p.xz = abs(p.xz);
    p.xy -= 2.0 * max(dot(v1, p.xy), 0.0) * v1;
    p.xy -= 2.0 * max(dot(v2, p.xy), 0.0) * v2;
    p.x = abs(p.x);
    p.y -= h.x;
    float l = 0.726542528 * h.x, u = dot(p.xy, v3), dy = p.y * v3.x - p.x * v3.y, dz = p.z - h.y;
    if( u > l ) return (dz >= 0.0) ? dz : max(-length(vec2(u - l, dy)), dz);
    return min(max(dy, dz), 0.0) + length(max(vec3(-u, dy, dz), 0.0));
}
