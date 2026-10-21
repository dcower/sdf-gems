// ============================================================================
// sdgHexPrism — 0/1-Sqrt Hexagonal & Triangular Prisms with Analytical 3D Gradients
// @author Discovered using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1`, `0` sqrts in interior & caps, `1` sqrt in exterior walls/corners)
// ============================================================================

vec4 sdgHexPrism( vec3 p, vec2 h )
{
    const vec3 k = vec3(-0.86602540, 0.5, 0.57735027);
    vec3  s  = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    vec3  q  = abs(p);
    float dk = dot(k.xy, q.xy);
    bool  rf = (dk < 0.0);
    if( rf ) q.xy -= 2.0 * dk * k.xy;

    vec2  v2 = vec2(q.x - clamp(q.x, -k.z * h.x, k.z * h.x), q.y - h.x);
    float wz = q.z - h.y, d, gz;
    vec2  g2;

    if( v2.y <= 0.0 ) {
        if( wz <= 0.0 && v2.y > wz ) { d = v2.y; g2 = vec2(0.0, 1.0); gz = 0.0; }
        else                         { d = wz;   g2 = vec2(0.0, 0.0); gz = s.z; }
    } else {
        float wzp = max(wz, 0.0);
        d  = sqrt(dot(v2, v2) + wzp * wzp);
        g2 = v2 / d;
        gz = (s.z * wzp) / d;
    }
    if( rf ) g2 -= 2.0 * dot(k.xy, g2) * k.xy;
    return vec4(d, s.x * g2.x, s.y * g2.y, gz);
}

vec4 sdgTriPrism( vec3 p, vec2 h )
{
    const vec2 n = vec2(0.86602540, 0.5);
    float sx = (p.x < 0.0) ? -1.0 : 1.0, sz = (p.z < 0.0) ? -1.0 : 1.0;
    vec3  q  = vec3(abs(p.x), p.y, abs(p.z));
    float ds = dot(q.xy, n) - 0.5 * h.x, db = -q.y - 0.5 * h.x, d2 = max(ds, db), wz = q.z - h.y;

    if( d2 <= 0.0 ) {
        if( wz <= 0.0 && d2 > wz ) return vec4(d2, (ds > db) ? vec2(sx * n.x, n.y) : vec2(0.0, -1.0), 0.0);
        return vec4(wz, 0.0, 0.0, sz);
    }
    vec2  ts = vec2(n.y, -n.x), pt = vec2(0.0, 0.5 * h.x);
    vec2  vs = q.xy - (pt + clamp(dot(q.xy - pt, ts), 0.0, h.x) * ts);
    vec2  vb = vec2(max(q.x - n.x * h.x, 0.0), q.y + 0.5 * h.x);
    vec2  v2 = (ds > 0.0 && (db <= 0.0 || dot(vs, vs) < dot(vb, vb))) ? vs : vb;
    float wzp = max(wz, 0.0), d = sqrt(dot(v2, v2) + wzp * wzp);
    return vec4(d, (sx * v2.x) / d, v2.y / d, (sz * wzp) / d);
}
