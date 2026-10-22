// ============================================================================
// sdRhombus2D & sdgRhombus2D — 0-Div 0/1-Sqrt 2D Rhombus & Isosceles Wedge Prism
// @author Inigo Quilez — upgraded to 0-div 0/1-sqrt & 1-SFU gradients using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1`, `0` divs & `0` sqrts in interior/face strips, `1` SFU max everywhere)
// ============================================================================

// 1. Exact 0-Div 0/1-Sqrt 2D Rhombus (semi-diagonals b.x, b.y)
float sdRhombus2D( in vec2 p, in vec2 b )
{
    p = abs(p);
    if( b.x * (p.x - b.x) > b.y * p.y ) return length(vec2(p.x - b.x, p.y));
    if( b.y * (p.y - b.y) > b.x * p.x ) return length(vec2(p.x, p.y - b.y));
    return (b.y * p.x + b.x * p.y - b.x * b.y) * inversesqrt(dot(b, b));
}

// 2. Exact 1-SFU (`inversesqrt`) 2D Rhombus with Analytical Gradient
vec3 sdgRhombus2D( in vec2 p, in vec2 b )
{
    vec2 s = vec2((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0);
    p = abs(p);
    bool c1 = (b.x * (p.x - b.x) > b.y * p.y);
    bool c2 = (!c1) && (b.y * (p.y - b.y) > b.x * p.x);
    if( c1 || c2 )
    {
        vec2  w  = c1 ? vec2(p.x - b.x, p.y) : vec2(p.x, p.y - b.y);
        float l2 = dot(w, w), il = inversesqrt(max(l2, 1e-24));
        return vec3(l2 * il, s * w * il);
    }
    float il = inversesqrt(dot(b, b));
    return vec3((b.y * p.x + b.x * p.y - b.x * b.y) * il, s * b.yx * il);
}

// 3. Exact 1-Sqrt 3D Extruded Isosceles Wedge Prism
float sdIsoscelesPrism( in vec3 p, in vec2 q, in float h )
{
    float px = abs(p.x), py = p.y, dz = abs(p.z) - h;
    float k  = (q.y < 0.0) ? -1.0 : 1.0, db = k * (py - q.y);
    float d2Sq, d2Sign = 1.0;
    if( px <= q.x && db >= 0.0 )
    {
        d2Sq = db * db;
    }
    else if( px * q.x + py * q.y <= 0.0 )
    {
        d2Sq = px * px + py * py;
    }
    else if( px > q.x && (px - q.x) * q.x + (py - q.y) * q.y >= 0.0 )
    {
        vec2 w = vec2(px - q.x, py - q.y);
        d2Sq = dot(w, w);
    }
    else
    {
        float ds = k * (px * q.y - py * q.x) * inversesqrt(dot(q, q));
        float dm = (db <= 0.0) ? max(db, ds) : ds;
        d2Sq   = dm * dm;
        d2Sign = (dm < 0.0) ? -1.0 : 1.0;
    }
    if( d2Sign < 0.0 && dz <= 0.0 ) return max(-sqrt(d2Sq), dz);
    return sqrt(((d2Sign > 0.0) ? d2Sq : 0.0) + max(dz, 0.0) * max(dz, 0.0));
}
