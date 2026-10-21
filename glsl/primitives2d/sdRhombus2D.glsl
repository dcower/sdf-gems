// ============================================================================
// sdRhombus2D — 0-Div 0/1-Sqrt 2D Rhombus & Isosceles Triangle + Gradients
// @author Inigo Quilez — upgraded to 0-div 0/1-sqrt & 1-SFU gradients using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1`, `0` divs & `0` sqrts in interior/face strips, `1` SFU max everywhere)
// ============================================================================
//
// Mathematical Discovery (Direct Vertex-Wedge & 0-Division Edge Decomposition):
//   1. In IQ's sdRhombus(p, b) (Shadertoy XdXcRB), every pixel evaluates a
//      runtime division `/ dot(b, b)` followed by `length(p - 0.5*b*vec2(1-h, 1+h))`
//      (2 SFU ops: RCP + SQRT everywhere, even inside the rhombus).
//      In the first quadrant p = abs(p), a point projects past the right vertex
//      (b.x, 0) iff `b.x * (p.x - b.x) > b.y * p.y`, and past the top vertex
//      (0, b.y) iff `b.y * (p.y - b.y) > b.x * p.x`. In those two vertex wedges,
//      the Euclidean distance is simply `length(p - vertex)` (1 sqrt, 0 divs,
//      independent of |b|!). Everywhere else (100% of the interior and all 4
//      exterior face-normal strips), the distance is linear:
//        d = (b.y * p.x + b.x * p.y - b.x * b.y) * inversesqrt(dot(b, b))
//      requiring 0 sqrts for constant b (1 RSQ for dynamic b) and 0 divisions!
//   2. In IQ's sdTriangleIsosceles(p, q) (Shadertoy MldcD7), every pixel
//      evaluates 2 runtime divisions (`/ dot(q, q)` and `/ q.x`) plus `sqrt(d.x)`
//      (3 SFU ops total). By testing the base strip (`p.x <= q.x && db >= 0`),
//      the apex wedge (`dot(p, q) <= 0`), and the base-corner wedge
//      (`p.x > q.x && dot(p - q, q) >= 0`), all runtime divisions are eliminated
//      and the interior + all 3 exterior face strips require 0 sqrts for fixed q!

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
    if( c1 || c2 ) {
        vec2  w  = c1 ? vec2(p.x - b.x, p.y) : vec2(p.x, p.y - b.y);
        float l2 = dot(w, w), il = inversesqrt(max(l2, 1e-24));
        return vec3(l2 * il, s * w * il);
    }
    float il = inversesqrt(dot(b, b));
    return vec3((b.y * p.x + b.x * p.y - b.x * b.y) * il, s * b.yx * il);
}

// 3. Exact 0-Div 0/1-Sqrt 2D Isosceles Triangle (apex at origin, base half-width q.x at y = q.y)
float sdTriangleIsosceles2D( in vec2 p, in vec2 q )
{
    p.x = abs(p.x);
    float k = (q.y < 0.0) ? -1.0 : 1.0, db = k * (p.y - q.y);
    if( p.x <= q.x && db >= 0.0 ) return db;
    if( dot(p, q) <= 0.0 ) return length(p);
    if( p.x > q.x && dot(p - q, q) >= 0.0 ) return length(p - q);
    float ds = k * (p.x * q.y - p.y * q.x) * inversesqrt(dot(q, q));
    return (db <= 0.0) ? max(db, ds) : ds;
}

// 4. Exact 1-SFU (`inversesqrt`) 2D Isosceles Triangle with Analytical Gradient
vec3 sdgTriangleIsosceles2D( in vec2 p, in vec2 q )
{
    float sx = (p.x < 0.0) ? -1.0 : 1.0, k = (q.y < 0.0) ? -1.0 : 1.0;
    p.x = abs(p.x);
    float db = k * (p.y - q.y);
    if( p.x <= q.x && db >= 0.0 ) return vec3(db, 0.0, k);
    bool c_apex = (dot(p, q) <= 0.0);
    bool c_corn = (!c_apex) && (p.x > q.x) && (dot(p - q, q) >= 0.0);
    if( c_apex || c_corn ) {
        vec2  w  = c_apex ? p : (p - q);
        float l2 = dot(w, w), il = inversesqrt(max(l2, 1e-24));
        return vec3(l2 * il, vec2(sx * w.x, w.y) * il);
    }
    float il = inversesqrt(dot(q, q));
    float ds = k * (p.x * q.y - p.y * q.x) * il;
    if( db <= 0.0 && db > ds ) return vec3(db, 0.0, k);
    return vec3(ds, vec2(sx * k * q.y, -k * q.x) * il);
}

// 5. Exact 1-Sqrt 3D Extruded Isosceles Wedge Prism
float sdIsoscelesPrism( in vec3 p, in vec2 q, in float h )
{
    float px = abs(p.x), py = p.y, dz = abs(p.z) - h;
    float k  = (q.y < 0.0) ? -1.0 : 1.0, db = k * (py - q.y);
    float d2Sq, d2Sign = 1.0;
    if( px <= q.x && db >= 0.0 ) {
        d2Sq = db * db;
    } else if( px * q.x + py * q.y <= 0.0 ) {
        d2Sq = px * px + py * py;
    } else if( px > q.x && (px - q.x) * q.x + (py - q.y) * q.y >= 0.0 ) {
        vec2 w = vec2(px - q.x, py - q.y);
        d2Sq = dot(w, w);
    } else {
        float ds = k * (px * q.y - py * q.x) * inversesqrt(dot(q, q));
        float dm = (db <= 0.0) ? max(db, ds) : ds;
        d2Sq   = dm * dm;
        d2Sign = (dm < 0.0) ? -1.0 : 1.0;
    }
    if( d2Sign < 0.0 && dz <= 0.0 ) return max(-sqrt(d2Sq), dz);
    return sqrt(((d2Sign > 0.0) ? d2Sq : 0.0) + max(dz, 0.0) * max(dz, 0.0));
}
