// ============================================================================
// sdSquareAntiprism, sdPentAntiprism, sdHexAntiprism, sdAntiprism & sdg* —
// Exact 0-Trig 0/1-Sqrt 3D Antiprisms
//
// @author  Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^3; 0 trig & 0 sqrts in Square,
//          Pentagonal, and Hexagonal Antiprism interior/face prisms, 1 sqrt in edge/vertex cone)
// ============================================================================

#ifndef SDF_GEMS_SD_SQUARE_ANTIPRISM_GLSL
#define SDF_GEMS_SD_SQUARE_ANTIPRISM_GLSL

// 1. Exact 0-Trig 0/1-Sqrt Square Antiprism (base circumradius r, half-height h)
float sdSquareAntiprism( vec3 p, float r, float h )
{
    float x = abs(p.x), y = abs(p.y), z = p.z;
    if( y > x ) { float t = x; x = y; y = t; }

    const float k = 0.70710678, k1 = 0.29289322;
    float u = (x + y) * k, v = (x - y) * k;

    float ub = (z < 0.0) ? x : u;
    float vb = (z < 0.0) ? y : v;
    float zb = abs(z);

    bool  swap_s = (h * (u - x) + r * k1 * z < 0.0);
    float xs = swap_s ? u : x, ys = swap_s ? v : y;
    float us = swap_s ? x : u, zs = swap_s ? -z : z;

    float ls2    = 4.0 * h * h + r * r * k1 * k1;
    float ms     = 2.0 * h * us + r * k1 * zs - r * h * (1.0 + k);
    float d_side = ms * inversesqrt(ls2);
    float d_base = zb - h;
    if( max(d_base, d_side) <= 0.0 || (zb > h && ub <= r * k) ||
        (ms > 0.0 && r * k1 * ms >= ls2 * (zs - h) && k1 * (2.0 * h * h + r * r * k1) * ms <= ls2 * (h * (us - xs) + r * k1 * zs)) )
        return max(d_base, d_side);

    vec3  db = vec3(ub - r * k, max(vb - r * k, 0.0), zb - h);
    vec3  qs = vec3(xs - r, ys, zs - h);
    vec3  es = vec3(-r * k1, r * k, -2.0 * h);
    float ts = clamp(dot(qs, es) / (2.0 * r * r * k1 + 4.0 * h * h), 0.0, 1.0);
    vec3  ds = qs - ts * es;
    return sqrt(min(dot(db, db), dot(ds, ds)));
}

// 2. General N-gonal Antiprism (n >= 3, base circumradius r, half-height h)
float sdAntiprism( vec3 p, float n, float r, float h )
{
    float sp  = 3.14159265 / n;
    float c   = cos(sp), s = sin(sp), c1 = 1.0 - c;
    float rad = length(p.xy);
    float ang = mod(atan(p.y, p.x), 2.0 * sp);
    ang = sp - abs(ang - sp);
    float x = rad * cos(ang), y = rad * sin(ang), z = p.z;
    float u = x * c + y * s,  v = x * s - y * c;

    float ub = (z < 0.0) ? x : u;
    float vb = (z < 0.0) ? y : v;
    float zb = abs(z);

    bool  swap_s = (h * (u - x) + r * c1 * z < 0.0);
    float xs = swap_s ? u : x, ys = swap_s ? v : y;
    float us = swap_s ? x : u, zs = swap_s ? -z : z;

    float ls2    = 4.0 * h * h + r * r * c1 * c1;
    float ms     = 2.0 * h * us + r * c1 * zs - r * h * (1.0 + c);
    float d_side = ms * inversesqrt(ls2);
    float d_base = zb - h;
    if( max(d_base, d_side) <= 0.0 || (zb > h && ub <= r * c) ||
        (ms > 0.0 && r * c1 * ms >= ls2 * (zs - h) && c1 * (2.0 * h * h + r * r * c1) * ms <= ls2 * (h * (us - xs) + r * c1 * zs)) )
        return max(d_base, d_side);

    vec3  db = vec3(ub - r * c, max(vb - r * s, 0.0), zb - h);
    vec3  qs = vec3(xs - r, ys, zs - h);
    vec3  es = vec3(-r * c1, r * s, -2.0 * h);
    float ts = clamp(dot(qs, es) / (2.0 * r * r * c1 + 4.0 * h * h), 0.0, 1.0);
    vec3  ds = qs - ts * es;
    return sqrt(min(dot(db, db), dot(ds, ds)));
}

// 3. Analytical 3D Distance and Unit Gradient vec4(d, grad.xyz) for Square Antiprism
vec4 sdgSquareAntiprism( vec3 p, float r, float h )
{
    vec2  sg = vec2((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0);
    float x = abs(p.x), y = abs(p.y), z = p.z;
    bool  swap_xy = (y > x);
    if( swap_xy ) { float t = x; x = y; y = t; }

    const float k = 0.70710678, k1 = 0.29289322;
    float u = (x + y) * k, v = (x - y) * k;

    bool  neg_z = (z < 0.0);
    float ub = neg_z ? x : u;
    float vb = neg_z ? y : v;
    float zb = abs(z);

    bool  swap_s = (h * (u - x) + r * k1 * z < 0.0);
    float xs = swap_s ? u : x, ys = swap_s ? v : y;
    float us = swap_s ? x : u, zs = swap_s ? -z : z;

    float ls2    = 4.0 * h * h + r * r * k1 * k1;
    float invLs  = inversesqrt(ls2);
    float ms     = 2.0 * h * us + r * k1 * zs - r * h * (1.0 + k);
    float d_side = ms * invLs;
    float d_base = zb - h;
    float d_face = max(d_base, d_side);

    bool in_int  = (d_face <= 0.0);
    bool in_base = (zb > h && ub <= r * k);
    bool in_side = (ms > 0.0 && r * k1 * ms >= ls2 * (zs - h) &&
                    k1 * (2.0 * h * h + r * r * k1) * ms <= ls2 * (h * (us - xs) + r * k1 * zs));
    float d;
    vec3  g;
    if( in_int || in_base || in_side )
    {
        d = d_face;
        g = (in_base || (in_int && d_base >= d_side))
          ? vec3(0.0, 0.0, neg_z ? -1.0 : 1.0)
          : (swap_s ? vec3(2.0 * h * invLs, 0.0, -r * k1 * invLs)
                    : vec3(2.0 * h * k * invLs, 2.0 * h * k * invLs, r * k1 * invLs));
    }
    else
    {
        vec3  db = vec3(ub - r * k, max(vb - r * k, 0.0), zb - h);
        vec3  qs = vec3(xs - r, ys, zs - h);
        vec3  es = vec3(-r * k1, r * k, -2.0 * h);
        float ts = clamp(dot(qs, es) / (2.0 * r * r * k1 + 4.0 * h * h), 0.0, 1.0);
        vec3  ds = qs - ts * es;

        float b2   = dot(db, db), s2 = dot(ds, ds);
        vec3  vb_w = neg_z ? vec3(db.x, db.y, -db.z) : vec3((db.x + db.y) * k, (db.x - db.y) * k, db.z);
        vec3  vs_w = swap_s ? vec3((ds.x + ds.y) * k, (ds.x - ds.y) * k, -ds.z) : ds;
        vec3  ve   = (b2 <= s2) ? vb_w : vs_w;
        float v2m  = min(b2, s2);
        float invV = inversesqrt(max(v2m, 1e-30));
        d = v2m * invV;
        g = ve * invV;
    }

    if( swap_xy ) g.xy = g.yx;
    g.xy *= sg;
    return vec4(d, g);
}

// 4. Exact 0-Trig 0/1-Sqrt Pentagonal Antiprism (base circumradius r, half-height h; uniform at h = 0.5 * r)
float sdPentAntiprism( vec3 p, float r, float h )
{
    const float c = 0.80901699, s = 0.58778525, c1 = 0.19098301;
    float x = p.x, y = abs(p.y), z = p.z;
    float d1 = min(s * x - c * y, 0.0); x -= 2.0 * d1 * s; y += 2.0 * d1 * c;
    float d2 = min(s * x + c * y, 0.0); x -= 2.0 * d2 * s; y -= 2.0 * d2 * c;
    y = abs(y);

    float u  = x * c + y * s, v = x * s - y * c;
    float ub = (z < 0.0) ? x : u, vb = (z < 0.0) ? y : v, zb = abs(z);

    bool  swap_s = (h * (u - x) + r * c1 * z < 0.0);
    float xs = swap_s ? u : x, ys = swap_s ? v : y;
    float us = swap_s ? x : u, zs = swap_s ? -z : z;

    float ls2    = 4.0 * h * h + r * r * c1 * c1;
    float ms     = 2.0 * h * us + r * c1 * zs - r * h * (1.0 + c);
    float d_side = ms * inversesqrt(ls2);
    float d_base = zb - h;
    if( max(d_base, d_side) <= 0.0 || (zb > h && ub <= r * c) ||
        (ms > 0.0 && r * c1 * ms >= ls2 * (zs - h) && c1 * (2.0 * h * h + r * r * c1) * ms <= ls2 * (h * (us - xs) + r * c1 * zs)) )
        return max(d_base, d_side);

    vec3  db = vec3(ub - r * c, max(vb - r * s, 0.0), zb - h);
    vec3  qs = vec3(xs - r, ys, zs - h);
    vec3  es = vec3(-r * c1, r * s, -2.0 * h);
    float ts = clamp(dot(qs, es) / (2.0 * r * r * c1 + 4.0 * h * h), 0.0, 1.0);
    vec3  ds = qs - ts * es;
    return sqrt(min(dot(db, db), dot(ds, ds)));
}

// 5. Analytical 3D Distance and Unit Gradient vec4(d, grad.xyz) for Pentagonal Antiprism
vec4 sdgPentAntiprism( vec3 p, float r, float h )
{
    const float c = 0.80901699, s = 0.58778525, c1 = 0.19098301;
    float sy0 = (p.y < 0.0) ? -1.0 : 1.0;
    float x = p.x, y = abs(p.y), z = p.z;
    float d1 = s * x - c * y; bool r1 = (d1 < 0.0); if( r1 ) { x -= 2.0 * d1 * s; y += 2.0 * d1 * c; }
    float d2 = s * x + c * y; bool r2 = (d2 < 0.0); if( r2 ) { x -= 2.0 * d2 * s; y -= 2.0 * d2 * c; }
    float sy1 = (y < 0.0) ? -1.0 : 1.0;
    y = abs(y);

    float u     = x * c + y * s, v = x * s - y * c;
    bool  neg_z = (z < 0.0);
    float ub    = neg_z ? x : u, vb = neg_z ? y : v, zb = abs(z);

    bool  swap_s = (h * (u - x) + r * c1 * z < 0.0);
    float xs = swap_s ? u : x, ys = swap_s ? v : y;
    float us = swap_s ? x : u, zs = swap_s ? -z : z;

    float ls2    = 4.0 * h * h + r * r * c1 * c1;
    float invLs  = inversesqrt(ls2);
    float ms     = 2.0 * h * us + r * c1 * zs - r * h * (1.0 + c);
    float d_side = ms * invLs;
    float d_base = zb - h;
    float d_face = max(d_base, d_side);

    bool in_int  = (d_face <= 0.0);
    bool in_base = (zb > h && ub <= r * c);
    bool in_side = (ms > 0.0 && r * c1 * ms >= ls2 * (zs - h) &&
                    c1 * (2.0 * h * h + r * r * c1) * ms <= ls2 * (h * (us - xs) + r * c1 * zs));
    float d;
    vec3  g;
    if( in_int || in_base || in_side )
    {
        d = d_face;
        g = (in_base || (in_int && d_base >= d_side))
          ? vec3(0.0, 0.0, neg_z ? -1.0 : 1.0)
          : (swap_s ? vec3(2.0 * h * invLs, 0.0, -r * c1 * invLs)
                    : vec3(2.0 * h * c * invLs, 2.0 * h * s * invLs, r * c1 * invLs));
    }
    else
    {
        vec3  db = vec3(ub - r * c, max(vb - r * s, 0.0), zb - h);
        vec3  qs = vec3(xs - r, ys, zs - h);
        vec3  es = vec3(-r * c1, r * s, -2.0 * h);
        float ts = clamp(dot(qs, es) / (2.0 * r * r * c1 + 4.0 * h * h), 0.0, 1.0);
        vec3  ds = qs - ts * es;

        float b2   = dot(db, db), s2 = dot(ds, ds);
        vec3  vb_w = neg_z ? vec3(db.x, db.y, -db.z) : vec3(db.x * c + db.y * s, db.x * s - db.y * c, db.z);
        vec3  vs_w = swap_s ? vec3(ds.x * c + ds.y * s, ds.x * s - ds.y * c, -ds.z) : ds;
        vec3  ve   = (b2 <= s2) ? vb_w : vs_w;
        float v2m  = min(b2, s2);
        float invV = inversesqrt(max(v2m, 1e-30));
        d = v2m * invV;
        g = ve * invV;
    }

    g.y *= sy1;
    if( r2 ) { float gd = s * g.x + c * g.y; g.x -= 2.0 * gd * s; g.y -= 2.0 * gd * c; }
    if( r1 ) { float gd = s * g.x - c * g.y; g.x -= 2.0 * gd * s; g.y += 2.0 * gd * c; }
    g.y *= sy0;
    return vec4(d, g);
}

// 6. Exact 0-Trig 0/1-Sqrt Hexagonal Antiprism (base circumradius r, half-height h)
float sdHexAntiprism( vec3 p, float r, float h )
{
    const float c = 0.86602540, s = 0.5, c1 = 0.13397460, k = 1.73205081;
    float x = abs(p.x), y = abs(p.y), z = p.z;
    if( x - k * y < 0.0 ) { float t = x; x = 0.5 * (t + k * y); y = 0.5 * abs(k * t - y); }

    float u  = x * c + y * s, v = x * s - y * c;
    float ub = (z < 0.0) ? x : u, vb = (z < 0.0) ? y : v, zb = abs(z);

    bool  swap_s = (h * (u - x) + r * c1 * z < 0.0);
    float xs = swap_s ? u : x, ys = swap_s ? v : y;
    float us = swap_s ? x : u, zs = swap_s ? -z : z;

    float ls2    = 4.0 * h * h + r * r * c1 * c1;
    float ms     = 2.0 * h * us + r * c1 * zs - r * h * (1.0 + c);
    float d_side = ms * inversesqrt(ls2);
    float d_base = zb - h;
    if( max(d_base, d_side) <= 0.0 || (zb > h && ub <= r * c) ||
        (ms > 0.0 && r * c1 * ms >= ls2 * (zs - h) && c1 * (2.0 * h * h + r * r * c1) * ms <= ls2 * (h * (us - xs) + r * c1 * zs)) )
        return max(d_base, d_side);

    vec3  db = vec3(ub - r * c, max(vb - r * s, 0.0), zb - h);
    vec3  qs = vec3(xs - r, ys, zs - h);
    vec3  es = vec3(-r * c1, r * s, -2.0 * h);
    float ts = clamp(dot(qs, es) / (2.0 * r * r * c1 + 4.0 * h * h), 0.0, 1.0);
    vec3  ds = qs - ts * es;
    return sqrt(min(dot(db, db), dot(ds, ds)));
}

// 7. Analytical 3D Distance and Unit Gradient vec4(d, grad.xyz) for Hexagonal Antiprism
vec4 sdgHexAntiprism( vec3 p, float r, float h )
{
    const float c = 0.86602540, s = 0.5, c1 = 0.13397460, k = 1.73205081;
    vec2  sg = vec2((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0);
    float x = abs(p.x), y = abs(p.y), z = p.z;
    bool  rf = (x - k * y < 0.0);
    if( rf ) { float t = x; x = 0.5 * (t + k * y); y = 0.5 * (k * t - y); }
    float sy1 = (y < 0.0) ? -1.0 : 1.0;
    y = abs(y);

    float u     = x * c + y * s, v = x * s - y * c;
    bool  neg_z = (z < 0.0);
    float ub    = neg_z ? x : u, vb = neg_z ? y : v, zb = abs(z);

    bool  swap_s = (h * (u - x) + r * c1 * z < 0.0);
    float xs = swap_s ? u : x, ys = swap_s ? v : y;
    float us = swap_s ? x : u, zs = swap_s ? -z : z;

    float ls2    = 4.0 * h * h + r * r * c1 * c1;
    float invLs  = inversesqrt(ls2);
    float ms     = 2.0 * h * us + r * c1 * zs - r * h * (1.0 + c);
    float d_side = ms * invLs;
    float d_base = zb - h;
    float d_face = max(d_base, d_side);

    bool in_int  = (d_face <= 0.0);
    bool in_base = (zb > h && ub <= r * c);
    bool in_side = (ms > 0.0 && r * c1 * ms >= ls2 * (zs - h) &&
                    c1 * (2.0 * h * h + r * r * c1) * ms <= ls2 * (h * (us - xs) + r * c1 * zs));
    float d;
    vec3  g;
    if( in_int || in_base || in_side )
    {
        d = d_face;
        g = (in_base || (in_int && d_base >= d_side))
          ? vec3(0.0, 0.0, neg_z ? -1.0 : 1.0)
          : (swap_s ? vec3(2.0 * h * invLs, 0.0, -r * c1 * invLs)
                    : vec3(2.0 * h * c * invLs, 2.0 * h * s * invLs, r * c1 * invLs));
    }
    else
    {
        vec3  db = vec3(ub - r * c, max(vb - r * s, 0.0), zb - h);
        vec3  qs = vec3(xs - r, ys, zs - h);
        vec3  es = vec3(-r * c1, r * s, -2.0 * h);
        float ts = clamp(dot(qs, es) / (2.0 * r * r * c1 + 4.0 * h * h), 0.0, 1.0);
        vec3  ds = qs - ts * es;

        float b2   = dot(db, db), s2 = dot(ds, ds);
        vec3  vb_w = neg_z ? vec3(db.x, db.y, -db.z) : vec3(db.x * c + db.y * s, db.x * s - db.y * c, db.z);
        vec3  vs_w = swap_s ? vec3(ds.x * c + ds.y * s, ds.x * s - ds.y * c, -ds.z) : ds;
        vec3  ve   = (b2 <= s2) ? vb_w : vs_w;
        float v2m  = min(b2, s2);
        float invV = inversesqrt(max(v2m, 1e-30));
        d = v2m * invV;
        g = ve * invV;
    }

    g.y *= sy1;
    if( rf ) g.xy = vec2(g.x + k * g.y, k * g.x - g.y) * 0.5;
    g.xy *= sg;
    return vec4(d, g);
}

#endif // SDF_GEMS_SD_SQUARE_ANTIPRISM_GLSL
