// ============================================================================
// sdSnubDodecahedron & sdPentagonalHexecontahedron (+ sdg*)
// Exact 0/1-Sqrt 3D Chiral Archimedean Snub Dodecahedron (92 Faces, 60 Vertices)
// & Catalan Dual Pentagonal Hexecontahedron (60 Pentagons, 92 Vertices, 150 Edges)
//
// @author  Discovered and derived using Gemini 4 Carbon (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^3; 0 sqrts in interior & all 92/60 face prisms,
//          1 sqrt in the exterior edge/vertex cones)
// @novelty Exact 0-trig, 0/1-sqrt Euclidean SDFs and analytical unit gradients for the two
//          60-element chiral icosahedral solids (symmetry group I subset SO(3), no mirror planes):
//            1. Icosahedral Chiral Parity Fold: folds p in R^3 via the 5-step I_h reflection
//               network while tracking determinant parity; skipping abs(p.x) into a single
//               parity sign flip (if (odd) p.x = -p.x) lands every point in the double chamber
//               D_1 = C_{I_h} u s_x(C_{I_h}), a C^0-continuous fundamental domain of chiral I.
//            2. Exact 6-Plane Interior Collapse (Snub Dodecahedron): LP proof over all 92 face
//               planes shows only 6 planes are active in D_1, pairing via |p.x| and |w| into
//               4 dot expressions (0 sqrts).
//            3. Nearest-Vertex Tangent-Cone & Intrinsic Bilateral Frame: LP proof over all
//               normal cones establishes that in D_1 only 3 of the 60 snub vertices (I, R_3^T, R_2)
//               are ever nearest, and the 3.3.3.3.5 vertex star has exact bilateral symmetry in
//               its orthonormal tangent frame (v_0, b_1, b_2), collapsing 5 face wedges + 5 edge
//               rays to 3 wedges + 3 shared-apex rays for sdSnubDodecahedron, and reducing all
//               60 faces of sdPentagonalHexecontahedron (via the Argmax-Face Theorem) to 1 fixed
//               mirror-symmetric pentagon with 3 edges.
//
// Algebraic constants: phi = (1 + sqrt(5))/2, xi = 0.943151259 (positive root of
//   xi^3 + 2*xi^2 - phi^2 = 0), unit snub generator v_0 = (phi*(1 - xi^2), xi^2 + 2*xi - phi - 1, 1)
//   / sqrt(2 - xi). Mirror enantiomorphs: pass vec3(-p.x, p.y, p.z) (and negate grad.x for sdg*).
//
// Parameters:
//   sdSnubDodecahedron(p, r):
//     r : inradius to the 12 regular pentagonal faces (insphere radius; 80 equilateral triangles
//         at inradius 1.048550123 * r, circumradius 1.088303306 * r, edge length 0.504816977 * r)
//   sdPentagonalHexecontahedron(p, r):
//     r : insphere radius tangent to all 60 congruent irregular pentagonal faces (12 five-valent
//         vertices at 1.088303306 * r, 80 three-valent at 1.037912526 * r; 2 long edges of length
//         0.500025330 * r and 3 equal short edges of length 0.285752834 * r per face)
// ============================================================================

#ifndef SDF_GEMS_SD_SNUB_DODECAHEDRON_GLSL
#define SDF_GEMS_SD_SNUB_DODECAHEDRON_GLSL

float sdSnubDodecahedron( vec3 p, float r )
{
    // 1. Icosahedral chiral parity fold into D_1 = C_{I_h} u s_x(C_{I_h}) (fundamental domain of I).
    const vec3 n3 = vec3(-0.5, -0.809016994, 0.309016994);
    bool odd = false;
    for( int i = 0; i < 3; i++ )
    {
        odd = odd ^^ (p.x < 0.0) ^^ (p.y < 0.0) ^^ (p.z < 0.0);
        p = abs(p);
        float dn = dot(p, n3);
        if( dn < 0.0 ) { p -= 2.0 * dn * n3; odd = !odd; }
    }
    if( p.y < 0.0 ) { p.y = -p.y; odd = !odd; }
    if( odd ) p.x = -p.x;

    // 2. Exact interior (6 active planes in D_1 collapsed to 4 expressions, 0 sqrts).
    float dp  = abs(p.x) * 0.525731112 + p.z * 0.850650808 - r;
    float dt0 = p.y * 0.356822090 + p.z * 0.934172359;
    float dt1 = abs(p.x * 0.091107505 - p.y * 0.103194006) + p.z * 0.990479894;
    float dt2 = dot(p, vec3(0.268143168, 0.389643726, 0.881065837));
    float dm  = max(dp, max(dt0, max(dt1, dt2)) - 1.048550123 * r);
    if( dm <= 0.0 ) return dm;

    // 3. Nearest snub vertex in D_1 among v_0, R_3^T v_0, and R_2 v_0; rotate it to v_0.
    float sb = dot(p, vec3(-0.437202286, 0.144773713, -0.055298638));
    float sc = dot(p.xy, vec2(-0.347727210, -0.306999989));
    if( max(sb, sc) > 0.0 )
    {
        if( sb >= sc ) { p.x = -p.x; p -= 2.0 * dot(p, n3) * n3; }
        else           { p.xy = -p.xy; }
    }

    // 4. Intrinsic bilateral vertex-star frame u = (dot(q, v_0), dot(q, b_1), |dot(q, b_2)|).
    vec3 u = vec3(dot(p, vec3(0.173863605,  0.153499994,  0.972732851)) - 1.088303306 * r,
                  dot(p, vec3(0.927503489, -0.357456413, -0.109371800)),
              abs(dot(p, vec3(0.330921025,  0.921228888, -0.204520441))));

    float hp = dot(u.xy, vec2(0.918861492, 0.394580230));
    if( hp > -1e-5 * r && dot(u, vec3(0.319222112, -0.743374563, 0.587785252)) <= 0.0 )
        return hp;

    float uxm = 0.133903947 * u.x;
    float ha  = dot(u, vec3(0.963472330, 0.013059319, 0.267489295));
    if( ha > -1e-5 * r && uxm + dot(u.yz, vec2( 0.841503864, -0.523393714)) <= 0.0
                       && uxm + dot(u.yz, vec2(-0.888486393, -0.438932412)) <= 0.0 )
        return ha;

    float hb = dot(u, vec3(0.963472330, -0.229720457, 0.137657474));
    if( hb > -1e-5 * r && uxm + dot(u.yz, vec2( 0.858374701,  0.495239947)) <= 0.0
                       && uxm + dot(u.yz, vec2(-0.031926683, -0.990479894)) <= 0.0 )
        return hb;

    // 5. Closest of the 3 folded edge rays (sharing u.x component -0.231928440) or vertex v_0.
    const vec3 e1 = vec3(-0.231928440,  0.540093234, 0.809016994);
    const vec3 e2 = vec3(-0.231928440, -0.458717106, 0.857780750);
    const vec3 e3 = vec3(-0.231928440, -0.972732851, 0.0);
    float tx = -0.231928440 * u.x;
    float t1 = tx + dot(u.yz, e1.yz);
    float t2 = tx + dot(u.yz, e2.yz);
    float t3 = tx + u.y * e3.y;
    float tu = 0.0;
    vec3  wu = vec3(0.0);
    if( t1 > tu ) { tu = t1; wu = e1; }
    if( t2 > tu ) { tu = t2; wu = e2; }
    if( t3 > tu ) { tu = t3; wu = e3; }
    return length(u - tu * wu);
}

vec4 sdgSnubDodecahedron( vec3 p, float r )
{
    const vec3 n3 = vec3(-0.5, -0.809016994, 0.309016994);
    vec3 s0 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    bool odd = (p.x < 0.0) ^^ (p.y < 0.0) ^^ (p.z < 0.0);
    p = abs(p);
    float dn0 = dot(p, n3); bool r0 = (dn0 < 0.0); if( r0 ) { p -= 2.0 * dn0 * n3; odd = !odd; }

    vec3 s1 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    odd = odd ^^ (p.x < 0.0) ^^ (p.y < 0.0) ^^ (p.z < 0.0);
    p = abs(p);
    float dn1 = dot(p, n3); bool r1 = (dn1 < 0.0); if( r1 ) { p -= 2.0 * dn1 * n3; odd = !odd; }

    vec3 s2 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    odd = odd ^^ (p.x < 0.0) ^^ (p.y < 0.0) ^^ (p.z < 0.0);
    p = abs(p);
    float dn2 = dot(p, n3); bool r2 = (dn2 < 0.0); if( r2 ) { p -= 2.0 * dn2 * n3; odd = !odd; }

    bool sy3 = (p.y < 0.0);
    if( sy3 ) { p.y = -p.y; odd = !odd; }
    if( odd ) p.x = -p.x;

    const float K   = 1.048550123;
    float w19 = p.x * 0.091107505 - p.y * 0.103194006;
    float dp  = abs(p.x) * 0.525731112 + p.z * 0.850650808 - r;
    float dt0 = p.y * 0.356822090 + p.z * 0.934172359 - K * r;
    float dt1 = abs(w19) + p.z * 0.990479894 - K * r;
    float dt2 = dot(p, vec3(0.268143168, 0.389643726, 0.881065837)) - K * r;
    float dm  = max(max(dp, dt0), max(dt1, dt2));

    float d;
    vec3  g;
    int   rot = 0;
    if( dm <= 0.0 )
    {
        d = dm;
        g = (dp >= dm)  ? vec3((p.x < 0.0) ? -0.525731112 : 0.525731112, 0.0, 0.850650808)
          : (dt0 >= dm) ? vec3(0.0, 0.356822090, 0.934172359)
          : (dt1 >= dm) ? vec3((w19 < 0.0) ? -0.091107505 : 0.091107505,
                               (w19 < 0.0) ?  0.103194006 : -0.103194006, 0.990479894)
                        : vec3(0.268143168, 0.389643726, 0.881065837);
    }
    else
    {
        float sb = dot(p, vec3(-0.437202286, 0.144773713, -0.055298638));
        float sc = dot(p.xy, vec2(-0.347727210, -0.306999989));
        if( max(sb, sc) > 0.0 )
        {
            if( sb >= sc ) { p.x = -p.x; p -= 2.0 * dot(p, n3) * n3; rot = 1; }
            else           { p.xy = -p.xy;                            rot = 2; }
        }

        const vec3 v0 = vec3(0.173863605,  0.153499994,  0.972732851);
        const vec3 b1 = vec3(0.927503489, -0.357456413, -0.109371800);
        const vec3 b2 = vec3(0.330921025,  0.921228888, -0.204520441);
        float uz_raw = dot(p, b2);
        vec3  u      = vec3(dot(p, v0) - 1.088303306 * r, dot(p, b1), abs(uz_raw));

        float hp  = dot(u.xy, vec2(0.918861492, 0.394580230));
        float uxm = 0.133903947 * u.x;
        float ha  = dot(u, vec3(0.963472330,  0.013059319, 0.267489295));
        float hb  = dot(u, vec3(0.963472330, -0.229720457, 0.137657474));
        vec3  gu;
        if( hp > -1e-5 * r && dot(u, vec3(0.319222112, -0.743374563, 0.587785252)) <= 0.0 )
        {
            d  = hp;
            gu = vec3(0.918861492, 0.394580230, 0.0);
        }
        else if( ha > -1e-5 * r && uxm + dot(u.yz, vec2( 0.841503864, -0.523393714)) <= 0.0
                                && uxm + dot(u.yz, vec2(-0.888486393, -0.438932412)) <= 0.0 )
        {
            d  = ha;
            gu = vec3(0.963472330, 0.013059319, 0.267489295);
        }
        else if( hb > -1e-5 * r && uxm + dot(u.yz, vec2( 0.858374701,  0.495239947)) <= 0.0
                                && uxm + dot(u.yz, vec2(-0.031926683, -0.990479894)) <= 0.0 )
        {
            d  = hb;
            gu = vec3(0.963472330, -0.229720457, 0.137657474);
        }
        else
        {
            const vec3 e1 = vec3(-0.231928440,  0.540093234, 0.809016994);
            const vec3 e2 = vec3(-0.231928440, -0.458717106, 0.857780750);
            const vec3 e3 = vec3(-0.231928440, -0.972732851, 0.0);
            float tx = -0.231928440 * u.x;
            float t1 = tx + dot(u.yz, e1.yz);
            float t2 = tx + dot(u.yz, e2.yz);
            float t3 = tx + u.y * e3.y;
            float tu = 0.0;
            vec3  wu = vec3(0.0);
            if( t1 > tu ) { tu = t1; wu = e1; }
            if( t2 > tu ) { tu = t2; wu = e2; }
            if( t3 > tu ) { tu = t3; wu = e3; }
            vec3  vu   = u - tu * wu;
            float v2   = dot(vu, vu);
            float invV = inversesqrt(max(v2, 1e-30));
            d  = v2 * invV;
            gu = (v2 > 1e-30) ? vu * invV : vec3(1.0, 0.0, 0.0);
        }
        if( uz_raw < 0.0 ) gu.z = -gu.z;
        g = gu.x * v0 + gu.y * b1 + gu.z * b2;
    }

    if( rot == 1 ) { g -= 2.0 * dot(g, n3) * n3; g.x = -g.x; }
    if( rot == 2 ) { g.xy = -g.xy; }
    if( odd ) g.x = -g.x;
    if( sy3 ) g.y = -g.y;
    if( r2 ) g -= 2.0 * dot(g, n3) * n3; g *= s2;
    if( r1 ) g -= 2.0 * dot(g, n3) * n3; g *= s1;
    if( r0 ) g -= 2.0 * dot(g, n3) * n3; g *= s0;
    return vec4(d, g);
}

// ============================================================================
// sdPentagonalHexecontahedron & sdgPentagonalHexecontahedron
// Exact 0/1-Sqrt Chiral Catalan Dual of the Snub Dodecahedron (60 Pentagons)
// ============================================================================

float sdPentagonalHexecontahedron( vec3 p, float r )
{
    // 1. Icosahedral chiral parity fold into D_1 = C_{I_h} u s_x(C_{I_h}).
    const vec3 n3 = vec3(-0.5, -0.809016994, 0.309016994);
    bool odd = false;
    for( int i = 0; i < 3; i++ )
    {
        odd = odd ^^ (p.x < 0.0) ^^ (p.y < 0.0) ^^ (p.z < 0.0);
        p = abs(p);
        float dn = dot(p, n3);
        if( dn < 0.0 ) { p -= 2.0 * dn * n3; odd = !odd; }
    }
    if( p.y < 0.0 ) { p.y = -p.y; odd = !odd; }
    if( odd ) p.x = -p.x;

    // 2. Farthest face plane = pentagon dual to the nearest snub vertex (v_0, R_3^T v_0, or R_2 v_0).
    float sb = dot(p, vec3(-0.437202286, 0.144773713, -0.055298638));
    float sc = dot(p.xy, vec2(-0.347727210, -0.306999989));
    if( max(sb, sc) > 0.0 )
    {
        if( sb >= sc ) { p.x = -p.x; p -= 2.0 * dot(p, n3) * n3; }
        else           { p.xy = -p.xy; }
    }

    float h = dot(p, vec3(0.173863605, 0.153499994, 0.972732851)) - r;
    if( h <= 0.0 ) return h;

    // 3. Argmax-Face Theorem: exact distance to the single mirror-symmetric pentagon F_0 (3 edges).
    vec2  c  = vec2(dot(p, vec3(0.927503489, -0.357456413, -0.109371800)),
                abs(dot(p, vec3(0.330921025,  0.921228888, -0.204520441))));
    vec2  e0 = vec2(-0.415868539,  0.277630489);          // long edge  C0 -> C1
    vec2  e1 = vec2(-0.251984170, -0.134754073);          // short edge C1 -> C2
    vec2  w0 = c - vec2(0.429422969, 0.0        ) * r;
    vec2  w1 = c - vec2(0.013554430, 0.277630489) * r;
    float xm = c.x + 0.238429740 * r;                     // half short base edge at x = -0.238429740 r
    if( e0.x * w0.y - e0.y * w0.x >= 0.0 &&
        e1.x * w1.y - e1.y * w1.x >= 0.0 && xm >= 0.0 ) return h;

    vec2 d0 = w0 - e0 * clamp(dot(w0, e0) *  3.999594754, 0.0, r);
    vec2 d1 = w1 - e1 * clamp(dot(w1, e1) * 12.246695165, 0.0, r);
    vec2 d2 = vec2(xm, max(c.y - 0.142876417 * r, 0.0));
    return sqrt(h * h + min(min(dot(d0, d0), dot(d1, d1)), dot(d2, d2)));
}

vec4 sdgPentagonalHexecontahedron( vec3 p, float r )
{
    const vec3 n3 = vec3(-0.5, -0.809016994, 0.309016994);
    vec3 s0 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    bool odd = (p.x < 0.0) ^^ (p.y < 0.0) ^^ (p.z < 0.0);
    p = abs(p);
    float dn0 = dot(p, n3); bool r0 = (dn0 < 0.0); if( r0 ) { p -= 2.0 * dn0 * n3; odd = !odd; }

    vec3 s1 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    odd = odd ^^ (p.x < 0.0) ^^ (p.y < 0.0) ^^ (p.z < 0.0);
    p = abs(p);
    float dn1 = dot(p, n3); bool r1 = (dn1 < 0.0); if( r1 ) { p -= 2.0 * dn1 * n3; odd = !odd; }

    vec3 s2 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    odd = odd ^^ (p.x < 0.0) ^^ (p.y < 0.0) ^^ (p.z < 0.0);
    p = abs(p);
    float dn2 = dot(p, n3); bool r2 = (dn2 < 0.0); if( r2 ) { p -= 2.0 * dn2 * n3; odd = !odd; }

    bool sy3 = (p.y < 0.0);
    if( sy3 ) { p.y = -p.y; odd = !odd; }
    if( odd ) p.x = -p.x;

    int   rot = 0;
    float sb  = dot(p, vec3(-0.437202286, 0.144773713, -0.055298638));
    float sc  = dot(p.xy, vec2(-0.347727210, -0.306999989));
    if( max(sb, sc) > 0.0 )
    {
        if( sb >= sc ) { p.x = -p.x; p -= 2.0 * dot(p, n3) * n3; rot = 1; }
        else           { p.xy = -p.xy;                            rot = 2; }
    }

    const vec3 v0 = vec3(0.173863605,  0.153499994,  0.972732851);
    const vec3 b1 = vec3(0.927503489, -0.357456413, -0.109371800);
    const vec3 b2 = vec3(0.330921025,  0.921228888, -0.204520441);
    float h = dot(p, v0) - r;
    float d;
    vec3  g = v0;
    if( h <= 0.0 )
    {
        d = h;
    }
    else
    {
        float cy = dot(p, b2);
        vec2  c  = vec2(dot(p, b1), abs(cy));
        vec2  e0 = vec2(-0.415868539,  0.277630489);
        vec2  e1 = vec2(-0.251984170, -0.134754073);
        vec2  w0 = c - vec2(0.429422969, 0.0        ) * r;
        vec2  w1 = c - vec2(0.013554430, 0.277630489) * r;
        float xm = c.x + 0.238429740 * r;
        if( e0.x * w0.y - e0.y * w0.x >= 0.0 &&
            e1.x * w1.y - e1.y * w1.x >= 0.0 && xm >= 0.0 )
        {
            d = h;
        }
        else
        {
            vec2  d0 = w0 - e0 * clamp(dot(w0, e0) *  3.999594754, 0.0, r);
            vec2  d1 = w1 - e1 * clamp(dot(w1, e1) * 12.246695165, 0.0, r);
            vec2  d2 = vec2(xm, max(c.y - 0.142876417 * r, 0.0));
            float q0 = dot(d0, d0), q1 = dot(d1, d1), q2 = dot(d2, d2);
            vec2  dd = (q0 <= q1 && q0 <= q2) ? d0 : (q1 <= q2) ? d1 : d2;
            dd.y = (cy < 0.0) ? -dd.y : dd.y;
            vec3  v    = h * v0 + dd.x * b1 + dd.y * b2;
            float v2   = dot(v, v);
            float invV = inversesqrt(max(v2, 1e-30));
            d = v2 * invV;
            g = (v2 > 1e-30) ? v * invV : v0;
        }
    }

    if( rot == 1 ) { g -= 2.0 * dot(g, n3) * n3; g.x = -g.x; }
    if( rot == 2 ) { g.xy = -g.xy; }
    if( odd ) g.x = -g.x;
    if( sy3 ) g.y = -g.y;
    if( r2 ) g -= 2.0 * dot(g, n3) * n3; g *= s2;
    if( r1 ) g -= 2.0 * dot(g, n3) * n3; g *= s1;
    if( r0 ) g -= 2.0 * dot(g, n3) * n3; g *= s0;
    return vec4(d, g);
}

#endif // SDF_GEMS_SD_SNUB_DODECAHEDRON_GLSL
