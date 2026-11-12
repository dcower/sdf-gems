// ============================================================================
// sdPentagonalIcositetrahedron & sdgPentagonalIcositetrahedron — Exact 0/1-Sqrt Chiral Catalan
// Pentagonal Icositetrahedron (24 Pentagons, 38 Vertices, 60 Edges)
//
// @author  Claude Opus 5.5 (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^3; 0 sqrts in interior & all 24 face prisms,
//          1 sqrt in the edge/vertex cones)
// @novelty Exact Euclidean SDF and analytical gradient for the chiral Catalan pentagonal
//          icositetrahedron (dual of the snub cube) via the Chiral Parity Fold (O_h abs + sort with
//          reflection-parity tracking, odd parity undone by one swap of y,z) followed by the
//          Argmax-Face Lemma: the closest surface point always lies on the pentagon whose plane is
//          farthest, i.e. the pentagon dual to the closest snub-cube vertex. In the folded domain
//          that is one of 3 pentagons, each mapped to a fixed pentagon by a rotation swizzle;
//          the pentagon's own in-plane mirror axis reduces it to 3 edges.
//
// Chirality: dual of sdSnubCube (same handedness). Mirror enantiomer:
//   sdPentagonalIcositetrahedron(p * vec3(-1.0, 1.0, 1.0), r).
//
// Parameters:
//   p : 3D query point
//   r : insphere radius tangent to all 24 pentagonal faces
//       (6 four-valent vertices at radius 1.17599990 * r, 32 three-valent at 1.10743557 * r;
//        long edge 0.72776796 * r, short edge 0.51264140 * r)
// ============================================================================

#ifndef SDF_GEMS_SD_PENTAGONAL_ICOSITETRAHEDRON_GLSL
#define SDF_GEMS_SD_PENTAGONAL_ICOSITETRAHEDRON_GLSL

float sdPentagonalIcositetrahedron( vec3 p, float r )
{
    // 1. Chiral parity fold into the rotation-group fundamental domain D = C u s_yz(C).
    bool odd = (p.x < 0.0) ^^ (p.y < 0.0) ^^ (p.z < 0.0);
    p = abs(p);
    if( p.x < p.y ) { p.xy = p.yx; odd = !odd; }
    if( p.y < p.z ) { p.yz = p.zy; odd = !odd; }
    if( p.x < p.y ) { p.xy = p.yx; odd = !odd; }
    if( odd ) p.yz = p.zy;

    // 2. Farthest face plane = pentagon dual to the nearest snub vertex v0, v16 or v8.
    float s16 = dot(p, vec3(-0.45631099, 0.70440226, -0.24809127));
    float s8  = dot(p, vec3( 0.0,        0.24809127, -0.83928675));
    if( max(s16, s8) > 0.0 ) p = (s16 >= s8) ? p.yzx : vec3(p.x, -p.z, p.y);

    float h = dot(p, vec3(0.85034021, 0.25135865, 0.46232063)) - r;
    if( h <= 0.0 ) return h;

    // 3. Exact distance to the fixed pentagon (mirror-folded in its plane: 3 edges).
    vec2 c = vec2(dot(p, vec3(0.52623334, -0.40617032, -0.74706368)),
                  abs(dot(p, vec3(0.0,         0.87854682, -0.47765625))));
    vec2  e0 = vec2(-0.55442194,  0.47144726);           // long edge  C0 -> C1
    vec2  e1 = vec2(-0.46531899, -0.21512657);           // short edge C1 -> C2
    vec2  w0 = c - vec2( 0.61885036, 0.0       ) * r;
    vec2  w1 = c - vec2( 0.06442842, 0.47144726) * r;
    float xm = c.x + 0.40089056 * r;                     // half short edge C2 -> (x2, 0)
    if( e0.x * w0.y - e0.y * w0.x >= 0.0 &&
        e1.x * w1.y - e1.y * w1.x >= 0.0 && xm >= 0.0 ) return h;

    vec2 d0 = w0 - e0 * clamp(dot(w0, e0) * 1.88805279, 0.0, r);
    vec2 d1 = w1 - e1 * clamp(dot(w1, e1) * 3.80515766, 0.0, r);
    vec2 d2 = vec2(xm, max(c.y - 0.25632070 * r, 0.0));
    return sqrt(h * h + min(min(dot(d0, d0), dot(d1, d1)), dot(d2, d2)));
}

vec4 sdgPentagonalIcositetrahedron( vec3 p, float r )
{
    vec3 sg = vec3((p.x < 0.0) ? -1.0 : 1.0,
                   (p.y < 0.0) ? -1.0 : 1.0,
                   (p.z < 0.0) ? -1.0 : 1.0);
    bool odd = (p.x < 0.0) ^^ (p.y < 0.0) ^^ (p.z < 0.0);
    p = abs(p);
    bool s1 = (p.x < p.y); if( s1 ) { p.xy = p.yx; odd = !odd; }
    bool s2 = (p.y < p.z); if( s2 ) { p.yz = p.zy; odd = !odd; }
    bool s3 = (p.x < p.y); if( s3 ) { p.xy = p.yx; odd = !odd; }
    if( odd ) p.yz = p.zy;

    int   rot = 0;
    float s16 = dot(p, vec3(-0.45631099, 0.70440226, -0.24809127));
    float s8  = dot(p, vec3( 0.0,        0.24809127, -0.83928675));
    if( max(s16, s8) > 0.0 )
    {
        if( s16 >= s8 ) { p = p.yzx;                rot = 1; }
        else            { p = vec3(p.x, -p.z, p.y); rot = 2; }
    }

    const vec3 n0 = vec3(0.85034021,  0.25135865,  0.46232063);
    const vec3 b1 = vec3(0.52623334, -0.40617032, -0.74706368);
    const vec3 b2 = vec3(0.0,         0.87854682, -0.47765625);
    float h = dot(p, n0) - r;
    float d;
    vec3  g = n0;
    if( h <= 0.0 )
    {
        d = h;
    }
    else
    {
        float cy = dot(p, b2);
        vec2  c  = vec2(dot(p, b1), abs(cy));
        vec2  e0 = vec2(-0.55442194,  0.47144726);
        vec2  e1 = vec2(-0.46531899, -0.21512657);
        vec2  w0 = c - vec2( 0.61885036, 0.0       ) * r;
        vec2  w1 = c - vec2( 0.06442842, 0.47144726) * r;
        float xm = c.x + 0.40089056 * r;
        if( e0.x * w0.y - e0.y * w0.x >= 0.0 &&
            e1.x * w1.y - e1.y * w1.x >= 0.0 && xm >= 0.0 )
        {
            d = h;
        }
        else
        {
            vec2  d0 = w0 - e0 * clamp(dot(w0, e0) * 1.88805279, 0.0, r);
            vec2  d1 = w1 - e1 * clamp(dot(w1, e1) * 3.80515766, 0.0, r);
            vec2  d2 = vec2(xm, max(c.y - 0.25632070 * r, 0.0));
            float q0 = dot(d0, d0), q1 = dot(d1, d1), q2 = dot(d2, d2);
            vec2  dd = (q0 <= q1 && q0 <= q2) ? d0 : (q1 <= q2) ? d1 : d2;
            dd.y = (cy < 0.0) ? -dd.y : dd.y;
            vec3  v    = h * n0 + dd.x * b1 + dd.y * b2;
            float v2   = dot(v, v);
            float invV = inversesqrt(max(v2, 1e-30));
            d = v2 * invV;
            g = (v2 > 1e-30) ? v * invV : n0;
        }
    }

    if( rot == 1 ) g = g.zxy;
    if( rot == 2 ) g = vec3(g.x, g.z, -g.y);
    if( odd ) g.yz = g.zy;
    if( s3 ) g.xy = g.yx;
    if( s2 ) g.yz = g.zy;
    if( s1 ) g.xy = g.yx;
    return vec4(d, sg * g);
}

#endif // SDF_GEMS_SD_PENTAGONAL_ICOSITETRAHEDRON_GLSL
