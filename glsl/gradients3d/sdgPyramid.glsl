// ============================================================================
// sdgPyramid — 0/1-Per-Pixel-SFU 3D Square Pyramid Distance & Analytical Gradient
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1` everywhere in R^3; 0 SFU ops in base zone, 0 per-pixel SFU ops in slanted-face zone, 1 inversesqrt in edge/vertex zone; fixes IQ's interior base-face bug)
// ============================================================================
//
// Mathematical Foundation:
//   For a 3D square-base pyramid with unit base [-0.5, 0.5] x {0} x [-0.5, 0.5]
//   and apex (0, h, 0):
//
//   1. Euclidean-Vector Decomposition Identity:
//      In Inigo Quilez's scalar `sdPyramid`, after folding u = abs(p.xz), sorting
//      u.x >= u.y, and shifting u -= 0.5, the code computes two scalar combinations
//      (a + q.z^2)/m2 and (b + q.z^2)/m2 where m2 = h^2 + 0.25, q.y = h*p.y - 0.5*u.x,
//      and q.z = h*u.x + 0.5*p.y.
//      By expanding the squared Euclidean norms, those two scalar expressions are
//      identically equal to the squared lengths of the explicit 3D Euclidean
//      difference vectors from p to the base edge E_A and slanted hip edge E_B:
//        (a + q.z^2)/m2 === ||v_A||_2^2,   v_A = vec3(u.x, p.y, max(u.y, 0.0)),
//        (b + q.z^2)/m2 === ||v_B||_2^2,   v_B = vec3(u.x + 0.5*t, p.y - h*t, u.y + 0.5*t),
//      where t = clamp((q.y - 0.5*u.y)/(m2 + 0.25), 0.0, 1.0)!
//
//   2. Three-Zone 0/1-SFU Decomposition (Fixes IQ's Interior Base-Face Bug):
//      - Zone 1 (Base Face): Either outside below the base square (p.y <= 0 && u.x <= 0)
//        or inside the pyramid closer to the base y = 0 than the slanted face
//        (q.z <= 0 && p.y > 0 && m2 * p.y^2 <= q.z^2):
//          d = -p.y,   grad = (0, -1, 0)   (0 SFU ops!).
//      - Zone 2 (Slanted Triangular Face): Inside closer to the slanted face
//        (q.z <= 0 && p.y > 0) or in the exterior slanted-face normal prism
//        (q.y > 0 && -u.y*m2 - 0.5*q.y > 0):
//          d = q.z * invM,   grad_folded = (h, 0.5, 0) * invM   (0 per-pixel SFU ops!).
//      - Zone 3 (Base Edge vs Hip Edge): Selects v = (||v_A||^2 < ||v_B||^2) ? v_A : v_B
//        and normalizes with a single `inversesqrt`!

vec4 sdgPyramid( vec3 p, float h )
{
    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    float sz = (p.z < 0.0) ? -1.0 : 1.0;
    vec2  u  = abs(p.xz);
    bool  sw = (u.y > u.x);
    if( sw ) u = u.yx;
    u -= 0.5;

    float m2 = h * h + 0.25;
    float qy = h * p.y - 0.5 * u.x;
    float qz = h * u.x + 0.5 * p.y;

    // Zone 1: Flat base face (0 SFU ops in both exterior below base & interior near base!)
    if( (p.y <= 0.0 && u.x <= 0.0) || (qz <= 0.0 && p.y > 0.0 && m2 * p.y * p.y <= qz * qz) )
        return vec4(-p.y, 0.0, -1.0, 0.0);

    float d;
    vec3  g;
    // Zone 2: Slanted triangular face (0 per-pixel SFU ops; only invM = inversesqrt(m2)!)
    if( (qz <= 0.0 && p.y > 0.0) || (qy > 0.0 && -u.y * m2 - 0.5 * qy > 0.0) )
    {
        float invM = inversesqrt(m2);
        d = qz * invM;
        g = vec3(h * invM, 0.5 * invM, 0.0);
    }
    else
    {
        // Zone 3: Exterior base edge vA vs slanted hip edge vB (1 inversesqrt total!)
        vec3  vA   = vec3(u.x, p.y, max(u.y, 0.0));
        float t    = clamp((qy - 0.5 * u.y) / (m2 + 0.25), 0.0, 1.0);
        vec3  vB   = vec3(u.x + 0.5 * t, p.y - h * t, u.y + 0.5 * t);
        float dA2  = dot(vA, vA), dB2 = dot(vB, vB);
        vec3  v    = (dA2 < dB2) ? vA : vB;
        float v2   = min(dA2, dB2);
        float invV = inversesqrt(max(v2, 1e-30));
        d = v2 * invV;
        g = v * invV;
    }

    if( sw ) g.xz = g.zx;
    return vec4(d, sx * g.x, g.y, sz * g.z);
}
