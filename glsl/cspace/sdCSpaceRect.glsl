// ============================================================================
// sdCSpaceRectRectSlice, sdgCSpaceRectRectPose, sdCSpaceRectRect, sdgCSpaceRectRect,
// stepCSpaceRectRect, stepCSpaceRectRectCA, sdgCSpaceRectRectPivot
// — Configuration-space obstacle of a rotating rectangle vs. a rectangle
//   ("piano mover" C-obstacle) as a sphere-traceable solid.
//
// Robot A: rectangle, half-extents a, centred at (x, y), rotated by theta.
// Obstacle B: axis-aligned rectangle, half-extents b, centred at the origin.
// F(x, y, theta) = signed separation of the two bodies:
//   F > 0: Euclidean distance between A and B; F < 0: -(penetration depth)
//   (minimum translation that separates them); F = 0: contact.
// For fixed theta, F(., ., theta) is the exact 2D SDF of the Minkowski sum
// B (+) A_theta, a centrally symmetric convex octagon.
//
// In pose space q = (x, y, z = rho * theta), F is Lipschitz with the TIGHT constant
//   L = sqrt(1 + (max(a.x, a.y) / rho)^2)
// (the torque arm |dF/dtheta| is bounded by the largest derivative of the robot's
// support function, max(a.x, a.y), not by its circumradius length(a)).
//
// @author Claude Fable 5.2 (sdf-gems, 2026)
// @exact  sdCSpaceRectRectSlice / sdgCSpaceRectRectPose: Exact (2D SDF in (x, y) for
//         fixed theta; exact partial derivatives). sdCSpaceRectRect / sdgCSpaceRectRect:
//         Conservative bound in q (F / L, tight Lipschitz constant, NOT a Euclidean
//         SDF of the 3D solid). stepCSpaceRectRect*: conservative ray step lengths.
// ============================================================================

#ifndef SDF_GEMS_HAS_CSPACE_RECT
#define SDF_GEMS_HAS_CSPACE_RECT

// One boundary segment of the octagon. w: query point, S: segment start, dir: unit
// direction, len: length, nrm: outward normal, bA: for robot-edge segments the
// obstacle corner they are attached to (useB = 1), otherwise aFix is the robot
// vertex (relative to the robot centre) the obstacle edge is attached to.
// Keeps the nearest candidate in (best = vec4(d2, sigma, 0, 0), diff, aRel, nBest).
void cspaceRectSeg( in vec2 w, in float sigma, in vec2 S, in vec2 dir, in float len, in vec2 nrm,
                    in vec2 bA, in vec2 aFix, in float useB,
                    inout float bestD2, inout vec2 bestDiff, inout vec2 bestA, inout vec2 bestN )
{
    float t  = clamp(dot(w - S, dir), 0.0, len);
    vec2  m  = S + t * dir;
    vec2  df = w - m;
    float d2 = dot(df, df);
    if( d2 < bestD2 )
    {
        bestD2   = d2;
        bestDiff = sigma * df;
        bestA    = sigma * ((useB > 0.5) ? (m - bA) : aFix);
        bestN    = sigma * nrm;
    }
}

// 1. Exact signed separation + exact pose-space partial derivatives.
//    Returns vec4(F, dF/dx, dF/dy, dF/dtheta):
//      (dF/dx, dF/dy) = unit contact normal n (points from the obstacle to the robot),
//      dF/dtheta      = r x n = r.x * n.y - r.y * n.x, the torque arm of the robot's
//                       witness point r (relative to the robot centre).
//    0 divisions on the distance path, 1 sqrt + 1 division for the normal.
vec4 sdgCSpaceRectRectPose( in vec2 p, in float theta, in vec2 a, in vec2 b )
{
    // Reduce theta to [0, pi/2): a rectangle turned by a quarter is the same
    // rectangle with swapped half-extents.
    float k  = floor(theta * 0.63661977236758134308);
    float th = theta - k * 1.57079632679489661923;
    if( mod(k, 2.0) > 0.5 ) a = a.yx;
    float c = cos(th), s = sin(th);
    vec2  u = vec2(c, s), v = vec2(-s, c);      // robot axes

    vec2 A0 = a.x * u - a.y * v;                // robot vertices (relative to its centre)
    vec2 A1 = a.x * u + a.y * v;
    vec2 B0 = b;                                // obstacle corners
    vec2 B1 = vec2(-b.x, b.y);

    float bestD2 = 3.0e38;
    vec2  bestDiff = vec2(0.0), bestA = vec2(0.0), bestN = vec2(1.0, 0.0);
    for( int i = 0; i < 2; i++ )
    {
        float sg = (i == 0) ? 1.0 : -1.0;       // central symmetry: test p and -p
        vec2  w  = sg * p;
        cspaceRectSeg(w, sg, B0 + A0, v,                2.0 * a.y, u,               B0, A0,  1.0, bestD2, bestDiff, bestA, bestN);
        cspaceRectSeg(w, sg, B0 + A1, vec2(-1.0, 0.0),  2.0 * b.x, vec2(0.0, 1.0),  B0, A1,  0.0, bestD2, bestDiff, bestA, bestN);
        cspaceRectSeg(w, sg, B1 + A1, -u,               2.0 * a.x, v,               B1, A1,  1.0, bestD2, bestDiff, bestA, bestN);
        cspaceRectSeg(w, sg, B1 - A0, vec2(0.0, -1.0),  2.0 * b.y, vec2(-1.0, 0.0), B1, -A0, 0.0, bestD2, bestDiff, bestA, bestN);
    }

    // Interior test and interior distance from the 4 pairs of supporting lines:
    // h(n) = b . |n| + a . |R^T n| is the support function of the octagon.
    float g = max(max(abs(p.x) - (b.x + a.x * c + a.y * s), abs(p.y) - (b.y + a.x * s + a.y * c)),
                  max(abs(dot(p, u)) - (a.x + b.x * c + b.y * s), abs(dot(p, v)) - (a.y + b.x * s + b.y * c)));
    bool  inside = g < 0.0;
    float dist   = sqrt(bestD2);
    float F      = inside ? g : dist;
    vec2  n      = (inside || bestD2 < 1e-12) ? bestN : (bestDiff / dist);
    // bestA is the robot's part of the nearest octagon point; the robot's witness point
    // relative to its centre is r = -bestA, and dF/dtheta = r x n.
    return vec4(F, n, n.x * bestA.y - n.y * bestA.x);
}

// 2. Exact 2D SDF of the C-obstacle slice at fixed theta (value only).
float sdCSpaceRectRectSlice( in vec2 xy, in float theta, in vec2 a, in vec2 b )
{
    return sdgCSpaceRectRectPose(xy, theta, a, b).x;
}

// 3. Tight Lipschitz constant of F in q = (x, y, rho * theta): |grad_xy F| = 1 and
//    |dF/dtheta| <= max(a.x, a.y), attained, so L = sqrt(1 + (max(a.x, a.y) / rho)^2).
float cspaceRectLipschitz( in vec2 a, in float rho )
{
    float m = max(a.x, a.y) / rho;
    return sqrt(1.0 + m * m);
}

// 4. Conservative sphere-traceable bound of the 3D C-obstacle in q = (x, y, rho * theta).
float sdCSpaceRectRect( in vec3 q, in vec2 a, in vec2 b, in float rho )
{
    return sdgCSpaceRectRectPose(q.xy, q.z / rho, a, b).x / cspaceRectLipschitz(a, rho);
}

//    Baseline for comparison only: the textbook constant with the circumradius
//    r_max = length(a) ("every robot point moves at most r_max |dtheta|").
float sdCSpaceRectRectRmax( in vec3 q, in vec2 a, in vec2 b, in float rho )
{
    float m = length(a) / rho;
    return sdgCSpaceRectRectPose(q.xy, q.z / rho, a, b).x / sqrt(1.0 + m * m);
}

// 5. Conservative bound + exact unit normal of the level set of F through q
//    (exact surface normal of the C-obstacle on F = 0).
vec4 sdgCSpaceRectRect( in vec3 q, in vec2 a, in vec2 b, in float rho )
{
    vec4 f = sdgCSpaceRectRectPose(q.xy, q.z / rho, a, b);
    vec3 g = vec3(f.yz, f.w / rho);
    return vec4(f.x / cspaceRectLipschitz(a, rho), g * inversesqrt(dot(g, g)));
}

// 6. Direction-aware conservative step along a unit ray direction rd in q-space:
//    F is 1-Lipschitz in the anisotropic norm |dxy| + max(a.x, a.y) |dtheta|, so the
//    ray cannot reach F = 0 before t = F / (|rd.xy| + max(a.x, a.y) |rd.z| / rho).
//    Always >= sdCSpaceRectRect. Valid inside and outside (negative inside).
float stepCSpaceRectRect( in vec3 q, in vec3 rd, in vec2 a, in vec2 b, in float rho )
{
    float F = sdgCSpaceRectRectPose(q.xy, q.z / rho, a, b).x;
    return F / max(length(rd.xy) + max(a.x, a.y) * abs(rd.z) / rho, 1e-30);
}

// 7. Conservative advancement along the current separating axis n (exterior only):
//    the gap along the fixed axis n closes at a rate of at most
//    -n . rd.xy + max(a.x, a.y) |rd.z| / rho, so the bodies cannot touch before
//    t = F / that rate. Unbounded when the motion opens the gap: capped at tMax.
//    Always >= stepCSpaceRectRect for F > 0; returns F (negative) inside.
float stepCSpaceRectRectCA( in vec3 q, in vec3 rd, in vec2 a, in vec2 b, in float rho, in float tMax )
{
    vec4  f    = sdgCSpaceRectRectPose(q.xy, q.z / rho, a, b);
    float rate = -dot(f.yz, rd.xy) + max(a.x, a.y) * abs(rd.z) / rho;
    return (f.x <= 0.0) ? f.x : ((f.x < rate * tMax) ? (f.x / rate) : tMax);
}

// 8. Offset pivot: the robot rotates about a pivot at local offset -c from its centre,
//    i.e. the pose (x, y) is the pivot and the robot centre is xy + R(theta) c.
//    Returns vec4(F, dF/dx, dF/dy, dF/dtheta) exactly as in (1).
vec4 sdgCSpaceRectRectPivot( in vec2 xy, in float theta, in vec2 a, in vec2 c, in vec2 b )
{
    float co = cos(theta), si = sin(theta);
    vec2  rc = vec2(co * c.x - si * c.y, si * c.x + co * c.y);
    vec4  f  = sdgCSpaceRectRectPose(xy + rc, theta, a, b);
    return vec4(f.xyz, f.w + rc.x * f.z - rc.y * f.y);
}

//    Tight bound of |dF/dtheta| for the offset pivot (reduces to max(a.x, a.y) at c = 0):
//    the largest tangential coordinate of an edge point, plus a whole vertex radius
//    when the pivot lies outside the slab of the opposite pair of edges.
float cspaceRectTorqueBound( in vec2 a, in vec2 c )
{
    vec2  e = a + abs(c);
    float L = max(e.x, e.y);
    if( abs(c.x) > a.x ) L = max(L, length(vec2(abs(c.x) - a.x, e.y)));
    if( abs(c.y) > a.y ) L = max(L, length(vec2(e.x, abs(c.y) - a.y)));
    return L;
}

#endif // SDF_GEMS_HAS_CSPACE_RECT
