// ============================================================================
// biPonceletCaustic — Integrable Dynamical Billiards, Confocal Poncelet Porism
//                     Polygons, Caustic Invariants & Whispering-Gallery SDFs
//
//   biEllipseCausticParam, biEllipsoidCaustics,
//   biPonceletCaustic3 / 4 / 6, biPonceletStep2D,
//   sdPonceletTriangle2D, sdPonceletPolygon2D,
//   sdWhisperingGallery3D, sdgWhisperingGallery3D,
//   biConfocalCoords2D, sdConfocalWaveguide3D
//
// @author Gemini 4 Argon (sdf-gems, 2026)
// @exact  biEllipseCausticParam / biEllipsoidCaustics: Exact algebraic invariants.
//         sdPonceletTriangle2D / sdPonceletPolygon2D / sdWhisperingGallery3D:
//         Exact Euclidean SDF everywhere outside the tubular segments (d >= -r, and
//         everywhere for r = 0); conservative lower bound inside vertex-sphere overlaps.
//         sdConfocalWaveguide3D: Globally 1-Lipschitz conservative bound (2 sqrts, 0 trig).
// ============================================================================

#ifndef SDF_GEMS_HAS_BI_PONCELET_CAUSTIC
#define SDF_GEMS_HAS_BI_PONCELET_CAUSTIC

// ---------------------------------------------------------------------------
// 1. Jacobi-Chasles Confocal Billiard Caustic Invariants (2D Ellipse & 3D Ellipsoid)
// ---------------------------------------------------------------------------

// Returns the conserved confocal caustic parameter lambda <= a^2 for a 2D line (p, v)
// inside the ellipse x^2/a^2 + y^2/b^2 = 1 (with ab = vec2(a, b), a >= b > 0):
//   (x + t vx)^2 / (a^2 - lambda) + (y + t vy)^2 / (b^2 - lambda) = 1
// has a double root in t iff:
//   lambda = (b^2 vx^2 + a^2 vy^2 - (x vy - y vx)^2) / (vx^2 + vy^2).
// Conserved across every specular bounce off E(a, b) without tracing bounces:
//   - 0 < lambda < b^2 : whispering-gallery orbit tangent to confocal ellipse E_c(sqrt(a^2-lambda), sqrt(b^2-lambda))
//   - lambda == b^2    : focal orbit passing through the foci (+-sqrt(a^2-b^2), 0)
//   - b^2 < lambda < a^2 : bouncing-ball orbit tangent to confocal hyperbola through the focal segment.
float biEllipseCausticParam( in vec2 p, in vec2 v, in vec2 ab )
{
    float a2 = ab.x * ab.x;
    float b2 = ab.y * ab.y;
    float L  = p.x * v.y - p.y * v.x;
    float v2 = max(dot(v, v), 1e-30);
    return (b2 * v.x * v.x + a2 * v.y * v.y - L * L) / v2;
}

// Returns the two conserved Jacobi-Chasles confocal quadric caustic parameters
// vec2(lambda1, lambda2) (with lambda1 <= lambda2) for a 3D line (p, v) inside
// the triaxial ellipsoid x^2/a^2 + y^2/b^2 + z^2/c^2 = 1 (abc = vec3(a, b, c)).
// By Cauchy-Binet on the discriminant of Q_lambda(p + t v) = 1, the cubic term cancels
// leaving the exact quadratic c2 * lambda^2 - c1 * lambda + c0 = 0.
vec2 biEllipsoidCaustics( in vec3 p, in vec3 v, in vec3 abc )
{
    vec3  A  = abc * abc;
    vec3  v2 = v * v;
    vec3  L  = cross(p, v);
    vec3  L2 = L * L;
    float c2 = max(dot(v2, vec3(1.0)), 1e-30);
    float c1 = v2.x * (A.y + A.z) + v2.y * (A.z + A.x) + v2.z * (A.x + A.y) - dot(L2, vec3(1.0));
    float c0 = v2.x * A.y * A.z + v2.y * A.z * A.x + v2.z * A.x * A.y - dot(A, L2);
    float hb = 0.5 * c1 / c2;
    float sq = sqrt(max(hb * hb - c0 / c2, 0.0));
    return vec2(hb - sq, hb + sq);
}

// ---------------------------------------------------------------------------
// 2. Closed-Form Cayley Caustic Invariants & Perimeters for Poncelet N-gons
//    Returns vec4(lambda_N, a_c, b_c, perimeter_N) for N = 3, 4, 6.
// ---------------------------------------------------------------------------

// N = 3 (Confocal Poncelet Triangle family):
// Cancellation-free radicals via D = sqrt(a^4 - a^2 b^2 + b^4).
vec4 biPonceletCaustic3( in vec2 ab )
{
    float a2  = ab.x * ab.x, b2 = ab.y * ab.y;
    float D   = sqrt(max(a2 * (a2 - b2) + b2 * b2, 0.0));
    float s   = a2 + b2 + 2.0 * D;
    float lam = 3.0 * a2 * b2 / max(s, 1e-30);
    float ac  = ab.x * a2 / max(D + b2, 1e-30);
    float bc  = ab.y * b2 / max(D + a2, 1e-30);
    float L3  = 2.0 * (a2 + b2 + D) * sqrt(3.0 / max(s, 1e-30));
    return vec4(lam, ac, bc, L3);
}

// N = 4 (Confocal Poncelet Parallelogram family):
// lambda_4 = a^2 b^2 / (a^2 + b^2), invariant perimeter L_4 = 4 sqrt(a^2 + b^2).
vec4 biPonceletCaustic4( in vec2 ab )
{
    float a2  = ab.x * ab.x, b2 = ab.y * ab.y;
    float s   = a2 + b2;
    float rs  = inversesqrt(max(s, 1e-30));
    float lam = a2 * b2 * (rs * rs);
    return vec4(lam, a2 * rs, b2 * rs, 4.0 * s * rs);
}

// N = 6 (Confocal Poncelet Hexagon family):
// lambda_6 = a^2 b^2 / (a + b)^2, rational invariant perimeter L_6 = 4 (a^2 + a b + b^2) / (a + b).
vec4 biPonceletCaustic6( in vec2 ab )
{
    float a = ab.x, b = ab.y;
    float inv_ab = 1.0 / max(a + b, 1e-30);
    float lam = (a * b * inv_ab) * (a * b * inv_ab);
    float ac  = a * sqrt(max(a * (a + 2.0 * b), 0.0)) * inv_ab;
    float bc  = b * sqrt(max(b * (b + 2.0 * a), 0.0)) * inv_ab;
    float L6  = 4.0 * (a * a + a * b + b * b) * inv_ab;
    return vec4(lam, ac, bc, L6);
}

// Exact 1-sqrt, 0-trig Poncelet vertex map: steps from p0 on E(a,b) along the tangent
// to confocal caustic E_c(ac, bc) in direction dir in {-1.0, +1.0} to the next bounce p1 on E(a,b).
vec2 biPonceletStep2D( in vec2 p0, in vec2 ab, in vec3 lam_ac_bc, in float dir )
{
    float a2 = ab.x * ab.x, b2 = ab.y * ab.y;
    float lam = lam_ac_bc.x, ac = lam_ac_bc.y, bc = lam_ac_bc.z;
    vec2  M  = vec2(b2 * p0.x, a2 * p0.y);
    vec2  T  = vec2(-a2 * p0.y, b2 * p0.x);
    float sp = ab.x * ab.y * sqrt(max(lam, 0.0));
    float cp = sqrt(max(b2 * (bc * bc) * (p0.x * p0.x) + a2 * (ac * ac) * (p0.y * p0.y), 0.0));
    vec2  V  = -sp * M + (dir * cp) * T;
    float den = max(b2 * V.x * V.x + a2 * V.y * V.y, 1e-30);
    float t   = -2.0 * (b2 * p0.x * V.x + a2 * p0.y * V.y) / den;
    return p0 + t * V;
}

// Helper: exact distance from p to segment a..b in 2D
float biSegDist2D( in vec2 p, in vec2 a, in vec2 b )
{
    vec2  pa = p - a, ba = b - a;
    float h  = clamp(dot(pa, ba) / max(dot(ba, ba), 1e-30), 0.0, 1.0);
    return length(pa - ba * h);
}

// ---------------------------------------------------------------------------
// 3. Exact 2D Poncelet Billiard Polygon SDFs (N = 3, 4, 6)
// ---------------------------------------------------------------------------

// Exact Euclidean SDF to the closed 3-bounce Poncelet triangle (tube radius r)
// inscribed in E(a,b) at continuous phase theta0 in [0, 2*pi).
float sdPonceletTriangle2D( in vec2 p, in vec2 ab, in float theta0, in float r )
{
    vec4  c3 = biPonceletCaustic3(ab);
    vec2  p0 = vec2(ab.x * cos(theta0), ab.y * sin(theta0));
    vec2  p1 = biPonceletStep2D(p0, ab, c3.xyz,  1.0);
    vec2  p2 = biPonceletStep2D(p0, ab, c3.xyz, -1.0);
    float d  = min(biSegDist2D(p, p0, p1), min(biSegDist2D(p, p1, p2), biSegDist2D(p, p2, p0)));
    return d - r;
}

// Exact Euclidean SDF to the closed N-bounce Poncelet billiard polygon (n in {3, 4, 6})
// inscribed in E(a,b) at continuous phase theta0, with tube radius r.
// Exploits central symmetry p_{k + N/2} = -p_k for N = 4 (0 steps) and N = 6 (1 step).
float sdPonceletPolygon2D( in vec2 p, in vec2 ab, in float theta0, in int n, in float r )
{
    float c0 = cos(theta0), s0 = sin(theta0);
    vec2  p0 = vec2(ab.x * c0, ab.y * s0);
    if (n <= 3)
    {
        return sdPonceletTriangle2D(p, ab, theta0, r);
    }
    if (n == 4)
    {
        float a2 = ab.x * ab.x, b2 = ab.y * ab.y;
        vec2  w  = vec2(-a2 * s0, b2 * c0);
        vec2  u1 = w * inversesqrt(max(dot(w, w), 1e-30));
        vec2  p1 = ab * u1;
        float d  = min(min(biSegDist2D(p, p0, p1), biSegDist2D(p, p1, -p0)),
                       min(biSegDist2D(p, -p0, -p1), biSegDist2D(p, -p1, p0)));
        return d - r;
    }
    vec4  c6 = biPonceletCaustic6(ab);
    vec2  p1 = biPonceletStep2D(p0, ab, c6.xyz,  1.0);
    vec2  p5 = biPonceletStep2D(p0, ab, c6.xyz, -1.0);
    vec2  p2 = -p5, p3 = -p0, p4 = -p1;
    float d  = min(min(biSegDist2D(p, p0, p1), biSegDist2D(p, p1, p2)),
                   min(min(biSegDist2D(p, p2, p3), biSegDist2D(p, p3, p4)),
                       min(biSegDist2D(p, p4, p5), biSegDist2D(p, p5, p0))));
    return d - r;
}

// ---------------------------------------------------------------------------
// 4. 3D Poncelet Whispering-Gallery Billiard Assembly (Exact outside tubes + sdg)
//    Combines the closed 3-bounce Poncelet triangle (radius r) in the z = 0 plane
//    with the two confocal focus spheres at (+-sqrt(a^2 - b^2), 0, 0) (radius 1.35 * r).
// ---------------------------------------------------------------------------

vec4 biCapsuleSdg3D( in vec3 p, in vec3 a, in vec3 b, in float r )
{
    vec3  pa = p - a, ba = b - a;
    float h  = clamp(dot(pa, ba) / max(dot(ba, ba), 1e-30), 0.0, 1.0);
    vec3  v  = pa - ba * h;
    float l2 = dot(v, v);
    if (l2 > 1e-24)
    {
        float rl = inversesqrt(l2);
        return vec4(l2 * rl - r, v * rl);
    }
    return vec4(-r, 0.0, 0.0, 1.0);
}

vec4 biSphereSdg3D( in vec3 p, in vec3 c, in float r )
{
    vec3  v  = p - c;
    float l2 = dot(v, v);
    if (l2 > 1e-24)
    {
        float rl = inversesqrt(l2);
        return vec4(l2 * rl - r, v * rl);
    }
    return vec4(-r, 0.0, 0.0, 1.0);
}

vec4 sdgWhisperingGallery3D( in vec3 p, in vec2 ab, in float theta0, in float r )
{
    vec4  c3 = biPonceletCaustic3(ab);
    vec2  p0 = vec2(ab.x * cos(theta0), ab.y * sin(theta0));
    vec2  p1 = biPonceletStep2D(p0, ab, c3.xyz,  1.0);
    vec2  p2 = biPonceletStep2D(p0, ab, c3.xyz, -1.0);

    vec4  dg = biCapsuleSdg3D(p, vec3(p0, 0.0), vec3(p1, 0.0), r);
    vec4  d1 = biCapsuleSdg3D(p, vec3(p1, 0.0), vec3(p2, 0.0), r);
    if (d1.x < dg.x) dg = d1;
    vec4  d2 = biCapsuleSdg3D(p, vec3(p2, 0.0), vec3(p0, 0.0), r);
    if (d2.x < dg.x) dg = d2;

    float c  = sqrt(max((ab.x - ab.y) * (ab.x + ab.y), 0.0));
    float rf = 1.35 * r;
    vec4  f1 = biSphereSdg3D(p, vec3(-c, 0.0, 0.0), rf);
    if (f1.x < dg.x) dg = f1;
    vec4  f2 = biSphereSdg3D(p, vec3( c, 0.0, 0.0), rf);
    if (f2.x < dg.x) dg = f2;
    return dg;
}

float sdWhisperingGallery3D( in vec3 p, in vec2 ab, in float theta0, in float r )
{
    return sdgWhisperingGallery3D(p, ab, theta0, r).x;
}

// ---------------------------------------------------------------------------
// 5. Confocal Elliptic Coordinates & 1-Lipschitz Whispering-Gallery Waveguide
//    mu = 0.5 * (r1 + r2) >= c,  nu = 0.5 * (r1 - r2) in [-c, c]  (2 sqrts, 0 trig).
//    Returns vec4(mu, nu, ||grad mu||, ||grad nu||).
// ---------------------------------------------------------------------------

vec4 biConfocalCoords2D( in vec2 p, in float c )
{
    float r1  = length(vec2(p.x + c, p.y));
    float r2  = length(vec2(p.x - c, p.y));
    float mu  = 0.5 * (r1 + r2);
    float nu  = 0.5 * (r1 - r2);
    float den = max((mu - nu) * (mu + nu), 1e-30);
    float gmu = sqrt(max((mu - c) * (mu + c), 0.0) / den);
    float gnu = sqrt(max((c - nu) * (c + nu), 0.0) / den);
    return vec4(mu, nu, gmu, gnu);
}

// Globally 1-Lipschitz conservative SDF for a 3D confocal whispering-gallery
// waveguide ring bounded between inner caustic ellipse E_c(sqrt(a^2 - lambda_c), sqrt(b^2 - lambda_c))
// and outer billiard ellipse E(a, b), extruded to |z| <= h and rounded by r.
// Since ||grad mu|| = sqrt((mu^2 - c^2)/(mu^2 - nu^2)) <= 1 everywhere, both (mu - a) and (a_c - mu)
// are globally 1-Lipschitz in R^2 with only 2 sqrts (no quartic ellipse solver required).
float sdConfocalWaveguide3D( in vec3 p, in vec2 ab, in float lambda_c, in float h, in float r )
{
    float c   = sqrt(max((ab.x - ab.y) * (ab.x + ab.y), 0.0));
    float r1  = length(vec2(p.x + c, p.y));
    float r2  = length(vec2(p.x - c, p.y));
    float mu  = 0.5 * (r1 + r2);
    float ac  = sqrt(max(ab.x * ab.x - lambda_c, c * c));
    float dxy = max(mu - ab.x, ac - mu);
    vec2  w   = vec2(dxy, abs(p.z) - h);
    return min(max(w.x, w.y), 0.0) + length(max(w, 0.0)) - r;
}

#endif // SDF_GEMS_HAS_BI_PONCELET_CAUSTIC
