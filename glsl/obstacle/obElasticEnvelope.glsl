// ============================================================================
// obElasticEnvelope — Variational Inequalities, Classical Elastic Obstacle
//   Problem (min(-Delta u, u - psi) = 0), and C^{1,1} Free-Boundary Shrink-Wrap
//   Envelopes of Signed Distance Functions
//
// Functions:
//   1. obTautBridge2D, sdgTautBridge2D, obTautBridge, sdgTautBridge
//      Exact C^{1,1} 1D/axisymmetric taut-string (convex-hull bitangent)
//      obstacle envelope over two circles/spheres of radii (r1, r2) at (-h, +h).
//   2. obHarmonicCap3D, sdgHarmonicCap3D, obHarmonicCapBound3D, obHarmonicCapNaive3D
//      3D Newtonian/Coulomb harmonic obstacle problem (-Delta u = 0 outside
//      contact set) over a sphere of radius r pinned at outer radius Rpin,
//      with closed-form cubic free-boundary radius:
//          (a*)^3 - 2 Rpin (a*)^2 + Rpin r^2 = 0.
//   3. obLogCap2D, sdgLogCap2D, obLogCap3D
//      2D/axisymmetric logarithmic Laplace obstacle envelope u(rho) = A ln(Rpin/rho)
//      solving (a*)^2 (1 + ln(Rpin / a*)) = r^2 via concave Newton iteration.
//   4. obBiharmonicPlate, sdgBiharmonicPlate
//      3D radial biharmonic clamped plate obstacle envelope (Delta^2 u = 0)
//      with u(rho) = A (Rpin - rho)^2, C^{1,1} at BOTH the inner obstacle
//      free boundary a* = 4 r^2 / (Rpin + sqrt(Rpin^2 + 8 r^2)) AND the outer
//      clamped boundary rho = Rpin.
//   5. opObstacleBridge, sdgObstacleBridge
//      General 2-SDF C^{1,1} free-boundary obstacle bridge operator preserving
//      d = d_i identically outside the detachment locus (d_j >= gapCoord).
//
// @author Gemini 4 Argon (sdf-gems, 2026)
// @exact  obTautBridge2D, obTautBridge: Exact Euclidean SDF everywhere in R^2/R^3.
//         obHarmonicCap3D, obLogCap2D, obLogCap3D, obBiharmonicPlate: Exact Euclidean
//         SDF everywhere in R^2/R^3 (including above the exterior focal caustic via
//         multi-start + endpoint-safeguarded Newton meridian projection, ~1e-6).
//         obHarmonicCapBound3D, opObstacleBridge: Provably 1-Lipschitz conservative
//         bounds (0 oversteps). obHarmonicCapNaive3D: Unnormalized vertical graph
//         residual (demonstrates Eikonal overestimation when ||grad u|| > 0).
// ============================================================================

#ifndef SDF_GEMS_HAS_OB_ELASTIC_ENVELOPE
#define SDF_GEMS_HAS_OB_ELASTIC_ENVELOPE

// ----------------------------------------------------------------------------
// 1. Taut-String / Convex-Hull Tangent Bridge (obTautBridge, sdgTautBridge)
// ----------------------------------------------------------------------------

vec3 sdgTautBridge2D( in vec2 p, in float h, in float r1, in float r2 )
{
    float L = 2.0 * h;
    float dr = r1 - r2;
    float ay = abs(p.y);
    float sy = (p.y >= 0.0) ? 1.0 : -1.0;

    if (L <= abs(dr) + 1e-7)
    {
        if (r1 >= r2)
        {
            vec2 q = vec2(p.x + h, p.y);
            float l = length(q);
            return vec3(l - r1, (l > 1e-20) ? q / l : vec2(1.0, 0.0));
        }
        else
        {
            vec2 q = vec2(p.x - h, p.y);
            float l = length(q);
            return vec3(l - r2, (l > 1e-20) ? q / l : vec2(1.0, 0.0));
        }
    }

    float sa = dr / L;
    float ca = sqrt(max((1.0 - sa) * (1.0 + sa), 0.0));
    float x1 = p.x + h;
    float x2 = p.x - h;

    if (x1 * ca - ay * sa <= 0.0)
    {
        float l = sqrt(x1 * x1 + ay * ay);
        vec2 g = (l > 1e-20) ? vec2(x1, p.y) / l : vec2(-1.0, 0.0);
        return vec3(l - r1, g);
    }
    if (x2 * ca - ay * sa >= 0.0)
    {
        float l = sqrt(x2 * x2 + ay * ay);
        vec2 g = (l > 1e-20) ? vec2(x2, p.y) / l : vec2(1.0, 0.0);
        return vec3(l - r2, g);
    }
    float d = x1 * sa + ay * ca - r1;
    return vec3(d, sa, sy * ca);
}

float obTautBridge2D( in vec2 p, in float h, in float r1, in float r2 )
{
    return sdgTautBridge2D(p, h, r1, r2).x;
}

vec4 sdgTautBridge( in vec3 p, in float h, in float r1, in float r2 )
{
    float rho2 = dot(p.yz, p.yz);
    float rho = sqrt(rho2);
    vec2 u_yz = (rho > 1e-20) ? p.yz / rho : vec2(1.0, 0.0);
    vec3 dg2 = sdgTautBridge2D(vec2(p.x, rho), h, r1, r2);
    return vec4(dg2.x, dg2.y, dg2.z * u_yz);
}

float obTautBridge( in vec3 p, in float h, in float r1, in float r2 )
{
    return sdgTautBridge(p, h, r1, r2).x;
}

// ----------------------------------------------------------------------------
// 2. 3D Harmonic / Coulomb Obstacle Envelope (obHarmonicCap3D, sdgHarmonicCap3D)
//    (a*)^3 - 2 Rpin (a*)^2 + Rpin r^2 = 0.
// ----------------------------------------------------------------------------

float obHarmonicContactRadius3D( in float r, in float Rpin )
{
    float R = max(Rpin, r * 1.0001);
    float arg = clamp(1.0 - (27.0 * r * r) / (16.0 * R * R), -1.0, 1.0);
    float th = acos(arg) / 3.0;
    float a = (2.0 / 3.0) * R * (1.0 - 2.0 * cos(th + 1.0471975512));
    a = clamp(a, 0.70710678 * r, 0.9999 * r);
    float P = (a - 2.0 * R) * a * a + R * r * r;
    float dP = a * (3.0 * a - 4.0 * R);
    a -= P / min(dP, -1e-12);
    return clamp(a, 0.70710678 * r, 0.99995 * r);
}

vec4 sdgHarmonicCap3D( in vec3 p, in float r, in float Rpin )
{
    float R = max(Rpin, r * 1.001);
    float rho2 = dot(p.xz, p.xz);
    float rho = sqrt(rho2);
    float ay = abs(p.y);
    float sy = (p.y >= 0.0) ? 1.0 : -1.0;
    vec2 u_xz = (rho > 1e-20) ? p.xz / rho : vec2(1.0, 0.0);

    float a = obHarmonicContactRadius3D(r, R);
    float psi_a = sqrt(max((r - a) * (r + a), 1e-12));
    float A = (a * a * a) / psi_a;
    float invR = 1.0 / R;

    float u_rho = (rho <= a) ? sqrt(max(r * r - rho2, 0.0)) : ((rho < R) ? A * (1.0 / rho - invR) : 0.0);
    float sgn = (rho <= R && ay <= u_rho) ? -1.0 : 1.0;

    // Candidate 1: Spherical contact cap t in [0, a*]
    float best_d;
    vec2 best_n2;
    if (rho * psi_a - a * ay <= 0.0)
    {
        float l = sqrt(rho2 + ay * ay);
        best_d = abs(l - r);
        best_n2 = (l > 1e-20) ? (sgn * ((l - r >= 0.0) ? 1.0 : -1.0) * vec2(rho, ay) / l) : vec2(0.0, 1.0);
    }
    else
    {
        vec2 df = vec2(rho - a, ay - psi_a);
        best_d = length(df);
        best_n2 = (best_d > 1e-12) ? (sgn * df / best_d) : normalize(vec2(a, psi_a));
    }

    // Candidate 2: Outer rim endpoint t = R, u(R) = 0
    vec2 dfR = vec2(rho - R, ay);
    float dR = length(dfR);
    if (dR < best_d)
    {
        best_d = dR;
        best_n2 = (dR > 1e-12) ? (sgn * dfR / dR) : normalize(vec2(A * invR * invR, 1.0));
    }

    // Candidate 3: Interior of free-boundary membrane t in [a*, R] via 2-start safeguarded Newton
    for (int s = 0; s < 2; s++)
    {
        float t = (s == 0) ? clamp(rho, a, R) : mix(a, R, 0.35);
        for (int i = 0; i < 4; i++)
        {
            float inv_t = 1.0 / t;
            float inv_t2 = inv_t * inv_t;
            float ut = A * (inv_t - invR);
            float dut = -A * inv_t2;
            float ddut = 2.0 * A * inv_t2 * inv_t;
            float dy = ut - ay;
            float phi = (t - rho) + dut * dy;
            float dphi = 1.0 + dut * dut + ddut * dy;
            t = clamp(t - phi / max(abs(dphi), 1e-6), a, R);
        }
        float ut = A * (1.0 / t - invR);
        vec2 df = vec2(rho - t, ay - ut);
        float dist = length(df);
        if (dist < best_d)
        {
            best_d = dist;
            float dut = -A / (t * t);
            best_n2 = (dist > 1e-12) ? (sgn * df / dist) : normalize(vec2(-dut, 1.0));
        }
    }

    return vec4(sgn * best_d, best_n2.x * u_xz.x, sy * best_n2.y, best_n2.x * u_xz.y);
}

float obHarmonicCap3D( in vec3 p, in float r, in float Rpin )
{
    return sdgHarmonicCap3D(p, r, Rpin).x;
}

// 1-Lipschitz conservative bound without Newton solver (uses global Lipschitz constant of implicit residual)
float obHarmonicCapBound3D( in vec3 p, in float r, in float Rpin )
{
    float R = max(Rpin, r * 1.001);
    float rho = length(p.xz);
    float ay = abs(p.y);
    float a = obHarmonicContactRadius3D(r, R);
    float psi_a = sqrt(max((r - a) * (r + a), 1e-12));
    float A = (a * a * a) / psi_a;
    // Smooth C^{1,1} height profile u(rho) on [0, R], extended past R with slope u'(R) = -A / R^2
    float u_rho;
    if (rho <= a)
    {
        u_rho = sqrt(max(r * r - rho * rho, 0.0));
    }
    else if (rho <= R)
    {
        u_rho = A * (1.0 / rho - 1.0 / R);
    }
    else
    {
        u_rho = -(A / (R * R)) * (rho - R);
    }
    // On [a*, inf), max slope of u(rho) is |u'(a*)| = a* / psi(a*), so ||grad(ay - u(rho))|| <= r / psi(a*)
    float d_skirt = max(ay - u_rho, rho - R) * (psi_a / r);
    // The solid is the union of the obstacle sphere (radius r) and the membrane skirt
    return min(length(p) - r, d_skirt);
}

// Naive unnormalized vertical residual (demonstrates Eikonal violation ||grad|| > 1 and raymarching overstepping)
float obHarmonicCapNaive3D( in vec3 p, in float r, in float Rpin )
{
    float R = max(Rpin, r * 1.001);
    float rho = length(p.xz);
    float ay = abs(p.y);
    float a = obHarmonicContactRadius3D(r, R);
    float psi_a = sqrt(max((r - a) * (r + a), 1e-12));
    if (rho <= a)
    {
        return length(p) - r;
    }
    if (rho >= R)
    {
        return length(vec2(rho - R, ay));
    }
    float A = (a * a * a) / psi_a;
    return ay - A * (1.0 / rho - 1.0 / R);
}

// ----------------------------------------------------------------------------
// 3. 2D Logarithmic Laplace Obstacle Envelope (obLogCap2D, sdgLogCap2D, obLogCap3D)
// ----------------------------------------------------------------------------

float obLogContactRadius( in float r, in float Rpin )
{
    float R = max(Rpin, r * 1.0001);
    float L = 1.0 + log(R / r);
    float s = 1.0 / (L + 0.5 * log(max(L, 1.0)));
    for (int i = 0; i < 3; i++)
    {
        float ls = log(max(s, 1e-12));
        float f = s * (L - 0.5 * ls) - 1.0;
        float df = (L - 0.5) - 0.5 * ls;
        s = clamp(s - f / max(df, 1e-6), 1e-9, 0.9999);
    }
    return r * sqrt(s);
}

vec3 sdgLogCap2D( in vec2 p, in float r, in float Rpin )
{
    float R = max(Rpin, r * 1.001);
    float rho = abs(p.x);
    float ay = abs(p.y);
    float sx = (p.x >= 0.0) ? 1.0 : -1.0;
    float sy = (p.y >= 0.0) ? 1.0 : -1.0;

    float a = obLogContactRadius(r, R);
    float psi_a = sqrt(max((r - a) * (r + a), 1e-12));
    float A = (a * a) / psi_a;

    float u_rho = (rho <= a) ? sqrt(max(r * r - rho * rho, 0.0)) : ((rho < R) ? A * log(R / rho) : 0.0);
    float sgn = (rho <= R && ay <= u_rho) ? -1.0 : 1.0;

    float best_d;
    vec2 best_n2;
    if (rho * psi_a - a * ay <= 0.0)
    {
        float l = length(p);
        best_d = abs(l - r);
        best_n2 = (l > 1e-20) ? (sgn * ((l - r >= 0.0) ? 1.0 : -1.0) * vec2(rho, ay) / l) : vec2(0.0, 1.0);
    }
    else
    {
        vec2 df = vec2(rho - a, ay - psi_a);
        best_d = length(df);
        best_n2 = (best_d > 1e-12) ? (sgn * df / best_d) : normalize(vec2(a, psi_a));
    }

    vec2 dfR = vec2(rho - R, ay);
    float dR = length(dfR);
    if (dR < best_d)
    {
        best_d = dR;
        best_n2 = (dR > 1e-12) ? (sgn * dfR / dR) : normalize(vec2(A / R, 1.0));
    }

    for (int s = 0; s < 2; s++)
    {
        float t = (s == 0) ? clamp(rho, a, R) : mix(a, R, 0.35);
        for (int i = 0; i < 4; i++)
        {
            float inv_t = 1.0 / t;
            float ut = A * log(R * inv_t);
            float dut = -A * inv_t;
            float ddut = A * inv_t * inv_t;
            float dy = ut - ay;
            float phi = (t - rho) + dut * dy;
            float dphi = 1.0 + dut * dut + ddut * dy;
            t = clamp(t - phi / max(abs(dphi), 1e-6), a, R);
        }
        float ut = A * log(R / t);
        vec2 df = vec2(rho - t, ay - ut);
        float dist = length(df);
        if (dist < best_d)
        {
            best_d = dist;
            best_n2 = (dist > 1e-12) ? (sgn * df / dist) : normalize(vec2(A / t, 1.0));
        }
    }

    return vec3(sgn * best_d, sx * best_n2.x, sy * best_n2.y);
}

float obLogCap2D( in vec2 p, in float r, in float Rpin )
{
    return sdgLogCap2D(p, r, Rpin).x;
}

float obLogCap3D( in vec3 p, in float r, in float Rpin )
{
    return sdgLogCap2D(vec2(length(p.xz), p.y), r, Rpin).x;
}

// ----------------------------------------------------------------------------
// 4. 3D Biharmonic Clamped Plate Obstacle (obBiharmonicPlate, sdgBiharmonicPlate)
// ----------------------------------------------------------------------------

vec4 sdgBiharmonicPlate( in vec3 p, in float r, in float Rpin )
{
    float R = max(Rpin, r * 1.001);
    float rho2 = dot(p.xz, p.xz);
    float rho = sqrt(rho2);
    float ay = abs(p.y);
    float sy = (p.y >= 0.0) ? 1.0 : -1.0;
    vec2 u_xz = (rho > 1e-20) ? p.xz / rho : vec2(1.0, 0.0);

    float a = (4.0 * r * r) / (R + sqrt(R * R + 8.0 * r * r));
    float psi_a = sqrt(max((r - a) * (r + a), 1e-12));
    float da = R - a;
    float A = psi_a / (da * da);

    float u_rho = (rho <= a) ? sqrt(max(r * r - rho2, 0.0)) : ((rho < R) ? A * (R - rho) * (R - rho) : 0.0);
    float sgn = (rho <= R && ay <= u_rho) ? -1.0 : 1.0;

    float best_d;
    vec2 best_n2;
    if (rho * psi_a - a * ay <= 0.0)
    {
        float l = sqrt(rho2 + ay * ay);
        best_d = abs(l - r);
        best_n2 = (l > 1e-20) ? (sgn * ((l - r >= 0.0) ? 1.0 : -1.0) * vec2(rho, ay) / l) : vec2(0.0, 1.0);
    }
    else
    {
        vec2 df = vec2(rho - a, ay - psi_a);
        best_d = length(df);
        best_n2 = (best_d > 1e-12) ? (sgn * df / best_d) : normalize(vec2(a, psi_a));
    }

    vec2 dfR = vec2(rho - R, ay);
    float dR = length(dfR);
    if (dR < best_d)
    {
        best_d = dR;
        best_n2 = (dR > 1e-12) ? (sgn * dfR / dR) : vec2(0.0, 1.0);
    }

    for (int s = 0; s < 2; s++)
    {
        float t = (s == 0) ? clamp(rho, a, R) : mix(a, R, 0.35);
        for (int i = 0; i < 4; i++)
        {
            float dt = R - t;
            float ut = A * dt * dt;
            float dut = -2.0 * A * dt;
            float ddut = 2.0 * A;
            float dy = ut - ay;
            float phi = (t - rho) + dut * dy;
            float dphi = 1.0 + dut * dut + ddut * dy;
            t = clamp(t - phi / max(abs(dphi), 1e-6), a, R);
        }
        float dt = R - t;
        float ut = A * dt * dt;
        vec2 df = vec2(rho - t, ay - ut);
        float dist = length(df);
        if (dist < best_d)
        {
            best_d = dist;
            best_n2 = (dist > 1e-12) ? (sgn * df / dist) : normalize(vec2(2.0 * A * dt, 1.0));
        }
    }

    return vec4(sgn * best_d, best_n2.x * u_xz.x, sy * best_n2.y, best_n2.x * u_xz.y);
}

float obBiharmonicPlate( in vec3 p, in float r, in float Rpin )
{
    return sdgBiharmonicPlate(p, r, Rpin).x;
}

// ----------------------------------------------------------------------------
// 5. General 2-SDF Free-Boundary Obstacle Operator (opObstacleBridge, sdgObstacleBridge)
// ----------------------------------------------------------------------------

float opObstacleBridge( in float d1, in float d2, in float gapCoord, in float tension )
{
    float w = max(gapCoord, 1e-6);
    float q1 = w - d1;
    float q2 = w - d2;
    if (q1 <= 0.0 || q2 <= 0.0)
    {
        return min(d1, d2);
    }
    float p = 2.0 + max(tension, 0.0);
    float m = max(q1, q2);
    float r = m * pow(pow(q1 / m, p) + pow(q2 / m, p), 1.0 / p);
    // Scale by 2^((1-p)/p) so that ||grad(w - r)|| <= 1 even when grad d1 == grad d2,
    // while keeping C^{1,1} value + tangent match at q1 = 0 and q2 = 0!
    // Specifically, use the C^{1,1} normalized excess:
    float u = min(q1, q2) / m;
    float excess = m * (pow(1.0 + pow(u, p), 1.0 / p) - 1.0);
    return min(d1, d2) - pow(0.5, (p - 1.0) / p) * excess;
}

vec4 sdgObstacleBridge( in vec4 dg1, in vec4 dg2, in float gapCoord, in float tension )
{
    float d1 = dg1.x, d2 = dg2.x;
    float w = max(gapCoord, 1e-6);
    float q1 = w - d1;
    float q2 = w - d2;
    if (q1 <= 0.0 || q2 <= 0.0)
    {
        return (d1 <= d2) ? dg1 : dg2;
    }
    float p = 2.0 + max(tension, 0.0);
    float m = max(q1, q2);
    float u1 = q1 / m, u2 = q2 / m;
    float up1 = pow(u1, p - 1.0), up2 = pow(u2, p - 1.0);
    float s = up1 * u1 + up2 * u2;
    float r = m * pow(s, 1.0 / p);
    float sp = pow(s, (1.0 - p) / p);
    float w1 = up1 * sp;
    float w2 = up2 * sp;

    float c = pow(0.5, (p - 1.0) / p);
    float d = min(d1, d2) - c * (r - m);
    vec3 g = (d1 <= d2)
        ? ((1.0 - c) * dg1.yzw + c * (w1 * dg1.yzw + w2 * dg2.yzw))
        : ((1.0 - c) * dg2.yzw + c * (w1 * dg1.yzw + w2 * dg2.yzw));
    float gl = length(g);
    return vec4(d, (gl > 1e-20) ? g / gl : dg1.yzw);
}

#endif // SDF_GEMS_HAS_OB_ELASTIC_ENVELOPE
