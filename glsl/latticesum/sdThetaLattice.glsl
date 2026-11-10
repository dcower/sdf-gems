// ============================================================================
// sdThetaLattice — Poisson-summed (Jacobi theta) smooth union of an ENTIRE lattice of spheres
//
//   S(x)  = sum_{n in Lattice} exp(-|x - n|^2 / (2 s^2))        (Gaussian lattice sum, all images)
//   E(x)  = -2 s^2 log( S(x) / S(0) )   >= 0                    (log-sum-exp soft squared distance, re-zeroed)
//   L(s)^2 = (1/2) d^2E/dx^2 (0) = -s^2 (d^2/dx^2 log S)(0)     (closed form in s)
//   f(x)  = sqrt(E(x)) / L(s) - r
//
// For the simple cubic lattice S = theta(x) theta(y) theta(z), so E / L^2 is a SUM of three 1D terms.
// theta is evaluated by its image series (s < 0.22) or its Fourier series (s >= 0.22), both written in
// cancellation-free, underflow-free "soft squared distance" form.
//
//   s -> 0   : f -> |x - round(x)| - r        (the exact sphere lattice SDF)
//   s -> inf : f -> sqrt(sin^2(pi x) + sin^2(pi y) + sin^2(pi z)) / pi - r   (Schwarz-P nodal family)
//
// Lattice period is 1 (conventional cubic cell); scale with  c * sdThetaLatticeSC(p / c, s, r / c).
// The surface is non-empty only for r < max sqrt(E)/L (SC: between sqrt(3)/pi = 0.551 and sqrt(3)/2).
//
// @author Claude Fable 5.2 (sdf-gems, 2026)
// @exact  Conservative bound, NOT an exact SDF. Simple cubic: globally 1-Lipschitz (|grad f| <= 1, equality
//         only at the lattice points; lemma verified numerically + proof sketch in the doc). BCC / FCC:
//         1-Lipschitz verified numerically only (no proof). C-infinity away from the lattice points.
// ============================================================================

#ifndef SDF_GEMS_SD_THETA_LATTICE_GLSL
#define SDF_GEMS_SD_THETA_LATTICE_GLSL

#ifndef SDF_GEMS_HAS_THETA_LATTICE_HELPERS
#define SDF_GEMS_HAS_THETA_LATTICE_HELPERS

// Crossover between the image series (s < SC) and the Fourier series (s >= SC).
#define SDF_GEMS_THETA_LATTICE_SC 0.22

// log(1 + x), accurate for tiny x.
float thetaLatticeLog1p( float x )
{
    return (abs(x) < 2e-2) ? x * (1.0 - x * (0.5 - x * (1.0 / 3.0 - x * (0.25 - 0.2 * x)))) : log(1.0 + x);
}

// -log(1 - x) / x  for 0 <= x < 1, accurate for tiny x (-> 1).
float thetaLatticeLogRatio( float x )
{
    return (x < 2e-2) ? 1.0 + x * (0.5 + x * (1.0 / 3.0 + x * (0.25 + 0.2 * x))) : -log(1.0 - x) / x;
}

// Image regime. u must be reduced to [-1/2, 1/2].
// Returns vec2( e(u), e'(u) ),  e(u) = -2 s^2 log(theta(u) / theta(0))  using the images n = -1, 0, 1.
// Never underflows: e = u^2 - 2 s^2 log(1 + (E+ - E-)^2 / (1 + 2 e1)).
vec2 thetaLatticeImageE( float u, float s )
{
    float i4 = 0.25 / (s * s);
    float ep = exp((2.0 * u - 1.0) * i4);          // sqrt of the n = +1 image weight relative to n = 0
    float em = exp((-2.0 * u - 1.0) * i4);         // sqrt of the n = -1 image weight
    float e1 = exp(-2.0 * i4);                     // exp(-1 / (2 s^2)) = ep * em
    float df = ep - em;
    float m  = df * (ep + em) / (1.0 + ep * ep + em * em);   // Gibbs mean lattice point E[n]
    return vec2(u * u - 2.0 * s * s * thetaLatticeLog1p(df * df / (1.0 + 2.0 * e1)), 2.0 * (u - m));
}

// Image regime: L^2 = e''(0) / 2 = 1 - Var_0[n] / s^2 for the simple cubic lattice.
float thetaLatticeImageL2( float s )
{
    float e1 = exp(-0.5 / (s * s));
    return 1.0 - 2.0 * e1 / (s * s * (1.0 + 2.0 * e1));
}

#endif // SDF_GEMS_HAS_THETA_LATTICE_HELPERS

// ----------------------------------------------------------------------------
// 1. Plain 1D lattice sum of Gaussians (Jacobi theta-3):
//    vec2( theta(u), dtheta/du ),  theta(u) = sum_n exp(-(u - n)^2 / (2 s^2)).
//    NOTE: for small s the plain value underflows / loses relative accuracy in FP32 between lattice points
//    (it is ~exp(-1 / (8 s^2))). The distance functions below use the log-domain form instead.
// ----------------------------------------------------------------------------
vec2 thetaLattice1D( float u, float s )
{
    s = max(s, 1e-3);
    u -= floor(u + 0.5);
    if( s < SDF_GEMS_THETA_LATTICE_SC )
    {
        float i2 = 0.5 / (s * s);
        float w0 = exp(-u * u * i2);
        float wp = exp(-(u - 1.0) * (u - 1.0) * i2);
        float wm = exp(-(u + 1.0) * (u + 1.0) * i2);
        return vec2(w0 + wp + wm, -2.0 * i2 * (u * w0 + (u - 1.0) * wp + (u + 1.0) * wm));
    }
    float q  = exp(-19.7392088 * s * s);           // exp(-2 pi^2 s^2)
    float q2 = q * q, q4 = q2 * q2, q9 = q4 * q4 * q, q16 = q4 * q4 * q4 * q4;
    float c1 = cos(6.28318531 * u), s1 = sin(6.28318531 * u);
    float c2 = 2.0 * c1 * c1 - 1.0, s2 = 2.0 * s1 * c1;
    float c3 = c1 * (2.0 * c2 - 1.0), s3 = s1 * (2.0 * c2 + 1.0);
    float c4 = 2.0 * c2 * c2 - 1.0, s4 = 2.0 * s2 * c2;
    float a  = 2.50662827 * s;                     // s sqrt(2 pi)
    return vec2(a * (1.0 + 2.0 * (q * c1 + q4 * c2 + q9 * c3 + q16 * c4)),
                -12.5663706 * a * (q * s1 + 2.0 * q4 * s2 + 3.0 * q9 * s3 + 4.0 * q16 * s4));
}

// ----------------------------------------------------------------------------
// 2. Normalised 1D soft squared distance to the integer lattice:
//    vec2( h2(u), dh2/du ),  h2(u) = -2 s^2 log(theta(u) / theta(0)) / L(s)^2,   h2 ~ u^2 near u = 0.
//    sqrt(h2) is 1-Lipschitz for every s;  s -> 0: h2 -> (u - round(u))^2;  s -> inf: h2 -> sin^2(pi u) / pi^2.
//    Finite for every s (no underflow at small s, no cancellation at large s).
// ----------------------------------------------------------------------------
vec2 thetaLatticeSoftSq1D( float u, float s )
{
    s = max(s, 1e-3);
    u -= floor(u + 0.5);
    if( s < SDF_GEMS_THETA_LATTICE_SC )
        return thetaLatticeImageE(u, s) / thetaLatticeImageL2(s);

    float q   = exp(-19.7392088 * s * s);          // exp(-2 pi^2 s^2)
    float q3  = q * q * q, q8 = q3 * q3 * q * q, q15 = q8 * q3 * q3 * q;
    float a   = sin(3.14159265 * u), b = cos(3.14159265 * u);
    float S1  = 2.0 * a * b, C1 = 1.0 - 2.0 * a * a;           // sin, cos of 2 pi u
    float S2  = 2.0 * S1 * C1, C2 = 2.0 * C1 * C1 - 1.0;       // sin, cos of 4 pi u
    float t3  = a * (3.0 - 4.0 * a * a);                       // sin(3 pi u)
    float S3  = S1 * (4.0 * C1 * C1 - 1.0), S4 = 2.0 * S2 * C2;
    // theta(u)/theta(0) = 1 - dl,  dl = 4 sum_k q^(k^2) sin^2(pi k u) / (1 + 2 sum_k q^(k^2))
    float sum = a * a + q3 * S1 * S1 + q8 * t3 * t3 + q15 * S2 * S2;
    float w   = 1.0 + 2.0 * q * (1.0 + q3 + q8 + q15);
    float dl  = 4.0 * q * sum / w;
    float M   = 1.0 + 4.0 * q3 + 9.0 * q8 + 16.0 * q15;
    return vec2(sum * thetaLatticeLogRatio(dl) / (9.8696044 * M),
                (S1 + 2.0 * q3 * S2 + 3.0 * q8 * S3 + 4.0 * q15 * S4) / (3.14159265 * M * (1.0 - dl)));
}

// ----------------------------------------------------------------------------
// 3. Simple cubic lattice (points Z^3): seamless smooth union of all spheres of radius r.
//    s = Gaussian softness in lattice units (0 = hard union, >~ 0.5 = Schwarz-P nodal limit).
//    The vec3 overload takes one softness per axis (still 1-Lipschitz, by the same 1D lemma).
// ----------------------------------------------------------------------------
float sdThetaLatticeSC( in vec3 p, in vec3 s, in float r )
{
    return sqrt(thetaLatticeSoftSq1D(p.x, s.x).x + thetaLatticeSoftSq1D(p.y, s.y).x + thetaLatticeSoftSq1D(p.z, s.z).x) - r;
}

float sdThetaLatticeSC( in vec3 p, in float s, in float r )
{
    return sdThetaLatticeSC(p, vec3(s), r);
}

// vec4( f, grad f ). The gradient is the TRUE gradient of the field (|grad f| <= 1, not normalised):
// normalise it for shading. At a lattice point (a cone apex of f) a unit fallback is returned.
vec4 sdgThetaLatticeSC( in vec3 p, in float s, in float r )
{
    vec2  hx = thetaLatticeSoftSq1D(p.x, s), hy = thetaLatticeSoftSq1D(p.y, s), hz = thetaLatticeSoftSq1D(p.z, s);
    float h2 = hx.x + hy.x + hz.x;
    float h  = sqrt(h2);
    vec3  g  = (h2 > 0.0) ? vec3(hx.y, hy.y, hz.y) * (0.5 / h) : vec3(1.0, 0.0, 0.0);
    return vec4(h - r, g);
}

// ----------------------------------------------------------------------------
// 4. BCC and FCC lattices as 2 and 4 cubic cosets (conventional cubic cell of side 1).
//    BCC: Z^3 + {(0,0,0), (1/2,1/2,1/2)}                nearest-neighbour distance sqrt(3)/2
//    FCC: Z^3 + {(0,0,0), (0,1/2,1/2), (1/2,0,1/2), (1/2,1/2,0)}   nearest-neighbour distance 1/sqrt(2)
//    Image regime: log-sum-exp over the cosets of sums of 1D soft squared distances.
//    Fourier regime: even/odd split theta(u) = Ce + Co, theta(u - 1/2) = Ce - Co, which gives
//      S_BCC = 2 (Ce Ce Ce + Ce Co Co + Co Ce Co + Co Co Ce),   S_FCC = 4 (Ce Ce Ce + Co Co Co),
//    evaluated as cancellation-free differences from the value at the lattice point.
//    type: 0 = BCC, 1 = FCC.
// ----------------------------------------------------------------------------
float thetaLatticeCosetSoftSq( in vec3 p, in float s, in int type )
{
    s = max(s, 1e-3);
    vec3 u = p - floor(p + 0.5);
    // Move the nearest lattice point (of any coset) to the origin. The field is invariant under the coset
    // shifts, and every formula below is cancellation-free only near the origin coset.
    {
        vec3 h = u - 0.5 * sign(u) + 0.5 * (1.0 - abs(sign(u)));      // u - 1/2 reduced to [-1/2, 1/2]
        if( type == 0 )
        {
            if( dot(h, h) < dot(u, u) ) u = h;
        }
        else
        {
            vec3 uu = u * u, hh = h * h;
            vec4 n2 = vec4(uu.x + uu.y + uu.z, uu.x + hh.y + hh.z, hh.x + uu.y + hh.z, hh.x + hh.y + uu.z);
            float m = min(min(n2.x, n2.y), min(n2.z, n2.w));
            u = (m == n2.x) ? u : (m == n2.y) ? vec3(u.x, h.y, h.z) : (m == n2.z) ? vec3(h.x, u.y, h.z) : vec3(h.x, h.y, u.z);
        }
    }
    if( s < SDF_GEMS_THETA_LATTICE_SC )
    {
        vec3 v  = u - 0.5 * sign(u) + 0.5 * (1.0 - abs(sign(u)));     // u - 1/2 reduced to [-1/2, 1/2]
        vec3 ea = vec3(thetaLatticeImageE(u.x, s).x, thetaLatticeImageE(u.y, s).x, thetaLatticeImageE(u.z, s).x);
        vec3 eb = vec3(thetaLatticeImageE(v.x, s).x, thetaLatticeImageE(v.y, s).x, thetaLatticeImageE(v.z, s).x);
        float eh  = thetaLatticeImageE(0.5, s).x;                     // e(1/2) = -2 s^2 log(theta(1/2) / theta(0))
        float i2  = 0.5 / (s * s);
        float e1  = exp(-i2);
        float La  = 1.0 - 2.0 * e1 / (s * s * (1.0 + 2.0 * e1));      // -s^2 theta''(0) / theta(0)
        float Lb  = 1.0 - (0.25 + 2.0 * e1 * e1 / (1.0 + e1 * e1)) / (s * s);   // -s^2 theta''(1/2) / theta(1/2)
        float T0  = ea.x + ea.y + ea.z;
        if( type == 0 )
        {
            float T1  = eb.x + eb.y + eb.z;
            float Tm  = min(T0, T1);
            float num = exp(-(T0 - Tm) * i2) + exp(-(T1 - Tm) * i2);
            float rho = exp(-3.0 * eh * i2);
            float E   = Tm - 2.0 * s * s * log(num / (1.0 + rho));
            return max(E, 0.0) * (1.0 + rho) / (La + Lb * rho);
        }
        vec3  T   = vec3(ea.x + eb.y + eb.z, eb.x + ea.y + eb.z, eb.x + eb.y + ea.z);
        float Tm  = min(T0, min(T.x, min(T.y, T.z)));
        float num = exp(-(T0 - Tm) * i2) + dot(exp(-(T - Tm) * i2), vec3(1.0));
        float rho = exp(-2.0 * eh * i2);
        float E   = Tm - 2.0 * s * s * log(num / (1.0 + 3.0 * rho));
        return max(E, 0.0) * (1.0 + 3.0 * rho) / (La * (1.0 + rho) + 2.0 * Lb * rho);
    }

    float q   = exp(-19.7392088 * s * s);          // exp(-2 pi^2 s^2)
    float q2  = q * q, q4 = q2 * q2, q8 = q4 * q4, q12 = q8 * q4, q16 = q8 * q8;
    vec3  a   = sin(3.14159265 * u), b = cos(3.14159265 * u);
    vec3  c1  = 1.0 - 2.0 * a * a;                 // cos(2 pi u)
    vec3  t2  = 2.0 * a * b;                       // sin(2 pi u)
    vec3  t3  = a * (3.0 - 4.0 * a * a);           // sin(3 pi u)
    vec3  t4  = 2.0 * t2 * c1;                     // sin(4 pi u)
    float we  = 1.0 + 2.0 * q4 + 2.0 * q16;        // even part of theta(0) / (s sqrt(2 pi))
    float wo1 = 1.0 + q8;                          // odd part / (2 q)
    float wo  = 2.0 * q * wo1;
    float Me  = 1.0 + 9.0 * q8, mo1 = 1.0 + 4.0 * q12;
    vec3  Ce  = we - 4.0 * q4 * (t2 * t2 + q12 * t4 * t4);            // even part of theta(u)
    vec3  Cq  = 2.0 * c1 * (1.0 + q8 * (4.0 * c1 * c1 - 3.0));        // odd part / q
    vec3  dq  = 4.0 * (a * a + q8 * t3 * t3);                         // (wo - Co) / q   >= 0
    float wq  = 2.0 * wo1;                                            // wo / q
    if( type == 0 )
    {
        vec3  d2  = 4.0 * q2 * (t2 * t2 + q12 * t4 * t4);            // (we - Ce) / q^2 >= 0
        float nE  = d2.x * we * we + Ce.x * (d2.y * we + Ce.y * d2.z);          // (we^3 - Ce Ce Ce) / q^2
        // sum over i of (we wo^2 - Ce_i Co_j Co_k) / q^2
        float nO  = d2.x * q2 * Cq.y * Cq.z + we * (dq.y * wq + Cq.y * dq.z)
                  + d2.y * q2 * Cq.z * Cq.x + we * (dq.z * wq + Cq.z * dq.x)
                  + d2.z * q2 * Cq.x * Cq.y + we * (dq.x * wq + Cq.x * dq.y);
        float N   = nE + nO;
        float S0  = we * we * we + 3.0 * we * wo * wo;
        return N * thetaLatticeLogRatio(min(N * q2 / S0, 0.999999)) / (157.913670 * (Me * we * wo1 + q2 * mo1 * (we * we + wo * wo)));
    }
    vec3  d3  = 4.0 * q * (t2 * t2 + q12 * t4 * t4);                 // (we - Ce) / q^3 >= 0
    float nE  = d3.x * we * we + Ce.x * (d3.y * we + Ce.y * d3.z);              // (we^3 - Ce Ce Ce) / q^3
    float nO  = dq.x * wq * wq + Cq.x * (dq.y * wq + Cq.y * dq.z);              // (wo^3 - Co Co Co) / q^3
    float N   = nE + nO;
    float S0  = we * we * we + wo * wo * wo;
    return N * thetaLatticeLogRatio(min(N * q2 * q / S0, 0.999999)) / (157.913670 * (Me * wo1 * wo1 + q * mo1 * we * we));
}

float sdThetaLatticeBCC( in vec3 p, in float s, in float r )
{
    return sqrt(max(thetaLatticeCosetSoftSq(p, s, 0), 0.0)) - r;
}

float sdThetaLatticeFCC( in vec3 p, in float s, in float r )
{
    return sqrt(max(thetaLatticeCosetSoftSq(p, s, 1), 0.0)) - r;
}

#endif // SDF_GEMS_SD_THETA_LATTICE_GLSL
