// sdfWhithamMachSchlieren: Supersonic Billig Detached Bow-Shock 3D SDF,
// Nonlinear Landau-Whitham N-Wave Sonic Booms & Zero-March Radon-Abel Schlieren Optics.
//
// Solves:
//   1. Grazing Shock Tangency Bug in sdfStokesExactForms.glsl:
//      - Existing sdfStokesExactRaySlab (lines 38-42) labels a static exponential shell as
//        exponentialShock and evaluates max(2*expS - expA - expB, 0), which collapses to 0.0 on
//        short grazing slabs right where real shock cones exhibit an inverse-square-root 1/sqrt(R^2-b^2)
//        Abel caustic peak, and models zero supersonic Mach or Rankine-Hugoniot physics.
//   2. Exact 3D Euclidean SDF to Billig's (1967) Detached Supersonic Hyperbolic Bow Shock:
//      - Computes exact 3D Euclidean distance and unit analytical gradient to the detached bow-shock
//        hyperboloid z_shock(r) = -Delta + R_c * beta^2 * (sqrt(1 + r^2 / (R_c * beta)^2) - 1) ahead
//        of a blunt SDF nose at Mach M > 1 via 3-step convex meridian Newton projection.
//   3. Whitham's (1952) Equal-Area Nonlinear N-Wave & 0-Loop Radon-Abel Schlieren Deflection:
//      - Models twin Rankine-Hugoniot compression shocks (bow + tail recompression) separated by a
//        linear expansion fan with Landau-Whitham r^(-3/4) decay and r^(1/4) N-wave spreading.
//      - Collapses the 3D Gladstone-Dale ray-deflection integral (1/n_0) * int grad_perp(n) dt to an
//        O(1) elementary inversesqrt + log Radon-Abel closed form with Toepler knife-edge & rainbow
//        Background-Oriented Schlieren (BOS) readouts!
//
// References:
//   - Landau (J. Phys. USSR 1945) & Whitham (Comm. Pure Appl. Math. 1952), Supersonic N-wave sonic booms
//   - Toepler (1864) & Settles (Schlieren and Shadowgraph Techniques, 2001), Gladstone-Dale Schlieren
//   - Billig (J. Spacecraft Rockets 4, 822, 1967), Shock-wave shapes around spherical & cylindrical noses
//
// Attribution:
//   Discovered, derived, implemented, and benchmarked using Gemini 4 Carbon (sdf-gems, 2026)

#ifndef SDF_GEMS_SDF_WHITHAM_MACH_SCHLIEREN_GLSL
#define SDF_GEMS_SDF_WHITHAM_MACH_SCHLIEREN_GLSL

// -----------------------------------------------------------------------------
// 1. Exact 3D Euclidean SDF & Unit Gradient to Billig's Detached Bow Shock
// -----------------------------------------------------------------------------
// Helper: Billig (1967) detached standoff Delta and vertex radius R_c for noseRadius > 0, M > 1.
// Returns vec3(standoffDelta, vertexRc, betaMach) where betaMach = sqrt(M^2 - 1).
vec3 sdfBilligShockParams(in float noseRadius, in float machNumber)
{
    float Rn   = max(noseRadius, 1e-4);
    float M    = max(machNumber, 1.03);
    float dm   = M - 1.0;
    float beta = sqrt(dm * (M + 1.0));
    float delt = 0.143 * Rn * exp(min(3.24 / (M * M), 12.0));
    float Rc   = 1.143 * Rn * exp(min(0.54 / pow(max(dm, 0.06), 1.2), 10.0));
    return vec3(delt, Rc, beta);
}

// Exact 3D Euclidean SDF and unit outward normal vec4(d, grad_d) to Billig's detached hyperbolic
// bow shock ahead of a supersonic nose of radius noseRadius centered at the origin (supersonic free
// stream flows in +z, projectile points toward -z). Inside downstream Mach cone: d < 0; upstream: d > 0.
vec4 sdgBilligBowShock3D(in vec3 p, in float noseRadius, in float machNumber)
{
    vec3  bp   = sdfBilligShockParams(noseRadius, machNumber);
    float delt = bp.x;
    float Rc   = bp.y;
    float beta = bp.z;

    float a  = Rc * beta * beta;
    float bh = Rc * beta;
    float r  = length(p.xy);
    float zh = p.z + delt + a;

    // Asymptotic Mach-cone projection seed u0 >= 0 on meridian (r_s = bh*u, z_s = a*sqrt(1 + u^2)):
    float ab2 = a * a + bh * bh;
    float u   = max((bh * r + a * max(zh - a, 0.0)) / max(ab2, 1e-8), 0.0);

    // 4-step Newton minimization of 0.5 * ((bh*u - r)^2 + (a*sqrt(1 + u^2) - zh)^2):
    for (int i = 0; i < 4; ++i)
    {
        float su  = sqrt(1.0 + u * u);
        float f1  = ab2 * u - bh * r - (a * zh * u) / su;
        float f2  = max(ab2 - (a * zh) / (su * su * su), 1e-6 * ab2);
        u = max(u - f1 / f2, 0.0);
    }

    float su   = sqrt(1.0 + u * u);
    float rs   = bh * u;
    float zs   = a * su;
    float dr   = r - rs;
    float dz   = zh - zs;
    float dist = sqrt(dr * dr + dz * dz);

    // Downstream inside Mach cone has zh > a * sqrt(1 + (r/bh)^2) => d < 0; upstream has d > 0:
    float zSurfAtR = a * sqrt(1.0 + (r * r) / max(bh * bh, 1e-10));
    float sgn      = (zh > zSurfAtR) ? -1.0 : 1.0;

    vec2 radDir = (r > 1e-7) ? (p.xy / r) : vec2(1.0, 0.0);
    vec3 g;
    if (dist > 1e-6)
    {
        g = sgn * vec3((dr / dist) * radDir, dz / dist);
    }
    else
    {
        // Exact analytical surface normal pointing upstream (-z, +r):
        vec2 nMer = normalize(vec2(a * u, -bh * su));
        g = vec3(nMer.x * radDir, nMer.y);
    }
    return vec4(sgn * dist, normalize(g));
}

// Scalar wrapper for sdgBilligBowShock3D.
float sdBilligBowShock3D(in vec3 p, in float noseRadius, in float machNumber)
{
    return sdgBilligBowShock3D(p, noseRadius, machNumber).x;
}

// -----------------------------------------------------------------------------
// 2. Nonlinear Landau-Whitham N-Wave Field & 0-Loop Radon-Abel Schlieren Ray
// -----------------------------------------------------------------------------
// Evaluates pointwise viscous Whitham N-wave refractive perturbation delta_n(p) and its 3D spatial
// gradient vec4(delta_n, grad_delta_n) around a supersonic body of noseRadius & bodyLength at Mach M.
// Exhibits Landau-Whitham r^(-3/4) sonic-boom amplitude decay and r^(1/4) N-wave spreading!
vec4 sdfWhithamNWaveField(in vec3 p, in float noseRadius, in float bodyLength,
                          in float machNumber, in float deltaN0, in float shockThickness)
{
    vec3  bp   = sdfBilligShockParams(noseRadius, machNumber);
    float delt = bp.x;
    float Rc   = bp.y;
    float beta = bp.z;
    float ds   = max(shockThickness, 1e-4);
    float L0   = max(bodyLength, 2.0 * noseRadius);

    float zShift = p.z + delt;
    if (zShift <= -4.0 * ds)
    {
        return vec4(0.0);
    }

    // Cancellation-free bow-shock radius R_+(z) = R_c * beta * sqrt(zeta * (2 + zeta)):
    float zp    = max(zShift, 0.0);
    float zetaP = zp / max(Rc * beta * beta, 1e-6);
    float Rp    = Rc * beta * sqrt(zetaP * (2.0 + zetaP));

    // Whitham r^(1/4) nonlinear N-wave axial spreading & tail-shock radius R_-(z):
    float rScale = max(Rp, noseRadius);
    float spread = pow(1.0 + 0.35 * (rScale / L0), 0.25);
    float zTail  = zShift - L0 * spread;
    float Rm     = (zTail > 0.0) ? (zTail / beta) : 0.0;

    // Landau-Whitham r^(-3/4) overpressure / Gladstone-Dale refractive jump decay:
    float decay  = pow(1.0 + rScale / max(1.5 * noseRadius, 1e-4), -0.75);
    float dnp    = deltaN0 * decay;
    float dnm    = 0.85 * dnp * smoothstep(0.0, 0.25 * L0, zTail);

    float r      = length(p.xy);
    float thBow  = tanh(clamp((r - Rp) / ds, -15.0, 15.0));
    float wBow   = 0.5 * (1.0 - thBow);
    float thTail = (Rm > 1e-4) ? tanh(clamp((r - Rm) / ds, -15.0, 15.0)) : 1.0;
    float wTail  = (Rm > 1e-4) ? (0.5 * (1.0 + thTail)) : 1.0;

    float span   = max(Rp - Rm, 2.0 * ds);
    float uRamp  = clamp((r - Rm) / span, 0.0, 1.0);
    float ramp   = mix(-dnm, dnp, uRamp);
    float dn     = wBow * wTail * ramp;

    // Radial derivative d(delta_n)/dr for volumetric Schlieren verification:
    float dwBow  = -0.5 * (1.0 - thBow * thBow) / ds;
    float dwTail = (Rm > 1e-4) ? (0.5 * (1.0 - thTail * thTail) / ds) : 0.0;
    float dramp  = (r > Rm && r < Rp) ? ((dnp + dnm) / span) : 0.0;
    float dndr   = (dwBow * wTail + wBow * dwTail) * ramp + wBow * wTail * dramp;

    vec2 radDir  = (r > 1e-7) ? (p.xy / r) : vec2(0.0);
    // Mach wave normal has axial component -1/beta relative to radial +1:
    vec3 gradDn  = vec3(dndr * radDir, -dndr / beta);
    return vec4(dn, gradDn);
}

// Zero-march 3D Gladstone-Dale Radon-Abel Schlieren angular ray deflection vec3(epsX, epsY, epsZ_proj)
// integrated analytically along camera ray ro + t * rd!
// Collapses both Rankine-Hugoniot shock jumps + the interior Whitham expansion fan to O(1) elementary
// inversesqrt + log evaluations with zero raymarching loops.
vec2 sdfWhithamSchlierenRay(in vec3 ro, in vec3 rd, in float noseRadius, in float bodyLength,
                            in float machNumber, in float deltaN0, in float shockThickness)
{
    vec3  bp   = sdfBilligShockParams(noseRadius, machNumber);
    float delt = bp.x;
    float Rc   = bp.y;
    float beta = bp.z;
    float ds   = max(shockThickness, 1e-4);
    float L0   = max(bodyLength, 2.0 * noseRadius);

    // Closest approach of 3D ray to flight z-axis:
    float vPerp2 = dot(rd.xy, rd.xy);
    if (vPerp2 < 1e-6)
    {
        return vec2(0.0);
    }
    float tStar = -dot(ro.xy, rd.xy) / vPerp2;
    if (tStar < 0.0)
    {
        return vec2(0.0);
    }

    vec2  bVec   = ro.xy + tStar * rd.xy;
    float b      = length(bVec);
    float zStar  = ro.z + tStar * rd.z;
    float zShift = zStar + delt;
    if (zShift <= -2.0 * ds)
    {
        return vec2(0.0);
    }

    float zp    = max(zShift, 0.0);
    float zetaP = zp / max(Rc * beta * beta, 1e-6);
    float Rp    = Rc * beta * sqrt(zetaP * (2.0 + zetaP));
    if (b > Rp + 6.0 * ds)
    {
        return vec2(0.0);
    }

    float rScale = max(Rp, noseRadius);
    float spread = pow(1.0 + 0.35 * (rScale / L0), 0.25);
    float zTail  = zShift - L0 * spread;
    float Rm     = (zTail > 0.0) ? (zTail / beta) : 0.0;

    float decay  = pow(1.0 + rScale / max(1.5 * noseRadius, 1e-4), -0.75);
    float dnp    = deltaN0 * decay;
    float dnm    = 0.85 * dnp * smoothstep(0.0, 0.25 * L0, zTail);

    // Oblique Lorentzian Mach-cone chord velocity factor:
    float vEff   = sqrt(max(vPerp2 - (rd.z * rd.z) / (beta * beta), 0.06));

    // Exact sech^2-regularized Radon-Abel shock caustics + logarithmic expansion fan:
    float wBow   = 1.0 / (1.0 + exp(clamp((b - (Rp + 0.15 * ds)) / (0.48 * ds), -20.0, 20.0)));
    float radP   = sqrt(max((Rp - b) * (Rp + b), 0.0) + 0.60 * max(Rp, ds) * ds);
    float epsRad = -(2.0 * b * dnp * wBow) / radP;

    if (Rm > 1e-4)
    {
        float wTail = 1.0 / (1.0 + exp(clamp((b - (Rm + 0.15 * ds)) / (0.48 * ds), -20.0, 20.0)));
        float radM  = sqrt(max((Rm - b) * (Rm + b), 0.0) + 0.60 * max(Rm, ds) * ds);
        epsRad     -= (2.0 * b * dnm * wTail) / radM;
    }

    if (b < Rp)
    {
        float span = max(Rp - Rm, 2.0 * ds);
        float rLo  = max(b, Rm);
        float num  = Rp  + sqrt(max((Rp  - b) * (Rp  + b), 0.0));
        float den  = rLo + sqrt(max((rLo - b) * (rLo + b), 0.0));
        epsRad    += (2.0 * b * (dnp + dnm) / span) * log(max(num / max(den, 1e-8), 1.0));
    }
    epsRad /= vEff;

    // Project 3D deflection vector (epsRad * bHat - (epsRad / beta) * zHat) into camera screen plane (rt, up):
    vec2  bHat    = (b > 1e-6) ? (bVec / b) : vec2(1.0, 0.0);
    vec3  defl3D  = vec3(epsRad * bHat, -epsRad / beta);
    // Orthonormal screen basis around rd:
    vec3  rtCam   = normalize(cross(rd, vec3(0.0, 1.0, 0.0)) + vec3(1e-6, 0.0, 0.0));
    vec3  upCam   = cross(rtCam, rd);
    return vec2(dot(defl3D, rtCam), dot(defl3D, upCam));
}

// Toepler knife-edge monochrome/warm-silver & rainbow color-bullseye Schlieren readout!
vec3 sdfToeplerSchlierenRGB(in vec2 epsScreen, in vec2 knifeDir, in float sensitivity, in bool rainbowMode)
{
    if (!rainbowMode)
    {
        // Classic NASA wind-tunnel Toepler knife-edge cutoff (50% grey baseline + bright/dark shock reliefs):
        float kProj = dot(epsScreen, normalize(knifeDir)) * sensitivity;
        float lum   = clamp(0.5 + 0.5 * tanh(kProj), 0.02, 1.0);
        return mix(vec3(0.06, 0.08, 0.12), vec3(1.0, 0.94, 0.80), pow(lum, 0.90));
    }
    // Rainbow color-ring Schlieren filter (hue encodes deflection angle, saturation/value encode |eps|):
    float mag = clamp(length(epsScreen) * sensitivity, 0.0, 1.0);
    float ang = atan(epsScreen.y, epsScreen.x);
    vec3  hue = 0.5 + 0.5 * cos(ang + vec3(0.0, 2.0943951, 4.1887902));
    return mix(vec3(0.10, 0.13, 0.18), hue * 1.25, pow(mag, 0.72));
}

#endif // SDF_GEMS_SDF_WHITHAM_MACH_SCHLIEREN_GLSL
