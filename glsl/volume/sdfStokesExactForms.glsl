// sdfStokesExactForms: Exterior differential forms Omega^k(R^3) & Stokes-exact ray/beam integrals.
// Collapses 1D ray and 2D pixel-cone volumetric integrals to 0D/1D boundary evaluations. See docs/volume/sdfStokesExactForms.md.
#ifndef SDF_GEMS_SDF_STOKES_EXACT_FORMS_GLSL
#define SDF_GEMS_SDF_STOKES_EXACT_FORMS_GLSL

// Estimate closest-approach parameter uStar in [0, 1] along slab [ta, tb]
// from endpoint SDF+gradient values dgA = (Fa, gradFa), dgB = (Fb, gradFb).
float sdfStokesTurningPointU(vec4 dgA, vec4 dgB, vec3 rd) {
    float dfa = dot(dgA.yzw, rd);
    float dfb = dot(dgB.yzw, rd);
    if (dfa < 0.0 && dfb > 0.0) {
        return clamp(-dfa / max(dfb - dfa, 1e-7), 0.0, 1.0);
    }
    return (dgA.x < dgB.x) ? 0.0 : 1.0;
}

// Estimate minimum SDF value F_* on [ta, tb] using cubic Hermite interpolation at uStar.
float sdfStokesMinDistance(vec4 dgA, vec4 dgB, vec3 rd, float dt) {
    float u = sdfStokesTurningPointU(dgA, dgB, rd);
    float dfa = dot(dgA.yzw, rd) * dt;
    float dfb = dot(dgB.yzw, rd) * dt;
    float u2 = u * u, u3 = u2 * u;
    float h00 = 2.0 * u3 - 3.0 * u2 + 1.0;
    float h10 = u3 - 2.0 * u2 + u;
    float h01 = -2.0 * u3 + 3.0 * u2;
    float h11 = u3 - u2;
    float fStar = h00 * dgA.x + h10 * dfa + h01 * dgB.x + h11 * dfb;
    return min(fStar, min(dgA.x, dgB.x));
}

// 0-step Stokes-exact slab integrals for 3 physical shell potentials:
// Returns vec3(exponentialShock, lorentzianRim, dipoleTangentLobe) >= 0.
vec3 sdfStokesExactRaySlab(vec4 dgA, vec4 dgB, vec3 rd, float dt, float beta, float gamma) {
    if (dt <= 0.0) return vec3(0.0);

    float fStar = sdfStokesMinDistance(dgA, dgB, rd, dt);

    // 1. Exponential shock shell: Phi_0 = exp(-beta * max(F, 0))
    float expA = exp(-beta * max(dgA.x, 0.0));
    float expB = exp(-beta * max(dgB.x, 0.0));
    float expS = exp(-beta * max(fStar, 0.0));
    float intExp = max(2.0 * expS - expA - expB, 0.0);

    // 2. Algebraic Lorentzian rim shell: Phi_0 = 1 / sqrt(F^2 + gamma^2)
    float g2 = max(gamma * gamma, 1e-8);
    float rimA = inversesqrt(dgA.x * dgA.x + g2);
    float rimB = inversesqrt(dgB.x * dgB.x + g2);
    float rimS = inversesqrt(fStar * fStar + g2);
    float intRim = max(2.0 * rimS - rimA - rimB, 0.0);

    // 3. Exact dipole tangent-lobe jump: Phi_0 = (gradF . d) / sqrt(F^2 + gamma^2)
    float dipA = dot(dgA.yzw, rd) * rimA;
    float dipB = dot(dgB.yzw, rd) * rimB;
    float intDipole = max(dipB - dipA, 0.0);

    return vec3(intExp, intRim, intDipole);
}

// Evaluate the 1-form vector potential A(x) = psi(F) * (gradF x swirlAxis)
// whose exterior derivative B = curl(A) is divergence-free and concentrated on {F = 0}.
vec3 sdfStokesVectorPotential1Form(vec4 dg, vec3 swirlAxis, float gamma) {
    float psi = inversesqrt(dg.x * dg.x + max(gamma * gamma, 1e-8));
    return psi * cross(dg.yzw, swirlAxis);
}

// Stokes' Theorem pixel-frustum triangle flux:
// Integrates iint_Triangle curl(A) . dS = oint_dTriangle A . dx
// given (F, gradF) and positions at the 3 corner rays (p0, p1, p2).
float sdfStokesTriangleBeamFlux(
    vec3 p0, vec4 dg0,
    vec3 p1, vec4 dg1,
    vec3 p2, vec4 dg2,
    vec3 swirlAxis, float gamma
) {
    vec3 A0 = sdfStokesVectorPotential1Form(dg0, swirlAxis, gamma);
    vec3 A1 = sdfStokesVectorPotential1Form(dg1, swirlAxis, gamma);
    vec3 A2 = sdfStokesVectorPotential1Form(dg2, swirlAxis, gamma);

    float e01 = 0.5 * dot(A0 + A1, p1 - p0);
    float e12 = 0.5 * dot(A1 + A2, p2 - p1);
    float e20 = 0.5 * dot(A2 + A0, p0 - p2);
    return e01 + e12 + e20;
}

#endif // SDF_GEMS_SDF_STOKES_EXACT_FORMS_GLSL
