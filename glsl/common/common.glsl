// ============================================================================
// sdf-gems: Common Math & Vector Helpers
// ============================================================================
#ifndef SDF_GEMS_COMMON_GLSL
#define SDF_GEMS_COMMON_GLSL

const float SDF_PI  = 3.14159265358979323846;
const float SDF_TAU = 6.28318530717958647692;
const float SDF_PHI = 1.61803398874989484820; // Golden ratio (1 + sqrt(5)) / 2

// @author Inigo Quilez (https://iquilezles.org/articles/distfunctions/)
float dot2( in vec2 v ) { return dot(v, v); }
float dot2( in vec3 v ) { return dot(v, v); }
float ndot( in vec2 a, in vec2 b ) { return a.x*b.x - a.y*b.y; }

// @author Mercury Demogroup (https://mercury.sexy/hg_sdf/)
float vmax( in vec2 v ) { return max(v.x, v.y); }
float vmax( in vec3 v ) { return max(v.x, max(v.y, v.z)); }
float vmin( in vec2 v ) { return min(v.x, v.y); }
float vmin( in vec3 v ) { return min(v.x, min(v.y, v.z)); }

#endif // SDF_GEMS_COMMON_GLSL
