#ifndef INTEGRATORS_H
#define INTEGRATORS_H

#include "bhtree.h"

// #if defined(_MSC_VER)
// #include <intrin.h>
// #elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
// #include <x86intrin.h>
// #endif


// Using Phase, we can jump through arrays of [x, xdot, y, ydot, ... ] in a single bound
typedef union PhaseVector {
	double coords[DIM*2];
	// Including this line should force this union to be half-word-aligned.
	// __m128d simd_vector[DIM*2*sizeof(double)/sizeof(__m128d)];
} Phase;

void accumulate_acceleration(Phase * differential, Phase * p, double * otherBody, double otherMass);

void finalize_acceleration(Phase * differential, Phase * p);

Phase * euler(Phase * y0_array, double * mass_array, double h, int n, double root_width);

Phase * leapfrog(Phase * y0_array, double * mass_array, double h, int n, double root_width);

Phase * RK4(Phase * y0_array, double * mass_array, double h, int n, double root_width);

void rk4_helper(int n, double * mass_array, Phase * destination_phases, Phase * source_phases, double root_width);

#endif