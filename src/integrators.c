#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*
Just in case you really want to get that 10th of a second of your life back on the array math.
#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
#include <x86intrin.h>
#endif
#include <omp.h>
*/
#include "bhtree.h"
#include "integrators.h"
#include "constants.h" // Need to create this


/*
	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	INTEGRATORS

	Authors:
		Ramsey Karim
		Harrison Agrusa
		Julian Marohnic

	Description:
		This file contains implementations for integrators using
			"phase vectors" such that particle phase data is stored
			as [x, xdot, y, ydot, ..]. It can be cast to a double array.
		The fundamental gravity differential is implemented here.
		The major assumption necessary for this to work is that
			in the phase vector, dimensional pairs of [position, velocity]
			are stored next to each other, and lined up one after
			the other. For DIM dimensions, there should be DIM*2 elements
			to the phase vector, alternating position-velocity.
		All the tree-calling is in here.
		The integrators supported are Runge-Kutta 4th order and
			Leapfrog (2nd order).

	Usage and Reminders:
		Used as a library, called elsewhere.

		typedef union PhaseVector {
			double coords[DIM*2];
		} Phase;

		For integrators/calculations, assume you have a pointer to
			the phase vector of the current particle and a pointer
			to the other particle's position

	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/

/*
	Reference material..

	void calc_acc(float * dest_x1, float * dest_y1, float * dest_z1,
		float * dest_x2, float * dest_y2, float * dest_z2,
		float * x1, float * y1, float * z1,
		float * x2, float * y2, float * z2,
		int n, float s) {
		float soft = s / pow((float) n, 1./3.);
		float xdiff, ydiff, zdiff;
		xdiff = (*x2) - (*x1);
		ydiff = (*y2) - (*y1);
		zdiff = (*z2) - (*z1);
		float total_dist = (xdiff*xdiff) + (ydiff*ydiff) + (zdiff*zdiff) + (soft*soft);
		total_dist = sqrt(total_dist) * total_dist;
		*dest_x1 = xdiff/total_dist;
		*dest_x2 = -(*dest_x1);
		*dest_y1 = ydiff/total_dist;
		*dest_y2 = -(*dest_y1);
		*dest_z1 = zdiff/total_dist;
		*dest_z2 = -(*dest_z1);
	}

*/


void accumulate_acceleration(Phase * differential, Phase * p, double * otherBody, double otherMass) {
	// Only modify the acceleration terms (in the velocity spots of Phase vector)
	int d;
	double differences[DIM], total_distance = 0, temp;
	for (d = 0; d < DIM; d++) {
		temp = otherBody[d] - p->coords[d*2];
		differences[d] = temp;
		// total_distance accumulated as r^2
		total_distance += temp * temp;
	}
	// Add in the softening parameter
	total_distance += SOFT_SQUARED;
	// Now multiply r^2 by a factor of r for r^3
	total_distance *= sqrt(total_distance);
	for (d = 0; d < DIM; d++) {
		// Units check out
		differential->coords[d*2 + 1] += otherMass * differences[d] / total_distance;
	}
}

void finalize_acceleration(Phase * differential, Phase * p) {
	// Assuming we've accumulated "accelerations" from every single particle already
	for (int d = 0; d < DIM; d++) {
		differential->coords[d*2] = p->coords[d*2 + 1];
		differential->coords[d*2 + 1] *= BIG_G;
	}
}

// Should we merge phasevector.c and integrators.c?
// Yes and we just did! integrators.c is below

/*
	So integrators are gonna look pretty ugly here.
	The integrators will call the trees.
	Higher order integrators will call multiple
		trees per step (RK4), which is kind of a
		bummer, but what can ya do
	These functions should primarily calling
		functions for PhaseVectors, defined in
		phasevector.h, so if you need more functionality,
		make it there and call it here
*/


Phase * leapfrog(Phase * y0_array, double * mass_array, double h, int n, double root_width) {
	// Leapfrog Step
	int i;
	// Allocate a differential array
	Phase * y_prime_array = (Phase *) calloc(n, sizeof(Phase));
	if (!y_prime_array) fatal("Error in LF.");
	// Cast arrays to SIMD-friendly type as well as double array
	double * d_y_prime_array = (double *) y_prime_array;
	double * d_y0_array = (double *) y0_array;

    //do midstep:
	for (i = 0; i < n*DIM*2; i+=2) {
		d_y0_array[i] += d_y0_array[i+1]*h/2.;
	}
	// Make the tree
    TNode_t * root = create_root(root_width);
    grow_tree(root, n, mass_array, y0_array);
	// This is the time consuming step
	for (i = 0; i < n; i++) {
		gather_forces(root, y_prime_array, y0_array, mass_array, i);
	}
	free(uproot(root));
	// Step the velocities using those accelerations
	// Then step the positions using those velocities
	for (i = 0; i < n*DIM*2; i+=2) {
		d_y0_array[i+1] += d_y_prime_array[i+1]*h;
		d_y0_array[i] += d_y0_array[i+1]*h/2.; // Oh boy you better make sure this makes any sense..
	}
	free(y_prime_array);
	return y0_array;
}

Phase * RK4(Phase * y0_array, double * mass_array, double h, int n, double root_width) {
    // RK4 Step
	/*
    Pseudo-code for RK4:
	    k1 = h*fprime(x, t)
	    k2 = h*fprime(x + 0.5*k1, t + 0.5*h)
	    k3 = h*fprime(x + 0.5*k2, t + 0.5*h)
	    k4 = h*fprime(x + k3, t + h)
	    return x + k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0

    */
    int i;
    // Allocate k1, k2, k3, k4 arrays and dummy array
    Phase * k1 = (Phase *) calloc(n, sizeof(Phase));
    Phase * k2 = (Phase *) calloc(n, sizeof(Phase));
    Phase * k3 = (Phase *) calloc(n, sizeof(Phase));
    Phase * k4 = (Phase *) calloc(n, sizeof(Phase));
    Phase * dummy = (Phase *) calloc(n, sizeof(Phase));
    if (!k1 || !k2 || !k3 || !k4 || !dummy) fatal("Error in RK4.");
    // Cast arrays to SIMD-friendly type as well as double array
    double * d_k1 = (double *) k1;
    double * d_k2 = (double *) k2;
    double * d_k3 = (double *) k3;
    double * d_k4 = (double *) k4;
    double * d_dummy = (double *) dummy;
    double * d_y0_array = (double *) y0_array;
    // K1
    rk4_helper(n, mass_array, k1, y0_array, root_width);
    // Prep for k2
    for (i = 0; i < n*DIM*2; i++) {
        d_k1[i] *= h;
        d_dummy[i] = 0.5*d_k1[i] + d_y0_array[i];
    }
    // K2
    // Use dummy variable
    rk4_helper(n, mass_array, k2, dummy, root_width);
    // Prep for k3
    for (i = 0; i < n*DIM*2; i++) {
        d_k2[i] *= h;
        d_dummy[i] = 0.5*d_k2[i] + d_y0_array[i];
    }
    // K3
    rk4_helper(n, mass_array, k3, dummy, root_width);
    // Prep for k4
    for (i = 0; i < n*DIM*2; i++) {
        d_k3[i] *= h;
        d_dummy[i] = d_k3[i] + d_y0_array[i];
    }
    // K4
    rk4_helper(n, mass_array, k4, dummy, root_width);
    // Compile all k's into final answer
    for (i = 0; i < n*DIM*2; i++) {
        d_k4[i] *= h;
        d_y0_array[i] += d_k1[i]/6.0 + d_k2[i]/3.0 + d_k3[i]/3.0 + d_k4[i]/6.0;
    }
    // Remember to freeeeeeee youuuuuuur malloccssssssss
    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(dummy);
    return y0_array;
}

void rk4_helper(int n, double * mass_array,
	Phase * destination_phases, Phase * source_phases,
	double root_width) {
	// Helper function for each k-step of the RK4 function above
	TNode_t * root = create_root(root_width);
	grow_tree(root, n, mass_array, source_phases);
	for (int i = 0; i < n; i++) {
		gather_forces(root, destination_phases, source_phases, mass_array, i);
	}
	free(uproot(root));
}


/*
	More heavily optimized version of RK4 -- speedup isn't much, unfortunately

Phase * RK4(Phase * y0_array, double * mass_array, double h, int n, double root_width) {
    // RK4 Step
    // Pseudo-code for RK4:
	   //  k1 = h*fprime(x, t)
	   //  k2 = h*fprime(x + 0.5*k1, t + 0.5*h)
	   //  k3 = h*fprime(x + 0.5*k2, t + 0.5*h)
	   //  k4 = h*fprime(x + k3, t + h)
	   //  return x + k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0

    int i;
    // Allocate k1, k2, k3, k4 arrays and dummy array
    Phase * k1 = (Phase *) calloc(n, sizeof(Phase));
    Phase * k2 = (Phase *) calloc(n, sizeof(Phase));
    Phase * k3 = (Phase *) calloc(n, sizeof(Phase));
    Phase * k4 = (Phase *) calloc(n, sizeof(Phase));
    Phase * dummy = (Phase *) calloc(n, sizeof(Phase));
    if (!k1 || !k2 || !k3 || !k4 || !dummy) fatal("Error in RK4.");
    // Cast arrays to SIMD-friendly type as well as double array
    // double * d_k1 = (double *) k1;
    // double * d_k2 = (double *) k2;
    // double * d_k3 = (double *) k3;
    // double * d_k4 = (double *) k4;
    // double * d_dummy = (double *) dummy;
    // double * d_y0_array = (double *) y0_array;
    __m128d * m_k1 = (__m128d *) k1;
    __m128d * m_k2 = (__m128d *) k2;
    __m128d * m_k3 = (__m128d *) k3;
    __m128d * m_k4 = (__m128d *) k4;
    __m128d * m_dummy = (__m128d *) dummy;
    __m128d * m_y0_array = (__m128d *) y0_array;
    __m128d m_h = _mm_set_pd(h, h);
    __m128d m_half = _mm_set_pd(0.5, 0.5);
    // K1
    rk4_helper(n, mass_array, k1, y0_array, root_width);
    // Prep for k2
    #pragma omp parallel for schedule(dynamic)
    for (i = 0; i < n*DIM*2*sizeof(double)/sizeof(__m128d); i++) {
        m_k1[i] = _mm_mul_pd(m_k1[i], m_h);
        m_dummy[i] = _mm_add_pd(_mm_mul_pd(m_k1[i], m_half), m_y0_array[i]);
    }
    // K2
    // Use dummy variable
    rk4_helper(n, mass_array, k2, dummy, root_width);
    // Prep for k3
    #pragma omp parallel for schedule(dynamic)
    for (i = 0; i < n*DIM*2*sizeof(double)/sizeof(__m128d); i++) {
        m_k2[i] = _mm_mul_pd(m_k2[i], m_h);
        m_dummy[i] = _mm_add_pd(_mm_mul_pd(m_k2[i], m_half), m_y0_array[i]);
    }
    // K3
    rk4_helper(n, mass_array, k3, dummy, root_width);
    // Prep for k4
    #pragma omp parallel for schedule(dynamic)
    for (i = 0; i < n*DIM*2*sizeof(double)/sizeof(__m128d); i++) {
        m_k3[i] = _mm_mul_pd(m_k3[i], m_h);
        m_dummy[i] = _mm_add_pd(m_k3[i], m_y0_array[i]);
    }
    // K4
    rk4_helper(n, mass_array, k4, dummy, root_width);
    // Compile all k's into final answer
    #pragma omp parallel for schedule(dynamic)
    for (i = 0; i < n*DIM*2*sizeof(double)/sizeof(__m128d); i++) {
        m_k4[i] = _mm_mul_pd(m_k4[i], m_h);
        m_y0_array[i] = _mm_add_pd(_mm_add_pd(
        	_mm_add_pd(_mm_div_pd(m_k1[i], _mm_set_pd(6.0, 6.0)), _mm_div_pd(m_k2[i], _mm_set_pd(3.0, 3.0))),
        	_mm_add_pd(_mm_div_pd(m_k3[i], _mm_set_pd(3.0, 3.0)), _mm_div_pd(m_k4[i], _mm_set_pd(6.0, 6.0)))
        	), m_y0_array[i]);
    }
    // Remember to freeeeeeee youuuuuuur malloccssssssss
    free(k1), free(k2), free(k3), free(k4), free(dummy);
    return y0_array;
}
*/


