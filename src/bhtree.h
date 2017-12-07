#ifndef BHTREE_H
#define BHTREE_H

#define DIM 3
#define N_CHILD 8 // 1 << DIM

#include <stdint.h>
#include "linkedlist.h"
#include "integrators.h"

// >>>>> GOLDEN RULES: <<<<<

// ALL ARRAYS DESCRIBING POSITIONS OR VELOCITIES ARE OF LENGTH DIM
// This implies that a vector [x, xdot, y, ydot, ..] is of length 2*DIM

// INFORMATION ABOUT THE iTH DIMENSION IS STORED IN THE DIMENSION ID AT
// THE (info << i) BIT
// TO RECOVER IT, ((dim_id >> i) & 1)

typedef struct TreeNode TNode_t;

struct TreeNode {
	// CHILD is either NULL or of size N_CHILD
	// If CHILD is NULL, node is a leaf
	TNode_t ** child;
	int occupancy;
	int occupant_if_leaf;
	double width;
	double mass;
	double center[DIM];
	double com[DIM];
};

// Basic safe alloc function with print statement
int safe_alloc(void * p);

void fatal(char *message);

void copy_position(double * destination, double * source);

void init_position_zero(double * destination);

TNode_t * create_node(double width, double * center);

TNode_t * create_root(double width);

TNode_t * remove_node(TNode_t * n);

TNode_t * uproot(TNode_t * root);

TNode_t ** create_child_array();

TNode_t * grow_tree(TNode_t * root, int n, double * mass_array, Phase * phase_array);

void grow_branch(TNode_t * parent_node, double * mass_array, Phase * phase_array, LL * parentList);

void accumulate_com(double * destination, double * position, double mass);

void finalize_com(double * destination, double total_mass);

void get_new_center(double * destination, double parent_width, double * parent_center, uint8_t dim_id);

void get_particle_position(double * destination, Phase * phase_array, int index);

void gather_forces(TNode_t * root, Phase * differential_array, Phase * phase_array, double * mass_array, int index);

void force_traverse(TNode_t * node, Phase * differential, Phase * p, int index);

double get_particle_distance(double * position, Phase * p);

void print_traverse(TNode_t * node, char * prelude);

void write_tree(const char * fname, TNode_t * root);

void write_traverse(FILE * output, TNode_t * node);

#endif