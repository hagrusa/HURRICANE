#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <stdint.h>
#include "bhtree.h"
#include "linkedlist.h"
#include "integrators.h"
#include "constants.h"

/*
	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	BHTREE

	Authors:
		Ramsey Karim
		Harrison Agrusa
		Julian Marohnic

	Description:
		Contains all Barnes-Hut tree functionality for the N-body code.

	Usage:
		Used as a library, called elsewhere.

	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/

// >>>>> GOLDEN RULES: <<<<<

// ALL ARRAYS DESCRIBING POSITIONS OR VELOCITIES ARE OF LENGTH DIM
// This implies that a vector [x, xdot, y, ydot, ..] is of length 2*DIM

// INFORMATION ABOUT THE iTH DIMENSION IS STORED IN THE DIMENSION ID AT
// THE (info << i) BIT
// TO RECOVER IT, ((dim_id >> i) & 1)

// Basic safe alloc function with print statement
int safe_alloc(void * p) {
	if (!p) {
		printf("ALLOCATION ERROR\n");
		return -1;
	} else {
		return 0;
	}
}

//A function to display an error message and then exit
void fatal(char *message) {
    char error_message[100];
    strcpy(error_message, "[!!] Fatal Error: ");
    strncat(error_message, message, 83); //83 = 100 - 17, max number of characters to append
    perror(error_message);
    exit(-1);
}

void copy_position(double * destination, double * source) {
	// Assuming DIM-element array of doubles
	for (int i = 0; i < DIM; i++) {
		destination[i] = source[i];
	}
}

void init_position_zero(double * destination) {
	// Assuming DIM-element array of doubles
	for (int i = 0; i < DIM; i++) {
		destination[i] = 0;
	}
}

TNode_t * create_node(double width, double * center) {
	// Check inputs (lite check)
	if (width < 0 || !center) return NULL;
	// Safe allocation
	TNode_t * newNode = (TNode_t*) malloc(sizeof(TNode_t));
	if (safe_alloc((void*) newNode) < 0) return NULL;
	newNode->child = NULL;
	newNode->occupancy = 0;
	newNode->occupant_if_leaf = -1;
	newNode->width = width;
	newNode->mass = 0;
	copy_position(newNode->center, center);
	// Not initializing COM because it should be dealt with immediately anyway
	// Or should we...?
	init_position_zero(newNode->com);
	return newNode;
}

TNode_t * create_root(double width) {
	if (width < 0) return NULL;
	double origin[3]; // Center on origin by default
	init_position_zero(origin);
	TNode_t * root = create_node(width, origin);
	return root; // NULL if error; manage upstream
}

TNode_t * remove_node(TNode_t * n) {
	if (n) {
		if (n->child) {
			for (int i = 0; i < N_CHILD; i++) {
				free(remove_node(n->child[i]));
			}
		}
		free(n->child);
	}
	return n;
}

TNode_t * uproot(TNode_t * root) {
	remove_node(root);
	return root;
}

TNode_t ** create_child_array() {
	TNode_t ** child_array;
	child_array = (TNode_t **) malloc(N_CHILD*sizeof(TNode_t *));
	if (safe_alloc((void *) child_array) < 0) return NULL;
	return child_array;
}

TNode_t * grow_tree(TNode_t * root, int n, double * mass_array, Phase * phase_array) {
	// Basic input checking
	if (!root || n < 1 || !mass_array || !phase_array) {
		return NULL;
	}
	// We should make a linked list
	LL * point_ll = int_array_to_ll(n);
	// Nice!

	// Ok now we should take that list and the current node (root) and start dividing
	// Call grow_branch (recursive)
	grow_branch(root, mass_array, phase_array, point_ll);
	// Delete the list
	ll_delete_list(point_ll);
	// That's it! Tree is created!
	return root;
}

void grow_branch(TNode_t * parent_node, double * mass_array, Phase * phase_array, LL * parentList) {
	/*
		Big ol' recursive function for growing trees
		This should handle child creation and sorting of particles
		It will also populate the COM, branch mass, occupancy attributes at every node,
			which will recursively fill the tree with all the information/shortcuts we need
		The parent_node should have center and width filled out; all else will be dealt with here
		COM, branch mass, occupancy of parent_node should be empty/meaningless at the start of this
			function call. They will be callously overwritten
		Child array will be created, since we know this branch isn't a leaf
	*/
	TNode_t * temp_node;
	uint8_t dimension, dim_id, position_check;
	int index, length;
	LL * branchLists[N_CHILD];
	double test_position[DIM], node_center[DIM];
	if (!parent_node) {
		printf("Growing tree failed. Exiting current call.\n");
		return;
	}
	// Allocate the child array
	if (!(parent_node->child = create_child_array())) return;
	// Get the parent center
	copy_position(node_center, parent_node->center);
	for (dimension = 0; dimension < N_CHILD; dimension++) {
		// Loop through each potential sub-division and initialize a list
		branchLists[dimension] = ll_create_list();
	}
	for (LLn * particleNode = ll_pop(parentList); particleNode != NULL; particleNode = ll_pop(parentList)) {
		// Initialize a clean id
		dim_id = 0;
		index = particleNode->value;
		// Get this particle's position
		get_particle_position(test_position, phase_array, index);
		// Loop through each spatial dimension
		for (dimension = 0; dimension < DIM; dimension++) {
			// Check if the particle is to the left or right of the center
			position_check = test_position[dimension] > node_center[dimension];
			// Record this information in the appropriate bit in the id
			dim_id |= position_check << dimension;
		}
		// Now use this id as an index
		ll_append_node(branchLists[dim_id], particleNode);
	}
	// parent_list is deleted in the call above this one
	// THIS PART IS RECURSIVE
	// Loop through sub-divisions again
	for (dimension = 0; dimension < N_CHILD; dimension++) {
		length = ll_quick_length(branchLists[dimension]);
		if (length == 0) {
			// BASE CASE 0: There are no particles in this cube
			// Delete this list and do not make a child node
			parent_node->child[dimension] = NULL; // Marked to skip
			ll_delete_list(branchLists[dimension]);
		} else if (length == 1) {
			// BASE CASE 1: There is exactly 1 particle in this cube: good!
			// Somehow make a leaf (temp_node->child remains NULL from initialization)
			// Delete this list
			index = branchLists[dimension]->head->value;
			get_new_center(test_position, parent_node->width, parent_node->center, dimension);
			// Now we have the center of the new node
			temp_node = create_node(parent_node->width/2, test_position);
			parent_node->child[dimension] = temp_node;
			temp_node->mass = mass_array[index];
			temp_node->occupancy = 1;
			temp_node->occupant_if_leaf = index;
			get_particle_position(temp_node->com, phase_array, index);
			accumulate_com(parent_node->com, temp_node->com, temp_node->mass);
			parent_node->mass += temp_node->mass;
			parent_node->occupancy += 1;
			// The new leaf node should be complete
			ll_delete_list(branchLists[dimension]); // #%(&@^$@*(&^)*(&Q^R)Q*&^R) THIS IS WHERE YOU LEFT OFF
			// ---------------------------------------------------------------------------^^^^^^^^^^^^^^^^^^
		} else {
			// RECURSIVE CASE: There are more than 1 particles in this cube
			// Create a new tree node
			// Recursively call grow_branch
			/*
				Fuck I'm getting back to this almost a week later
				Ok so here's where the code left off I think
				We have base cases for empty nodes (leave NULL) and 1 node (set as leaf, COM=Position)
				We might have leftover len(list) > 1 lists, so we must handle those
				!! We should be calculating COM and total branch mass in this loop!!!!! !!
				Let's have this part of the function accumulate mass and set the COM of this node
				Ok if we set the COM and mass of *this* node (parent_node), then the traversal is done recursively
				Need to write that into base cases
			*/
			// Copied from above for initializing child node
			get_new_center(test_position, parent_node->width, parent_node->center, dimension);
			// Now we have the center of the new node
			temp_node = create_node(parent_node->width/2, test_position);
			parent_node->child[dimension] = temp_node;
			// Attributes dealt with: center, width			
			// Attributes left: child, occupancy, occupant_if_leaf, mass, com
			// OMG recursion
			grow_branch(temp_node, mass_array, phase_array, branchLists[dimension]);
			// This entire branch is now populated
			// We've handled child (potentially), occupancy, occupant_if_leaf, mass, and com
			accumulate_com(parent_node->com, temp_node->com, temp_node->mass);
			parent_node->mass += temp_node->mass;
			parent_node->occupancy += temp_node->occupancy;
			// List has been used, can delete (this handles parent_list for the recursive call)
			ll_delete_list(branchLists[dimension]);
		}
	}
	// End of sub-cube loop
	finalize_com(parent_node->com, parent_node->mass);
	// That's it I think
	/*
		At this point we've dealt with:
			WIDTH: we came into this function with our width
			CENTER: we came into this function with our center
			CHILD: we knew this node wasn't a leaf, so we allocated & populated a child array
			MASS: we accumulated masses of children, which we assigned as a base case
			OCCUPANCY: we accumulated occupancies as we did masses, assigned to 1 as a base case
			OCCUPANT_IF_LEAF: ours should be -1; we assign indices to our leaf children
			COM: we accumulated our com from our children's com, which were assigned as a base case
		DONE!
	*/
}

void accumulate_com(double * destination, double * position, double mass) {
	// Assuming position and destination are of the correct length (DIM)
	for (int i = 0; i < DIM; i++) {
		destination[i] += position[i] * mass;
	}
}

void finalize_com(double * destination, double total_mass) {
	// Assuming position and destination are of the correct length (DIM)
	for (int i = 0; i < DIM; i++) {
		destination[i] /= total_mass;
	}
}

void get_new_center(double * destination, double parent_width, double * parent_center, uint8_t dim_id) {
/*
	Given a parent center and width and a child (dimension) ID,
	fill destination array with the center of the child node.
*/
	uint8_t temp, mask = 1;
	double width = parent_width/2.;
	for (uint8_t dimension = 0; dimension < DIM; dimension++) {
		// Retreive the relevant information
		temp = (dim_id >> dimension) & mask;
		// Sub-cell center should be parent center +/- (new width / 2)
		// Accomplish this by using temp to toggle adding width
		destination[dimension] = parent_center[dimension] - (width / 2.) + (temp * width);
	}
}

void get_particle_position(double * destination, Phase * phase_array, int index) {
	// It's a union, so casting like this should be safe
	// We know it's really just a double array anyway
	double * this_particle = (double *) (phase_array + index);
	for (int dimension = 0; dimension < DIM; dimension++) {
		// This grabs [x, y, ..] from [x, xdot, y, ydot, ..]
		destination[dimension] = this_particle[dimension * 2];
	}
}

void gather_forces(TNode_t * root, Phase * differential_array, Phase * phase_array, double * mass_array, int index) {
	// Gather things relevant to this particle, since we'll reuse these a lot
	Phase * p = phase_array + index, * differential = differential_array + index;
	// Call helper function for traversal
	if (!root) {
		printf("Your root is missing, we can't get forces from this.\n");
		return;
	}
	force_traverse(root, differential, p, index);
	finalize_acceleration(differential, p);
}

void force_traverse(TNode_t * node, Phase * differential, Phase * p, int index) {
	double apparent_angle;
	TNode_t * child;
	// Recursive traversal method, 2 base cases and 1 recursive case
	// BASE CASE 2: LEAF NODE
	//	Get force from this particle, as long as this particle isn't the test particle
	if (node->occupancy == 1) {
		if (node->occupant_if_leaf == index) return;
		accumulate_acceleration(differential, p, node->com, node->mass);
	} else {
		// ???????? Hey should we use COM or CENTER ???????
		// Consensus is center-of-mass
		apparent_angle = node->width / get_particle_distance(node->com, p);
		if (apparent_angle <= THETA) {
			accumulate_acceleration(differential, p, node->com, node->mass);
		} else {
			for (int i = 0; i < N_CHILD; i++) {
				// NULL check the child before calling
				if ((child = node->child[i])) force_traverse(child, differential, p, index);
			}
		}
	}
}

double get_particle_distance(double * position, Phase * p) {
	double * other_position = (double *) p;
	double distance = 0, temp;
	for (int dimension = 0; dimension < DIM; dimension++) {
		temp = position[dimension] - other_position[dimension * 2];
		distance += temp*temp;
	}
	// Consider using the square to save time, and saving
	//	the square of the critical angle so it's all easier
	return sqrt(distance);
}

void print_traverse(TNode_t * node, char * prelude) {
	if (!node) {
		printf("%sEmpty {}\n", prelude);
		return;
	}
	int leaf = node->occupancy == 1;
	printf("%s%s m(%lf) %c(%d) {",
		prelude,
		leaf ? "Leaf" : "Node",
		node->mass,
		leaf ? 'p' : 'n',
		leaf ? node->occupant_if_leaf : node->occupancy
	);
	if (!leaf) {
		printf("\n");
		// add 2 spaces to prelude
		char * newPrelude = (char *) malloc(sizeof(char)*(strlen(prelude) + 2));
		strcpy(newPrelude, prelude);
		strcpy(newPrelude + strlen(prelude), "|\t");
		for (int i = 0; i < N_CHILD; i++) {
			print_traverse(node->child[i], newPrelude);
		}
		free(newPrelude);
		printf("%s}\n", prelude);
	} else {
		printf("}\n");		
	}
}

void write_tree(const char * fname, TNode_t * root) {
	FILE * output;
	output = fopen(fname, "w");
	// Now traverse
	write_traverse(output, root);
	fclose(output);
}

void write_traverse(FILE * output, TNode_t * node) {
	if (!node) return;
	int i;
	fprintf(output, "\n");
	if (node->occupancy > 1) {
		for (i = 0; i < N_CHILD; i++) {
			write_traverse(output, node->child[i]);
		}
	} else {
		fprintf(output, "%lf ", node->width);
		for (i = 0; i < DIM; i++) {
			fprintf(output, " %lf %lf", node->center[i], node->com[i]);
		}
	}
}
