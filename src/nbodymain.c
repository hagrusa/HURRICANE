#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <stdint.h>
#include <time.h>
#include "bhtree.h"
#include "linkedlist.h"
#include "integrators.h"
#include "constants.h"

/*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    NBODYMAIN

    Authors:
        Ramsey Karim
        Harrison Agrusa
        Julian Marohnic

    Description:
        Main function for the N-body code.
        This handles input and output and maintains the main
            integration loop.
        This file doesn't do any of the heavy lifting. That's all
            spread throughout the other files, especially
            integrators.c and bhtree.c.

    Usage:
        Call nbodymain with the following 8 or 9 arguments:
            <input_file>: string
                Path to the input file. Inputs are expected as text
                files including only numbers. Each row indicates a
                unique particle. The first column is mass. The next
                three columns are x, y, z position, followed by
                the x, y, and z components of velocity, for a total
                of 7 columns.
                This program can support other spatial dimensions,
                but needs small tweaks in input/output. This should
                not be attempted without express consent of the authors.

            <step_size>: float/double
                Integration time step.

            <n_steps>: integer
                Number of integration steps.

            <softening>: float/double
                Gravitational force softening parameter.
                Set to 0 for no softening.

            <output_frequency>: integer
                Frequency at which output files are written.
                Reduce frequency for slightly better runtime,
                especially for Python-based analysis or graphing.
                This number is read in as:
                    "We will save an output file every <output_frequency>
                        steps."
                If you want every step saved, use 1. If you want every
                fourth step, use 4.
                Do not input 0 or negative numbers. The program will ask
                you to try again.

            <integrator_name>: string
                Integrator to use.
                The string must be EXACTLY "RK4" or "LF2".
                No other inputs are acceptable, and the program will
                ask you to try again.

            <critical_angle>: double
                Critical angle to use for the tree approximation.
                Set to a larger number for more approximation and
                potentially quicker integration, at the risk of
                accuracy loss.
                Set to 0 if you would like the program to act as a
                particle-particle simulation.

            <output_directory>: string
                Directory to which to output text files.
                Output is in the same format as input.
                This directory MUST exist prior to calling this
                program.

            <index (optional)>: integer
                You do not need this argument.
                You may input an integer to be used as the starting
                step number. Do not input a negative number.
                If present and nonzero, the program will integrate for
                <n_steps> - <index> steps and label the output files
                accordingly.
                This is useful for starting your program where it left
                off from a past output file. Simply use the step number
                of the previous output as this argument, and the rest
                will be handled.

    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/

double SOFT_SQUARED; //these are global variables!
double THETA;

void write_particles(char * fname, double * mass_array, Phase * phase_array, int n) {
    if (!mass_array || !phase_array) return;
    FILE * output;
    double * d_phase = (double *) phase_array;
    output = fopen(fname, "w");
    for (int i = 0; i < n; i++) {
        fprintf(output, "%lf %lf %lf %lf %lf %lf %lf\n", mass_array[i],
        	// Order has been edited to reflect x y z xdot ydot zdot in/output
            d_phase[i*DIM*2 + 0],
            d_phase[i*DIM*2 + 2],
            d_phase[i*DIM*2 + 4],
            d_phase[i*DIM*2 + 1],
            d_phase[i*DIM*2 + 3],
            d_phase[i*DIM*2 + 5]);
    }
    fclose(output);
}

double appropriate_width(Phase * phase_array, int n) {
    double current_max = 0, width = 1;
    for (int i = 0; i < n; i++) {
        for (int d = 0; d < DIM; d++) {
            current_max = fmax(current_max, fabs(phase_array[i].coords[d*2]));
        }
    }
    while (width < current_max) {
        width *= 2;
    }
    return width * 2;
}


int main(int argc, char ** argv) {
    char *fname = (char *) malloc(50); //text files cannot be more than 50 characters
    char *integrator = (char *) malloc(10); //integrator string
    strcpy(fname, argv[1]);
    FILE *file;
    file = fopen(fname, "r");

    int n = 0, i;
    int start_index = 0; //by default
    double dummy[7];

    if (argc < 9) {
        printf("Incorrect number of command line arguments.\n");
        printf("Usage: ./nbodymain data.txt h Nsteps epsilon outfreq integrator theta_critical output_dir index(optional)\n");
        exit(1);
    }
    else if (argc == 10) {
        start_index = atoi(argv[9]);
        // This lets you start the code at an index other than zero
        // Input the last index that successfully saved; the program will restart
        //  beginning from the index right after that (index + 1)
    }
    else if (argc > 10) {
        printf("Incorrect number of command line arguments.\n");
        printf("Usage: ./nbodymain data.txt h Nsteps epsilon outfreq integrator theta_critical output_dir index(optional)\n");
        exit(1);
    }

    double dt = atof(argv[2]);
    int steps = atoi(argv[3]);
    double epsilon = atof(argv[4]);
    SOFT_SQUARED = pow(epsilon, 2.0); //this is a global!!
    int outfreq = atoi(argv[5]);
    strcpy(integrator, argv[6]);
    THETA = atof(argv[7]); //this is a global

    char output_dir[50]; // Directory names cannot be more than 50 characters
    strcpy(output_dir, argv[8]);  //output directory
    int dir_name_len = strlen(output_dir);
    // ramsey: can you change the file saving stuff so that the output files will be saved to the user specified output directory?


    //error check inputs
    //double hmin = 0.0001; //arbitrary for now
    //double hmax = 100000000000.0;
    if (file == NULL) {
        fatal("Cannot open file");
    }
    //if (h < hmin || h > hmax) {   //should we have a max or min timestep??
    //    fatal("Timestep (h) must be between hmin and hmax");
    //}
    if (steps < 1) {
        fatal("Number of steps must be at least 1");
    }
    if (epsilon < 0.0) {
        fatal("Softening Parameter (epsilon) cannot be less than 0");
    }
    if (outfreq < 1) {
        fatal("Output frequency (outfreq) cannot be less than 1");
    }
    if (strcmp(integrator, "RK4") != 0 && strcmp(integrator, "LF2") != 0) {
        fatal("Integrator must be RK4 or LF2");
    }
    //find number of particles
    while(!feof(file)) {
		if (fscanf(file, "%lf %lf %lf %lf %lf %lf %lf",
			dummy + 0,
			dummy + 1, dummy + 2, dummy + 3,
			dummy + 4, dummy + 5, dummy + 6) == 7) {
        	n++;
		}
    }
    printf("    Number of Particles:            %d\n", n);
    printf("\nReading Input File (%s)......\n", fname);
    free(fname);
    rewind(file); //go back to begining on input file, start reading
    //initialize vectors
    double * mass_array;
    Phase * phase_array;
    mass_array = (double *) malloc(n * sizeof(double));
    phase_array = (Phase *) malloc(n * sizeof(Phase));
    if (safe_alloc(mass_array) + safe_alloc(phase_array) < 0) {
	    free(mass_array);
	    free(phase_array);
        return -1;
    }
    int error_count = n;
    for(i = 0; i < n; i++) {
        if (fscanf(file, "%lf %lf %lf %lf %lf %lf %lf", mass_array + i,
        	// Order has been edited to reflect x y z xdot ydot zdot in/output
        	phase_array[i].coords + 0,
			phase_array[i].coords + 2,
			phase_array[i].coords + 4,
			phase_array[i].coords + 1,
			phase_array[i].coords + 3,
			phase_array[i].coords + 5) == 7) {
			error_count--;
        }
    }
    printf("Successfully read particles with %d errors\n", error_count);
    fclose(file);
    // Test the nbody tree!

    double width; //width of tree
    char output_fname[100]; // (output_dir)+nbodyout_00000.txt
    strcpy(output_fname, output_dir);
    snprintf(output_fname + dir_name_len, sizeof(output_fname), "nbodyout_%06d.txt", 0);
    char * output_fname_number = output_fname + dir_name_len + 9;
    // result/
    // This lets you start the code at an index other than zero
    // Input the last index that successfully saved; the program will restart
    //  beginning from the index right after that (index + 1)
    if (start_index) {
        i = start_index;
    } else {
        i = 0;
    }
    // void get_particle_position(double * destination, Phase * phase_array, int index)

    printf("\nRunning Simulation....\n");
    clock_t start, end;
    double cpu_time_used;
    start = clock();

    //do computation
    for (; i < steps; i++) {
        width = appropriate_width(phase_array, n);
        printf("Iteration %03d/%03d", i+1, steps);

        //choose integrator
        if (strcmp("LF2\0", integrator) == 0) {
            phase_array = leapfrog(phase_array, mass_array, dt, n, width);
        }
        else if (strcmp("RK4\0", integrator) == 0) {
            phase_array = RK4(phase_array, mass_array, dt, n, width);
        }
        //write particles
        if (i % outfreq == 0) {
            snprintf(output_fname_number, 11, "%06d.txt", i);
            printf(" and saved most recently to %s", output_fname);
            write_particles(output_fname, mass_array, phase_array, n);
        }
        printf("\r");
    }


    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\nFinished integrating %d particles over %d steps.\n", n, steps);
    printf("\nTotal Computation Time: %f seconds\n", cpu_time_used);

    free(integrator);
    free(mass_array);
    free(phase_array);
    return 0;
}

