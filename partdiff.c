/*
 * partdiff_futhark - a partdiff port to Futhark by Ruben Felgenhauer.
 *
 * Original partdiff:
 * Copyright (C) 1997 Thomas Ludwig
 * Copyright (C) 1997 Thomas A. Zochler
 * Copyright (C) 1997 Andreas C. Schmidt
 * Copyright (C) 2007-2010 Julian M. Kunkel
 * Copyright (C) 2010-2021 Michael Kuhn
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <malloc.h>
#include <string.h>
#include <sys/time.h>

#include "partdiff_futhark.h"

#define MAX_INTERLINES    100000
#define MAX_ITERATION     200000
#define MAX_THREADS       1024
#define METH_GAUSS_SEIDEL 1
#define METH_JACOBI       2
#define FUNC_F0           1
#define FUNC_FPISIN       2
#define TERM_ACC          1
#define TERM_ITER         2

struct calculation_arguments
{
	uint64_t               N;      /* number of spaces between lines (lines=N+1) */
	double                 h;      /* length of a space between two lines */
	double*                M;      /* host copy of the result matrix */
	struct futhark_f64_2d* matrix; /* current matrix on the Futhark side */
};

struct futhark_state
{
	struct futhark_context_config* cfg;
	struct futhark_context*        ctx;
};

struct calculation_results
{
	uint64_t stat_iteration; /* number of current iteration */
	double   stat_accuracy;  /* actual accuracy of all slaves in iteration */
};

struct options
{
	uint64_t number;         /* Number of threads */
	uint64_t method;         /* Gauss Seidel or Jacobi method of iteration */
	uint64_t interlines;     /* matrix size = interlines*8+9 */
	uint64_t pert_func;      /* perturbation function */
	uint64_t termination;    /* termination condition */
	uint64_t term_iteration; /* terminate if iteration number reached */
	double   term_accuracy;  /* terminate if accuracy reached */
};

struct timeval start_time; /* time when program started */
struct timeval comp_time;  /* time when calculation completed */

static void
usage(char* name)
{
	printf("Usage: %s [num] [method] [lines] [func] [term] [acc/iter]\n", name);
	printf("\n");
	printf("  - num:       number of threads (1 .. %d)\n", MAX_THREADS);
	printf("  - method:    calculation method (1 .. 2)\n");
	printf("                 %1d: Gauß-Seidel\n", METH_GAUSS_SEIDEL);
	printf("                 %1d: Jacobi\n", METH_JACOBI);
	printf("  - lines:     number of interlines (0 .. %d)\n", MAX_INTERLINES);
	printf("                 matrixsize = (interlines * 8) + 9\n");
	printf("  - func:      perturbation function (1 .. 2)\n");
	printf("                 %1d: f(x,y) = 0\n", FUNC_F0);
	printf("                 %1d: f(x,y) = 2 * pi^2 * sin(pi * x) * sin(pi * y)\n", FUNC_FPISIN);
	printf("  - term:      termination condition (1 .. 2)\n");
	printf("                 %1d: sufficient accuracy\n", TERM_ACC);
	printf("                 %1d: number of iterations\n", TERM_ITER);
	printf("  - acc/iter:  depending on term:\n");
	printf("                 accuracy:   1e-4 .. 1e-20\n");
	printf("                 iterations:    1 .. %d\n", MAX_ITERATION);
	printf("\n");
	printf("Example: %s 1 2 100 1 2 100 \n", name);
}

static void
askParams(struct options* options, int argc, char** argv)
{
	int ret;

	if (argc < 7 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-?") == 0)
	{
		usage(argv[0]);
		exit(0);
	}

	ret = sscanf(argv[1], "%" SCNu64, &(options->number));

	if (ret != 1 || !(options->number >= 1 && options->number <= MAX_THREADS))
	{
		usage(argv[0]);
		exit(1);
	}

	ret = sscanf(argv[2], "%" SCNu64, &(options->method));

	if (ret != 1 || !(options->method == METH_GAUSS_SEIDEL || options->method == METH_JACOBI))
	{
		usage(argv[0]);
		exit(1);
	}

	ret = sscanf(argv[3], "%" SCNu64, &(options->interlines));

	if (ret != 1 || !(options->interlines <= MAX_INTERLINES))
	{
		usage(argv[0]);
		exit(1);
	}

	ret = sscanf(argv[4], "%" SCNu64, &(options->pert_func));

	if (ret != 1 || !(options->pert_func == FUNC_F0 || options->pert_func == FUNC_FPISIN))
	{
		usage(argv[0]);
		exit(1);
	}

	ret = sscanf(argv[5], "%" SCNu64, &(options->termination));

	if (ret != 1 || !(options->termination == TERM_ACC || options->termination == TERM_ITER))
	{
		usage(argv[0]);
		exit(1);
	}

	if (options->termination == TERM_ACC)
	{
		ret = sscanf(argv[6], "%lf", &(options->term_accuracy));

		options->term_iteration = MAX_ITERATION;

		if (ret != 1 || !(options->term_accuracy >= 1e-20 && options->term_accuracy <= 1e-4))
		{
			usage(argv[0]);
			exit(1);
		}
	}
	else
	{
		ret = sscanf(argv[6], "%" SCNu64, &(options->term_iteration));

		options->term_accuracy = 0;

		if (ret != 1 || !(options->term_iteration >= 1 && options->term_iteration <= MAX_ITERATION))
		{
			usage(argv[0]);
			exit(1);
		}
	}
}

static void
initVariables(struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options)
{
	arguments->N = (options->interlines * 8) + 9 - 1;
	arguments->h = 1.0 / arguments->N;

	results->stat_iteration = 0;
	results->stat_accuracy  = 0;
}

static void
checkFuthark(struct futhark_context* ctx, int ret)
{
	if (ret != 0)
	{
		char* error = futhark_context_get_error(ctx);
		fprintf(stderr, "Futhark error: %s\n", error != NULL ? error : "unknown error");
		free(error);
		exit(1);
	}
}

static void
initFuthark(struct futhark_state* state)
{
	state->cfg = futhark_context_config_new();
	state->ctx = futhark_context_new(state->cfg);

	char* error = futhark_context_get_error(state->ctx);

	if (error != NULL)
	{
		fprintf(stderr, "Futhark error: %s\n", error);
		free(error);
		exit(1);
	}
}

static void
deinitFuthark(struct futhark_state* state)
{
	futhark_context_free(state->ctx);
	futhark_context_config_free(state->cfg);
}

static void
freeMatrices(struct futhark_state* state, struct calculation_arguments* arguments)
{
	free(arguments->M);
	futhark_free_f64_2d(state->ctx, arguments->matrix);
}

static void*
allocateMemory(size_t size)
{
	void* p;

	if ((p = malloc(size)) == NULL)
	{
		printf("Memory error! (%" PRIu64 " Bytes requested)\n", size);
		exit(1);
	}

	return p;
}

static void
allocateResultMatrix(struct calculation_arguments* arguments)
{
	arguments->M = allocateMemory((arguments->N + 1) * (arguments->N + 1) * sizeof(double));
}

static void
initMatrices(struct futhark_state* state, struct calculation_arguments* arguments, struct options const* options)
{
	int const ret = futhark_entry_init_matrices(state->ctx, &arguments->matrix, (int64_t)options->interlines,
	                                            (int64_t)options->pert_func);
	checkFuthark(state->ctx, ret);
}

static void
calculate(struct futhark_state* state, struct calculation_arguments* arguments, struct calculation_results* results,
          struct options const* options)
{
	double const acc_iter = (options->termination == TERM_ITER) ? (double)options->term_iteration : options->term_accuracy;

	struct futhark_opaque_tup3_arr2d_t_t_i64* calc_result;
	int                                       ret = futhark_entry_calculate(state->ctx, &calc_result, (int64_t)options->method, (int64_t)options->pert_func,
	                                                                        (int64_t)options->termination, acc_iter, arguments->matrix);
	checkFuthark(state->ctx, ret);
	ret = futhark_context_sync(state->ctx);
	checkFuthark(state->ctx, ret);

	struct futhark_f64_2d* result_matrix;
	double                 residuum;
	int64_t                iterations;
	ret = futhark_project_opaque_tup3_arr2d_t_t_i64_0(state->ctx, &result_matrix, calc_result);
	checkFuthark(state->ctx, ret);
	ret = futhark_project_opaque_tup3_arr2d_t_t_i64_1(state->ctx, &residuum, calc_result);
	checkFuthark(state->ctx, ret);
	ret = futhark_project_opaque_tup3_arr2d_t_t_i64_2(state->ctx, &iterations, calc_result);
	checkFuthark(state->ctx, ret);

	results->stat_iteration = (uint64_t)iterations;
	results->stat_accuracy  = residuum;

	ret = futhark_values_f64_2d(state->ctx, result_matrix, arguments->M);
	checkFuthark(state->ctx, ret);
	ret = futhark_context_sync(state->ctx);
	checkFuthark(state->ctx, ret);

	futhark_free_f64_2d(state->ctx, result_matrix);
	futhark_free_opaque_tup3_arr2d_t_t_i64(state->ctx, calc_result);
}

static void
displayStatistics(struct calculation_results const* results, struct options const* options)
{
	double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

	printf("Calculation time:       %f s\n", time);
	// Note: We can't get the actual memory consumption from the Futhark library; print 0 instead.
	printf("Memory usage:           %f MiB\n", 0.0);
	printf("Calculation method:     ");

	if (options->method == METH_GAUSS_SEIDEL)
	{
		printf("Gauß-Seidel");
	}
	else if (options->method == METH_JACOBI)
	{
		printf("Jacobi");
	}

	printf("\n");
	printf("Interlines:             %" PRIu64 "\n", options->interlines);
	printf("Perturbation function:  ");

	if (options->pert_func == FUNC_F0)
	{
		printf("f(x,y) = 0");
	}
	else if (options->pert_func == FUNC_FPISIN)
	{
		printf("f(x,y) = 2 * pi^2 * sin(pi * x) * sin(pi * y)");
	}

	printf("\n");
	printf("Termination:            ");

	if (options->termination == TERM_ACC)
	{
		printf("Required accuracy");
	}
	else if (options->termination == TERM_ITER)
	{
		printf("Number of iterations");
	}

	printf("\n");
	printf("Number of iterations:   %" PRIu64 "\n", results->stat_iteration);
	printf("Residuum:               %e\n", results->stat_accuracy);
	printf("\n");
}

static void
displayMatrix(struct calculation_arguments* arguments, struct calculation_results* results, struct options* options)
{
	int x, y;

	int const interlines = options->interlines;
	int const N          = arguments->N;

	(void)results;

	typedef double (*matrix)[N + 1];

	matrix Matrix = (matrix)arguments->M;

	printf("Matrix:\n");

	for (y = 0; y < 9; y++)
	{
		for (x = 0; x < 9; x++)
		{
			printf("%7.4f", Matrix[y * (interlines + 1)][x * (interlines + 1)]);
		}

		printf("\n");
	}

	fflush(stdout);
}

int
main(int argc, char** argv)
{
	struct options               options;
	struct calculation_arguments arguments;
	struct calculation_results   results;
	struct futhark_state         futhark;

	askParams(&options, argc, argv);

	initVariables(&arguments, &results, &options);

	initFuthark(&futhark);

	allocateResultMatrix(&arguments);
	initMatrices(&futhark, &arguments, &options);

	gettimeofday(&start_time, NULL);
	calculate(&futhark, &arguments, &results, &options);
	gettimeofday(&comp_time, NULL);

	displayStatistics(&results, &options);
	displayMatrix(&arguments, &results, &options);

	freeMatrices(&futhark, &arguments);
	deinitFuthark(&futhark);

	return 0;
}
