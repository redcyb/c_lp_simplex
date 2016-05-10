#ifndef __LINEAR_PROGRAM_H__
#define __LINEAR_PROGRAM_H__

#include "matrix.h"

struct lp;

typedef struct lp * LP1;

Matrix solve(LP1 p);

LP1 read_LP();

void print_LP(LP1 P);

void destroy_LP(LP1 P);

void vis_init(int *argc, char **argv);

void visualize(LP1 P);

#endif
