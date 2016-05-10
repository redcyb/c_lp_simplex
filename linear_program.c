#include <float.h>
#include <stdbool.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "linear_program.h"

typedef enum {
    OPTIMAL, UNBOUNDED, INFEASIBLE
} simp;

typedef enum {
    EQ, LE, GE
} eq;

struct constraint {
    double *coefficients;
    double constant;
    eq type;
};

typedef struct constraint *Constraint;

struct lp {
    bool is_SF;

    // non-SEF data
    Constraint *constraint;
    int num_constraints;
    int variables;
    bool *ge0;
    int num_not_ge0;
    bool obj_max;

    // SEF data
    Matrix A;
    Matrix b;
    Matrix c;
    double z;
};

Constraint create_constraint(int vars) {
    Constraint con = malloc(sizeof(struct constraint));
    con->coefficients = malloc(sizeof(double) * vars);
    return con;
}

void destroy_constraint(Constraint con) {
    free(con->coefficients);
    free(con);
}

void destroy_LP(LP1 P) {
    int i = 0;
    for (i; i < P->num_constraints; i++) {
        destroy_constraint(P->constraint[i]);
    }
    free(P->constraint);
    if (P->is_SF) {
        destroy_matrix(P->A);
        destroy_matrix(P->b);
    }
    destroy_matrix(P->c);
}

LP1 copy_LP(LP1 P) {
    LP1 Q = malloc(sizeof(struct lp));
    Q->A = copy_matrix(P->A);
    Q->b = copy_matrix(P->b);
    Q->c = copy_matrix(P->c);
    Q->z = P->z;
    Q->obj_max = P->obj_max;
    return Q;
}

void canonical_form(LP1 P, int *B) {
    Matrix A_B = take_columns(P->A, B, P->A->rows);
    invert(A_B);
    Matrix A_BT = copy_matrix(A_B);
    transpose(A_BT);
    Matrix c_B = take_rows(P->c, B, P->A->rows);
    Matrix y = mult_new(A_BT, c_B);
    destroy_matrix(A_BT);
    destroy_matrix(c_B);

    transpose(y);

    Matrix z_diff = mult_new(y, P->b);
    P->z += single_to_num(z_diff);
    destroy_matrix(z_diff);


    Matrix c_diff = mult_new(y, P->A);

    multiply_scalar(c_diff, -1);
    transpose(c_diff);
    add_matrix(P->c, c_diff);
    destroy_matrix(c_diff);
    destroy_matrix(y);

    Matrix new_A = mult_new(A_B, P->A);
    Matrix new_b = mult_new(A_B, P->b);
    destroy_matrix(A_B);
    destroy_matrix(P->A);
    destroy_matrix(P->b);

    P->A = new_A;
    P->b = new_b;
}

bool check_basis(Matrix mtr, int *B) {
    Matrix A_B = take_columns(mtr, B, mtr->rows);
    bool res = rref(A_B) == mtr->rows;
    destroy_matrix(A_B);
    return res;
}

// subsets of size n out of numbers 0..m-1
// s : starting index in B
bool basis(Matrix mtr, int *B, int sn, int sm, int n, int m) {
    if (n == 0) {
        return check_basis(mtr, B);
    }
    int i = sm;
    for (i; i < m - n + 1; i++) {
        B[sn] = i;
        if (basis(mtr, B, sn + 1, i, n - 1, m)) {
            return true;
        }
    }
    return false;
}

simp simplex(LP1 P, int *B) {
    canonical_form(P, B);

//    printf("First Simplex :: \n");
//    print_LP(P);

    int k = -1;
    int i = 0;
    for (i; k == -1 && i < P->c->rows; i++) {
        if (P->c->table[i][0] > 0)
            k = i;
    }
    if (k == -1) return OPTIMAL;

    Matrix A_k = take_columns(P->A, &k, 1);
    bool unbounded = true;
    int r;
    double t = DBL_MAX;
    i = 0;
    for (i; i < A_k->rows; i++) {
        if (A_k->table[i][0] > 0) {
            unbounded = false;
            double temp = P->b->table[i][0] / A_k->table[i][0];
            if (temp < t) {
                t = temp;
                r = i;
            }
        }
    }
    if (unbounded) return UNBOUNDED;

    B[r] = k;
    return simplex(P, B);
}

void cast_to_SF(LP1 P) {
    assert(!P->is_SF);

    if (!P->obj_max) {
        P->z *= -1;
        multiply_scalar(P->c, -1);
        P->obj_max = true;
    }

    int ineqs = 0;
    int i = 0;
    for (i; i < P->num_constraints; i++) {
        if (P->constraint[i]->type != EQ) {
            ineqs++;
        }
    }

    P->A = create_matrix(P->num_constraints, P->variables + P->num_not_ge0 + ineqs);

    P->b = create_matrix(P->num_constraints, 1);

    Matrix new_c = create_matrix(P->variables + P->num_not_ge0 + ineqs, 1);

    i = 0;
    for (i; i < P->variables; i++) {

        new_c->table[i][0] = P->c->table[i][0];

        if (!P->ge0[i]) {
            i++;
            new_c->table[i][0] = -1 * P->c->table[i][0];
        }

        int j = 0;
        for (j; j < P->num_constraints; j++) {
            P->A->table[j][i] = P->constraint[j]->coefficients[i];
            if (!P->ge0[i])
                P->A->table[j][i] = -1 * P->constraint[j]->coefficients[i];
        }
    }

    destroy_matrix(P->c);
    P->c = new_c;

    int count_ineq = 0;
    i = 0;
    for (i; i < P->num_constraints; i++) {
        P->b->table[i][0] = P->constraint[i]->constant;
        if (P->constraint[i]->type != EQ) {
            bool le = P->constraint[i]->type == LE;
            P->A->table[i][P->variables + count_ineq] = le ? 1 : -1;
            count_ineq++;
        }
    }

    P->is_SF = true;
}

// P must be in canonical form for basis B
Matrix basic_solution(LP1 P, int *B) {
    Matrix x = create_matrix(P->A->cols, 1);
    int i = 0;
    for (i; i < P->A->rows; i++) {
        x->table[B[i]][0] = P->b->table[i][0];
    }
    return x;
}

Matrix solve(LP1 P) {
    assert(!P->is_SF);
    cast_to_SF(P);
    LP1 Q = copy_LP(P);
    Matrix im = create_identity(P->num_constraints);
    join_right(Q->A, im);
    destroy_matrix(im);

    destroy_matrix(Q->c);
    Q->c = create_matrix(Q->A->cols, 1);

    int *QB = malloc(sizeof(int) * P->A->rows);

    int i = 0;
    for (i; i < P->A->rows; i++) {
        Q->c->table[i + P->A->cols][0] = -1;
        QB[i] = P->A->cols + i;
    }


    simplex(Q, QB);
    Matrix x = basic_solution(Q, QB);
    bool feasible = true;
    i = 0;
    for (i; i < P->A->rows; i++) {
        if (x->table[i + P->A->cols][0] != 0)
            feasible = false;
    }
    destroy_matrix(x);
    if (!feasible) return NULL;

    simp final = simplex(P, QB);
    if (final == UNBOUNDED) return NULL;

    x = basic_solution(P, QB);
    Matrix rx = create_matrix(P->variables, 1);
    int v = 0;
    i = 0;
    for (v, i; v < P->variables; v++, i++) {
        if (P->ge0[v]) {
            rx->table[v][0] = x->table[i][0];
        }
        else {
            rx->table[v][0] = x->table[i][0] - x->table[i + 1][0];
            i++;
        }
    }
    destroy_matrix(x);
    return rx;
}

LP1 read_LP() {
    LP1 P = malloc(sizeof(struct lp));

    P->is_SF = false;
    P->z = 0;
    P->num_not_ge0 = 0;

    printf("Minimization (0) or Maximization (1)? : ");

    int maxmin;
    scanf("%d", &maxmin);
    if (maxmin < 1) P->obj_max = false;
    else P->obj_max = true;

    printf("\nNum of vars: ");
    scanf("%d", &P->variables);
    printf("\n");

    P->c = create_matrix(P->variables, 1);
    int i = 0;
    for (i; i < P->variables; i++) {
        printf("Var X%d: ", i + 1);
        scanf("%lf", &P->c->table[i][0]);
    }

    printf("\nNum of constraints: ");
    scanf("%d", &P->num_constraints);

    printf("\nEnter constraints:\n");
    P->constraint = malloc(sizeof(struct constraint *) * P->num_constraints);

    i = 0;
    for (i; i < P->num_constraints; i++) {
        P->constraint[i] = create_constraint(P->variables);

        printf("Constraint (%d) type is \"=\" or \">\" or \"<\"? : ", i + 1);
        char t[3];
        scanf("%s", t);

        if (t[0] == '=') P->constraint[i]->type = EQ;
        else if (t[0] == '>') P->constraint[i]->type = GE;
        else P->constraint[i]->type = LE;

        int v = 0;
        for (v; v < P->variables; v++) {
            printf("Row: %d -> Var: %d = ", i + 1, v + 1);
            scanf("%lf", &P->constraint[i]->coefficients[v]);
        }

        printf("Row: %d -> CONST. = ", i + 1);
        scanf("%lf", &P->constraint[i]->constant);
        printf("\n");
    }

    // read in which vars >= 0
    printf("Vars gte 0 from: %d\n", P->variables);

    P->ge0 = malloc(sizeof(bool) * P->variables);
    i = 0;

    for (i; i < P->variables; i++) {
//        int g0 = 1;
//        scanf("%d", &g0);
//        if (g0 == 0) {
//            P->ge0[i] = false;
//            P->num_not_ge0++;
//        }
//        else

        P->ge0[i] = true;
    }
    return P;
}

void print_LP(LP1 P) {
    if (P->obj_max)
        printf("max");
    else
        printf("min");
    printf(" %6.3f + [ ", P->z);
    int i = 0;
    for (i; i < P->c->rows; i++) {
        printf("%6.3f ", P->c->table[i][0]);
    }
    printf("]x\n");
    i = 0;
    for (i; i < P->A->rows; i++) {
        printf("[ ");
        int j = 0;
        for (j; j < P->A->cols; j++) {
            printf("%6.3f ", P->A->table[i][j]);
        }
        if (i == P->A->rows / 2)
            printf("] x = ");
        else
            printf("]     ");
        printf("[ %6.3f ]\n", P->b->table[i][0]);
    }
    printf("\n");
}
