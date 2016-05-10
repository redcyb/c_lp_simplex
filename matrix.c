#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

void destroy_data(double **d, int n);

double **create_data(int n, int m);

Matrix create_matrix(int n, int m) {
    Matrix mtr = malloc(sizeof(struct matrix));
    mtr->table = create_data(n, m);
    mtr->rows = n;
    mtr->cols = m;
    return mtr;
}

double **create_data(int n, int m) {
    double **d = malloc(sizeof(double *) * n);
    int i = 0;
    for (i; i < n; i++) {
        d[i] = malloc(sizeof(double) * m);
        int j = 0;
        for (j; j < m; j++) {
            d[i][j] = 0;
        }
    }
    return d;
}

void destroy_data(double **d, int n) {
    int i = 0;
    for (i; i < n; i++) {
        free(d[i]);
    }
    free(d);
}

void destroy_matrix(Matrix mtr) {
    assert(mtr != NULL);
    destroy_data(mtr->table, mtr->rows);
    free(mtr);
    mtr = NULL;
}

Matrix copy_matrix(Matrix mtr) {
    Matrix c = create_matrix(mtr->rows, mtr->cols);
    int i = 0;
    for (i; i < mtr->rows; i++) {
        int j = 0;
        for (j; j < mtr->cols; j++) {
            c->table[i][j] = mtr->table[i][j];
        }
    }
    return c;
}

void transpose(Matrix mtr) {
    double **nd = create_data(mtr->cols, mtr->rows);
    int i = 0;
    for (i; i < mtr->rows; i++) {
        int j = 0;
        for (j; j < mtr->cols; j++) {
            nd[j][i] = mtr->table[i][j];
        }
    }
    destroy_data(mtr->table, mtr->rows);
    mtr->table = nd;
    int temp = mtr->cols;
    mtr->cols = mtr->rows;
    mtr->rows = temp;
}

void mult_row(double *row, int m, double factor) {
    int i = 0;
    for (i; i < m; i++) {
        row[i] *= factor;
    }
}

void add_row(double *row1, double *row2, int m, double factor) {
    int i = 0;
    for (i; i < m; i++) {
        row1[i] += factor * row2[i];
    }
}

void swap_rows(double *row1, double *row2, int m) {
    double temp;
    int i = 0;
    for (i; i < m; i++) {
        temp = row1[i];
        row1[i] = row2[i];
        row2[i] = temp;
    }
}

int rref(Matrix mtr) {
    int cols = 0;
    int leading_ones = 0;
    int total_cols = mtr->cols < mtr->rows ? mtr->cols : mtr->rows;

    while (cols < total_cols) {
        if (mtr->table[leading_ones][cols] == 0) {
            int row2;
            for (row2 = leading_ones; row2 < mtr->rows; row2++) {
                if (mtr->table[row2][cols] != 0) break;
            }
            if (row2 < mtr->rows)
                swap_rows(mtr->table[leading_ones], mtr->table[row2], mtr->cols);
            else {
                cols++;
                continue;
            }
        }

        if (mtr->table[leading_ones][cols] != 1) {
            mult_row(mtr->table[leading_ones], mtr->cols,
                     1 / mtr->table[leading_ones][cols]);
        }

        int i = 0;
        for (i; i < mtr->rows; i++) {
            if (i != leading_ones && mtr->table[i][cols] != 0) {
                add_row(mtr->table[i], mtr->table[leading_ones], mtr->cols,
                        -mtr->table[i][cols]);
            }
        }
        cols++;
        leading_ones++;
    }
    return leading_ones;
}

void join_right(Matrix to, Matrix from) {
    assert(to->rows == from->rows);
    double **data = create_data(to->rows, to->cols + from->cols);
    int i = 0;
    for (i; i < to->rows; i++) {
        int j = 0;
        for (j; j < to->cols; j++) {
            data[i][j] = to->table[i][j];
        }
        j = 0;
        for (j; j < from->cols; j++) {
            data[i][to->cols + j] = from->table[i][j];
        }
    }
    destroy_data(to->table, to->rows);
    to->table = data;
    to->cols += from->cols;
}

void invert(Matrix mtr) {
    assert(mtr->cols == mtr->rows);
    Matrix im = create_identity(mtr->rows);
    join_right(mtr, im);
    destroy_matrix(im);
    rref(mtr);
    int *cols = malloc(sizeof(int) * mtr->rows);
    int i = 0;
    for (i; i < mtr->rows; i++) {
        cols[i] = i + mtr->rows;
    }
    Matrix inv = take_columns(mtr, cols, mtr->rows);
    mtr->table = inv->table;
    mtr->cols = inv->cols;
    free(cols);
    free(inv);
}

Matrix mult_new(const Matrix mtr1, const Matrix mtr2) {
    assert(mtr1->cols == mtr2->rows);
    Matrix mtr = malloc(sizeof(struct matrix));
    mtr->table = create_data(mtr1->rows, mtr2->cols);
    mtr->rows = mtr1->rows;
    mtr->cols = mtr2->cols;
    int i = 0;
    for (i; i < mtr1->rows; i++) {
        int j = 0;
        for (j; j < mtr2->cols; j++) {
            int p = 0; for (p; p < mtr1->cols; p++) {
                mtr->table[i][j] += mtr1->table[i][p] * mtr2->table[p][j];
            }
        }
    }
    return mtr;
}

void add_matrix(Matrix to, Matrix from) {
    assert(to->rows == from->rows);
    assert(to->cols == from->cols);
    int i = 0;
    for (i; i < to->rows; i++) {
        int j = 0;
        for (j; j < to->cols; j++) {
            to->table[i][j] += from->table[i][j];
        }
    }
}

void multiply_scalar(Matrix mtr, double s) {
    int i = 0;
    for (i; i < mtr->rows; i++) {
        int j = 0;
        for (j; j < mtr->cols; j++) {
            mtr->table[i][j] *= s;
        }
    }
}

Matrix create_identity(int n) {
    Matrix im = malloc(sizeof(struct matrix));
    im->table = create_data(n, n);
    im->rows = im->cols = n;
    int i = 0;
    for (i; i < n; i++) {
        im->table[i][i] = 1;
    }
    return im;
}

Matrix take_columns(Matrix mtr, int *cols, int m) {
    int i = 0;
    for (i; i < m; i++) {
        assert(cols[i] < mtr->cols);
    }
    Matrix res = malloc(sizeof(struct matrix));
    res->table = create_data(mtr->rows, m);
    res->rows = mtr->rows;
    res->cols = m;
    i = 0;
    for (i; i < mtr->rows; i++) {
        int j = 0;
        for (j; j < m; j++) {
            res->table[i][j] = mtr->table[i][cols[j]];
        }
    }
    return res;
}

Matrix take_rows(Matrix mtr, int *rows, int n) {
    Matrix res = malloc(sizeof(struct matrix));
    res->table = create_data(n, mtr->cols);
    res->rows = n;
    res->cols = mtr->cols;
    int i = 0;
    for (i; i < n; i++) {
        int j = 0;
        for (j; j < mtr->cols; j++) {
            res->table[i][j] = mtr->table[rows[i]][j];
        }
    }
    return res;
}

void print_matrix(Matrix mtr) {
    int i = 0;
    for (i; i < mtr->rows; i++) {
        printf("[ ");
        int j = 0;
        for (j; j < mtr->cols; j++) {
            printf("%6.3f ", mtr->table[i][j]);
        }
        printf("]\n");
    }
    printf("\n");
}

Matrix read_matrix() {
    int n, m;
    scanf("%d", &n);
    scanf("%d", &m);
    Matrix mtr = create_matrix(n, m);
    int i = 0;
    for (i; i < n; i++) {
        int j = 0;
        for (j; j < m; j++) {
            scanf("%lf", &mtr->table[i][j]);
        }
    }
    return mtr;
}

double single_to_num(Matrix mtr) {
    assert(mtr->cols == 1 && mtr->rows == 1);
    return mtr->table[0][0];
}
