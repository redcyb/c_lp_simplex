#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include "linear_program.h"

int main() {
    LP P = read_LP();
    Matrix mtr = solve(P);

    printf("\n=== RESULT ===\n\n");

    if(mtr != NULL) {
        print_matrix(mtr);
        destroy_matrix(mtr);
    }
    else
        printf("\nAnswer is NULL\n");

    printf("=== Press any key to exit the program ===\n");
    getch();
    destroy_LP(P);

    return 0;
}
