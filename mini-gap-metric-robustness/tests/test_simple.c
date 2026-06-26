#include <stdio.h>
#include "gap_core.h"
int main(void) {
    printf("Starting simple test...\n");
    fflush(stdout);
    GapMatrix* m = gap_mat_eye(2);
    if (m) {
        printf("Matrix created: %dx%d\n", m->rows, m->cols);
        gap_mat_free(m);
    }
    printf("Simple test done.\n");
    return 0;
}
