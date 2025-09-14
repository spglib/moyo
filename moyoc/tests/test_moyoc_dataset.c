#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "moyoc.h"

int main(void) {
    // hcp
    double a = 3.17;
    double c = 5.14;
    double basis[3][3] = {
        {a, 0.0, 0.0},
        {-a / 2.0, a * sqrt(3.0) / 2.0, 0.0},
        {0.0, 0.0, c},
    };
    double positions[][3] = {
        {1.0 / 3.0, 2.0 / 3.0, 1.0 / 4.0},
        {2.0 / 3.0, 1.0 / 3.0, 3.0 / 4.0},
    };
    int numbers[] = {0, 0};
    int num_atoms = 2;

    double symprec = 1e-4;
    double angle_tolerance = -1;
    MoyoSetting setting = MOYO_SETTING_SPGLIB;
    int hall_number = -1;

    MoyoDataset *dataset = moyo_dataset(
        &basis, positions, numbers, num_atoms,
        symprec, angle_tolerance, setting, hall_number
    );
    assert(dataset != NULL);

    // Identification
    printf("dataset->number: %d\n", dataset->number);
    printf("dataset->hall_number: %d\n", dataset->hall_number);
    printf("dataset->hm_symbol: %s\n", dataset->hm_symbol);
    assert(dataset->number == 194);
    assert(dataset->hall_number == 488);
    assert(strcmp(dataset->hm_symbol, "P 6_3/m m c") == 0);

    // Symmetry operations in the input cell
    for (int i = 0; i < (int)dataset->operations.num_operations; i++) {
        printf("Operation %d\n", i);
        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                printf("%2d ", dataset->operations.rotations[i][a][b]);
            }
            printf("\n");
        }
        for (int a = 0; a < 3; a++) {
            printf("%.2f ", dataset->operations.translations[i][a]);
        }
        printf("\n");
    }
    assert(dataset->operations.num_operations == 24);

    free_moyo_dataset(dataset);

    return 0;
}
