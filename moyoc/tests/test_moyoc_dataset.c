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

    printf("dataset->number: %d\n", dataset->number);
    printf("dataset->hall_number: %d\n", dataset->hall_number);
    printf("dataset->hm_symbol: %s\n", dataset->hm_symbol);
    // for (int i = 0; i < dataset->operations.len; i++) {
    //     printf("operation %d\n", i);
    //     MoyocOperation_t operation = *(dataset->operations.ptr + i);
    //     for (int j = 0; j < 3; j++) {
    //         printf("  ");
    //         for (int k = 0; k < 3; k++) {
    //             printf("%d ", operation.rotation.idx[j].idx[k]);
    //         }
    //         printf("\n");
    //     }
    //     printf("  ");
    //     for (int j = 0; j < 3; j++) {
    //         printf("%f ", operation.translation.idx[j]);
    //     }
    //     printf("\n");
    // }

    assert(dataset != NULL);
    assert(dataset->number == 194);
    assert(dataset->hall_number == 488);
    assert(strcmp(dataset->hm_symbol, "P 6_3/m m c") == 0);
    // assert(dataset->operations.len == 24);

    free_moyo_dataset(dataset);

    return 0;
}
