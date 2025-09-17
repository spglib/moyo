#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "moyoc.h"

int main() {
    // Wurtzite (P6_3mc, No. 186)
    double a = 3.81;
    double c = 6.24;
    double basis[3][3] = {
        {a, 0.0, 0.0},
        {-a / 2.0, a * sqrt(3.0) / 2.0, 0.0},
        {0.0, 0.0, c},
    };
    double z1_2b = 0.00014;
    double z2_2b = 0.37486;
    double positions[][3] = {
        {1.0 / 3.0, 2.0 / 3.0, z1_2b},
        {2.0 / 3.0, 1.0 / 3.0, z1_2b + 0.5},
        {1.0 / 3.0, 2.0 / 3.0, z2_2b},
        {2.0 / 3.0, 1.0 / 3.0, z2_2b + 0.5},
    };
    int numbers[] = {0, 0, 1, 1};
    int num_atoms = 4;

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
    printf("dataset->hm_symbol: %s\n", dataset->hm_symbol);
    assert(dataset->number == 186);
    assert(strcmp(dataset->hm_symbol, "P 6_3 m c") == 0);

    free_moyo_dataset(dataset);
}
