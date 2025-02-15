#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "moyoc.h"

int main(void) {
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
    MoyocSetting setting = STANDARD;

    const MoyocDataset *dataset = moyoc_dataset(&basis, positions, numbers, num_atoms, symprec, angle_tolerance, setting);
    printf("dataset->number: %d\n", dataset->number);
    assert(dataset != NULL);
    assert(dataset->number == 194);

    // TODO: free dataset

    return 0;
}
