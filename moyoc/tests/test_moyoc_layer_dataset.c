#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "moyoc.h"

// LG 61 (p4/mmm): square in-plane lattice with a single atom at (0, 0, z)
void test_p4_per_mmm(void) {
    double basis[3][3] = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 5.0},
    };
    double positions[][3] = {
        {0.0, 0.0, 0.1},
    };
    int32_t numbers[] = {1};
    int32_t num_atoms = 1;

    double symprec = 1e-4;
    double angle_tolerance = -1;
    MoyoLayerSetting setting = MOYO_LAYER_SETTING_STANDARD;
    int32_t hall_number = -1;
    bool rotate_basis = true;

    MoyoLayerDataset *dataset = moyo_layer_dataset_new(
        basis, positions, numbers, num_atoms,
        symprec, angle_tolerance, setting, hall_number, rotate_basis
    );
    assert(dataset != NULL);

    // Identification
    printf("dataset->number: %d\n", dataset->number);
    printf("dataset->hall_number: %d\n", dataset->hall_number);
    printf("dataset->hm_symbol: %s\n", dataset->hm_symbol);
    assert(dataset->number == 61);

    // Input cell
    assert(dataset->num_atoms == num_atoms);

    // Symmetry operations in the input cell
    printf("dataset->operations.num_operations: %d\n",
           dataset->operations.num_operations);
    assert(dataset->operations.num_operations == 16);

    // Site symmetry
    assert(dataset->orbits[0] == 0);
    printf("dataset->wyckoffs: %s\n", dataset->wyckoffs);
    assert(dataset->wyckoffs[0] == 'a');
    printf("dataset->site_symmetry_symbols[0]: %s\n",
           dataset->site_symmetry_symbols[0]);
    assert(strstr(dataset->site_symmetry_symbols[0], "4/mmm") != NULL);

    // Standardized layer cell: square with c along z
    assert(dataset->std_cell.num_atoms == 1);
    double a = 0.0;
    double b = 0.0;
    for (int i = 0; i < 3; i++) {
        a += dataset->std_cell.basis[0][i] * dataset->std_cell.basis[0][i];
        b += dataset->std_cell.basis[1][i] * dataset->std_cell.basis[1][i];
    }
    assert(fabs(sqrt(a) - sqrt(b)) < 1e-8);
    printf("dataset->pearson_symbol: %s\n", dataset->pearson_symbol);
    assert(strcmp(dataset->pearson_symbol, "tp1") == 0);

    // Primitive standardized layer cell
    assert(dataset->prim_std_cell.num_atoms == 1);
    assert(dataset->mapping_std_prim[0] == 0);

    moyo_layer_dataset_free(dataset);
}

// LG 2 (p-1): triclinic-oblique with two atoms at inversion-related positions
void test_p_minus_1(void) {
    double gamma = 80.0 * M_PI / 180.0;
    double a = 1.0;
    double b = 1.5;
    double basis[3][3] = {
        {a, 0.0, 0.0},
        {b * cos(gamma), b * sin(gamma), 0.0},
        {0.0, 0.0, 5.0},
    };
    double positions[][3] = {
        {0.1, 0.2, 0.1},
        {0.3, 0.5, 0.2},
    };
    int32_t numbers[] = {1, 1};
    int32_t num_atoms = 2;

    double symprec = 1e-4;
    double angle_tolerance = -1;
    MoyoLayerSetting setting = MOYO_LAYER_SETTING_STANDARD;
    int32_t hall_number = -1;
    bool rotate_basis = true;

    MoyoLayerDataset *dataset = moyo_layer_dataset_new(
        basis, positions, numbers, num_atoms,
        symprec, angle_tolerance, setting, hall_number, rotate_basis
    );
    assert(dataset != NULL);

    printf("dataset->number: %d\n", dataset->number);
    printf("dataset->hall_number: %d\n", dataset->hall_number);
    assert(dataset->number == 2);
    assert(dataset->hall_number == 2);
    assert(dataset->operations.num_operations == 2);
    assert(dataset->orbits[0] == 0);
    assert(dataset->orbits[1] == 0);
    printf("dataset->wyckoffs: %s\n", dataset->wyckoffs);
    assert(strcmp(dataset->wyckoffs, "ee") == 0);
    printf("dataset->pearson_symbol: %s\n", dataset->pearson_symbol);
    assert(strcmp(dataset->pearson_symbol, "mp2") == 0);

    moyo_layer_dataset_free(dataset);
}

int main(void) {
    test_p4_per_mmm();
    test_p_minus_1();
    return 0;
}
