#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "moyoc.h"

void show_cell(const MoyoCell *cell) {
    printf("Basis:\n");
    printf("a: %f %f %f\n", cell->basis[0][0], cell->basis[0][1], cell->basis[0][2]);
    printf("b: %f %f %f\n", cell->basis[1][0], cell->basis[1][1], cell->basis[1][2]);
    printf("c: %f %f %f\n", cell->basis[2][0], cell->basis[2][1], cell->basis[2][2]);

    printf("Positions:\n");
    for (int i = 0; i < cell->num_atoms; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%f ", cell->positions[i][j]);
        }
        printf("\n");
    }

    printf("Atomic numbers:\n");
    for (int i = 0; i < cell->num_atoms; i++) {
        printf("%d ", cell->numbers[i]);
    }
    printf("\n");
}

void show_operations(const MoyoOperations *operations) {
    for (int i = 0; i < (int)operations->num_operations; i++) {
        printf("Operation %d\n", i);
        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                printf("%2d ", operations->rotations[i][a][b]);
            }
            printf("\n");
        }
        for (int a = 0; a < 3; a++) {
            printf("%.2f ", operations->translations[i][a]);
        }
        printf("\n");
    }
}

void show_matrix3f(const double matrix[3][3]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
}

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
    printf("dataset->operations:\n");
    show_operations(&dataset->operations);
    assert(dataset->operations.num_operations == 24);

    // Site symmetry
    printf("dataset->orbits:\n");
    for (int i = 0; i < num_atoms; i++) {
        printf("%lu ", dataset->orbits[i]);
    }
    printf("\n");
    assert(dataset->orbits[0] == 0);
    assert(dataset->orbits[1] == 0);
    printf("dataset->wyckoffs:\n");
    for (int i = 0; i < num_atoms; i++) {
        printf("%c ", dataset->wyckoffs[i]);
    }
    printf("\n");
    assert(dataset->wyckoffs[0] == 'c');
    assert(dataset->wyckoffs[1] == 'c');
    printf("dataset->site_symmetry_symbols:\n");
    for (int i = 0; i < num_atoms; i++) {
        printf("%s ", dataset->site_symmetry_symbols[i]);
    }
    printf("\n");
    assert(strcmp(dataset->site_symmetry_symbols[0], "-6m2") == 0);
    assert(strcmp(dataset->site_symmetry_symbols[1], "-6m2") == 0);

    // Standardized cell
    printf("dataset->std_cell:\n");
    show_cell(&dataset->std_cell);
    assert(dataset->std_cell.num_atoms == 2);
    printf("dataset->std_linear:\n");
    show_matrix3f(dataset->std_linear);
    printf("dataset->std_origin_shift:\n");
    for (int i = 0; i < 3; i++) {
        printf("%f ", dataset->std_origin_shift[i]);
    }
    printf("\n");
    printf("dataset->std_rotation_matrix:\n");
    show_matrix3f(dataset->std_rotation_matrix);
    printf("dataset->pearson_symbol: %s\n", dataset->pearson_symbol);
    assert(strcmp(dataset->pearson_symbol, "hP2") == 0);

    // Primitive standardized cell
    printf("dataset->prim_std_cell:\n");
    show_cell(&dataset->prim_std_cell);
    printf("dataset->prim_std_linear:\n");
    show_matrix3f(dataset->prim_std_linear);
    printf("dataset->prim_std_origin_shift:\n");
    for (int i = 0; i < 3; i++) {
        printf("%f ", dataset->prim_std_origin_shift[i]);
    }
    printf("\n");
    printf("dataset->mapping_std_prim:\n");
    for (int i = 0; i < num_atoms; i++) {
        printf("%lu ", dataset->mapping_std_prim[i]);
    }
    printf("\n");
    assert(dataset->mapping_std_prim[0] == 0);
    assert(dataset->mapping_std_prim[1] == 1);

    // Final parameters
    printf("dataset->symprec: %f\n", dataset->symprec);
    printf("dataset->angle_tolerance: %f\n", dataset->angle_tolerance);

    free_moyo_dataset(dataset);

    return 0;
}
