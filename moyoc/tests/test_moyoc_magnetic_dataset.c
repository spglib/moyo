#include <assert.h>
#include <stdio.h>

#include "moyoc.h"

void show_magnetic_operations(const MoyoMagneticOperations *magnetic_operations) {
    for (int i = 0; i < magnetic_operations->num_operations; i++) {
        printf("Magnetic operation %d (time reversal: %d)\n", i,
               magnetic_operations->time_reversals[i]);
        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                printf("%2d ", magnetic_operations->rotations[i][a][b]);
            }
            printf("\n");
        }
        for (int a = 0; a < 3; a++) {
            printf("%.2f ", magnetic_operations->translations[i][a]);
        }
        printf("\n");
    }
}

void test_collinear(void) {
    // Rutile-type MnF2, type-III magnetic space group
    double basis[3][3] = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
    };
    double positions[][3] = {
        // Ti (2a)
        {0.0, 0.0, 0.0},
        {0.5, 0.5, 0.5},
        // O (4f)
        {0.3, 0.3, 0.0},
        {0.7, 0.7, 0.0},
        {0.2, 0.8, 0.5},
        {0.8, 0.2, 0.5},
    };
    int32_t numbers[] = {0, 0, 1, 1, 1, 1};
    double magnetic_moments[] = {0.7, -0.7, 0.0, 0.0, 0.0, 0.0};
    int32_t num_atoms = 6;

    double symprec = 1e-4;
    double angle_tolerance = -1;
    double mag_symprec = -1;
    bool is_axial = false;
    bool rotate_basis = true;

    MoyoCollinearMagneticDataset *dataset = moyo_collinear_magnetic_dataset_new(
        basis, positions, numbers, magnetic_moments, num_atoms,
        symprec, angle_tolerance, mag_symprec, is_axial, rotate_basis
    );
    assert(dataset != NULL);

    // Magnetic space-group type
    printf("dataset->uni_number: %d\n", dataset->uni_number);
    assert(dataset->uni_number == 1158);  // BNS 136.498

    // Input magnetic cell
    assert(dataset->num_atoms == num_atoms);

    // Magnetic symmetry operations in the input cell
    printf("dataset->magnetic_operations:\n");
    show_magnetic_operations(&dataset->magnetic_operations);
    assert(dataset->magnetic_operations.num_operations == 16);

    // Site symmetry
    printf("dataset->orbits:\n");
    for (int i = 0; i < num_atoms; i++) {
        printf("%d ", dataset->orbits[i]);
    }
    printf("\n");
    int32_t expected_orbits[] = {0, 0, 2, 2, 2, 2};
    for (int i = 0; i < num_atoms; i++) {
        assert(dataset->orbits[i] == expected_orbits[i]);
    }

    // Standardized magnetic cell
    assert(dataset->std_mag_cell.cell.num_atoms == 6);
    printf("dataset->std_mag_cell.magnetic_moments:\n");
    for (int i = 0; i < dataset->std_mag_cell.cell.num_atoms; i++) {
        printf("%f ", dataset->std_mag_cell.magnetic_moments[i]);
    }
    printf("\n");

    // Primitive standardized magnetic cell
    assert(dataset->prim_std_mag_cell.cell.num_atoms == 6);
    printf("dataset->mapping_std_prim:\n");
    for (int i = 0; i < num_atoms; i++) {
        printf("%d ", dataset->mapping_std_prim[i]);
        assert(dataset->mapping_std_prim[i] == i);
    }
    printf("\n");

    // Final parameters
    printf("dataset->symprec: %f\n", dataset->symprec);
    printf("dataset->angle_tolerance: %f\n", dataset->angle_tolerance);
    printf("dataset->mag_symprec: %f\n", dataset->mag_symprec);

    moyo_collinear_magnetic_dataset_free(dataset);
}

void test_noncollinear(void) {
    // Rutile-type structure with non-collinear moments along c
    double basis[3][3] = {
        {4.603, 0.0, 0.0},
        {0.0, 4.603, 0.0},
        {0.0, 0.0, 2.969},
    };
    double x_4f = 0.3046;
    double positions[][3] = {
        // Ti (2a)
        {0.0, 0.0, 0.0},
        {0.5, 0.5, 0.5},
        // O (4f)
        {x_4f, x_4f, 0.0},
        {1.0 - x_4f, 1.0 - x_4f, 0.0},
        {0.5 - x_4f, 0.5 + x_4f, 0.5},
        {0.5 + x_4f, 0.5 - x_4f, 0.5},
    };
    int32_t numbers[] = {0, 0, 1, 1, 1, 1};
    double magnetic_moments[][3] = {
        {0.0, 0.0, 0.7},
        {0.0, 0.0, -0.7},
        {0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0},
    };
    int32_t num_atoms = 6;

    double symprec = 1e-4;
    double angle_tolerance = -1;
    double mag_symprec = -1;
    bool is_axial = true;
    bool rotate_basis = true;

    MoyoNonCollinearMagneticDataset *dataset = moyo_noncollinear_magnetic_dataset_new(
        basis, positions, numbers, magnetic_moments, num_atoms,
        symprec, angle_tolerance, mag_symprec, is_axial, rotate_basis
    );
    assert(dataset != NULL);

    // Magnetic space-group type
    printf("dataset->uni_number: %d\n", dataset->uni_number);
    assert(dataset->uni_number == 1159);  // BNS 136.499

    // Input magnetic cell
    assert(dataset->num_atoms == num_atoms);

    // Standardized magnetic cell
    assert(dataset->std_mag_cell.cell.num_atoms == 6);
    printf("dataset->std_mag_cell.magnetic_moments:\n");
    for (int i = 0; i < dataset->std_mag_cell.cell.num_atoms; i++) {
        printf("%f %f %f\n", dataset->std_mag_cell.magnetic_moments[i][0],
               dataset->std_mag_cell.magnetic_moments[i][1],
               dataset->std_mag_cell.magnetic_moments[i][2]);
    }

    // Primitive standardized magnetic cell
    assert(dataset->prim_std_mag_cell.cell.num_atoms == 6);

    moyo_noncollinear_magnetic_dataset_free(dataset);
}

int main(void) {
    test_collinear();
    test_noncollinear();
    return 0;
}
