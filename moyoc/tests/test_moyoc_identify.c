#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "moyoc.h"

void test_point_group(void) {
    // P3c1 (no. 158) -> arithmetic crystal class 45 (3m1P)
    MoyoOperations *operations =
        moyo_operations_from_number(158, MOYO_SETTING_STANDARD, -1, true);
    assert(operations != NULL);

    MoyoPointGroup *point_group =
        moyo_point_group_new(operations->rotations, operations->num_operations, NULL);
    assert(point_group != NULL);
    printf("arithmetic_number (no. 158): %d\n", point_group->arithmetic_number);
    assert(point_group->arithmetic_number == 45);
    moyo_point_group_free(point_group);

    // Identity basis must give the same result as no basis
    const double identity_basis[3][3] = {
        {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    MoyoPointGroup *point_group_with_basis = moyo_point_group_new(
        operations->rotations, operations->num_operations, identity_basis);
    assert(point_group_with_basis != NULL);
    assert(point_group_with_basis->arithmetic_number == 45);
    moyo_point_group_free(point_group_with_basis);

    moyo_operations_free(operations);

    // Invalid arguments return NULL
    const int32_t identity_rotation[1][3][3] = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    assert(moyo_point_group_new(NULL, 1, NULL) == NULL);
    assert(moyo_point_group_new(identity_rotation, 0, NULL) == NULL);
}

void test_space_group(void) {
    // P31c (no. 159)
    MoyoOperations *operations =
        moyo_operations_from_number(159, MOYO_SETTING_STANDARD, -1, true);
    assert(operations != NULL);

    MoyoSpaceGroup *space_group = moyo_space_group_new(
        operations->rotations, operations->translations, operations->num_operations,
        NULL, MOYO_SETTING_STANDARD, -1, -1.0);
    assert(space_group != NULL);
    printf("number: %d, hall_number: %d\n", space_group->number,
           space_group->hall_number);
    assert(space_group->number == 159);
    assert(1 <= space_group->hall_number && space_group->hall_number <= 530);
    moyo_space_group_free(space_group);

    // Invalid arguments return NULL
    assert(moyo_space_group_new(NULL, operations->translations,
                                operations->num_operations, NULL,
                                MOYO_SETTING_STANDARD, -1, -1.0) == NULL);
    assert(moyo_space_group_new(operations->rotations, operations->translations, 0,
                                NULL, MOYO_SETTING_STANDARD, -1, -1.0) == NULL);
    assert(moyo_space_group_new(operations->rotations, operations->translations,
                                operations->num_operations, NULL,
                                MOYO_SETTING_HALL_NUMBER, 0, -1.0) == NULL);

    moyo_operations_free(operations);
}

void test_layer_group(void) {
    // p6 (LG 65)
    MoyoOperations *operations =
        moyo_operations_from_layer_number(65, MOYO_LAYER_SETTING_STANDARD, -1, true);
    assert(operations != NULL);

    MoyoLayerGroup *layer_group = moyo_layer_group_new(
        operations->rotations, operations->translations, operations->num_operations,
        NULL, MOYO_LAYER_SETTING_STANDARD, -1, -1.0);
    assert(layer_group != NULL);
    printf("layer number: %d, layer hall_number: %d\n", layer_group->number,
           layer_group->hall_number);
    assert(layer_group->number == 65);
    assert(1 <= layer_group->hall_number && layer_group->hall_number <= 116);
    moyo_layer_group_free(layer_group);

    // Identity basis exercises the from_lattice path
    const double identity_basis[3][3] = {
        {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    MoyoLayerGroup *layer_group_with_basis = moyo_layer_group_new(
        operations->rotations, operations->translations, operations->num_operations,
        identity_basis, MOYO_LAYER_SETTING_STANDARD, -1, -1.0);
    assert(layer_group_with_basis != NULL);
    assert(layer_group_with_basis->number == 65);
    moyo_layer_group_free(layer_group_with_basis);

    moyo_operations_free(operations);
}

void test_magnetic_space_group(void) {
    // R31'_c (UNI 1242)
    MoyoMagneticOperations *magnetic_operations =
        moyo_magnetic_operations_from_uni_number(1242, true);
    assert(magnetic_operations != NULL);

    MoyoMagneticSpaceGroup *magnetic_space_group = moyo_magnetic_space_group_new(
        magnetic_operations->rotations, magnetic_operations->translations,
        magnetic_operations->time_reversals, magnetic_operations->num_operations,
        NULL, -1.0);
    assert(magnetic_space_group != NULL);
    printf("uni_number: %d\n", magnetic_space_group->uni_number);
    assert(magnetic_space_group->uni_number == 1242);
    moyo_magnetic_space_group_free(magnetic_space_group);

    // Invalid arguments return NULL
    assert(moyo_magnetic_space_group_new(
               magnetic_operations->rotations, magnetic_operations->translations,
               NULL, magnetic_operations->num_operations, NULL, -1.0) == NULL);

    moyo_magnetic_operations_free(magnetic_operations);
}

int main(void) {
    test_point_group();
    test_space_group();
    test_layer_group();
    test_magnetic_space_group();
    return 0;
}
