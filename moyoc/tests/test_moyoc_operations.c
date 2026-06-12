#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#include "moyoc.h"

void test_operations_from_number(void) {
    // Ia-3d (no. 230): 96 operations in the conventional cell, 48 in the primitive cell
    MoyoOperations *operations =
        moyo_operations_from_number(230, MOYO_SETTING_SPGLIB, -1, false);
    assert(operations != NULL);
    printf("num_operations (no. 230, conventional): %d\n", operations->num_operations);
    assert(operations->num_operations == 96);
    moyo_operations_free(operations);

    MoyoOperations *prim_operations =
        moyo_operations_from_number(230, MOYO_SETTING_SPGLIB, -1, true);
    assert(prim_operations != NULL);
    printf("num_operations (no. 230, primitive): %d\n", prim_operations->num_operations);
    assert(prim_operations->num_operations == 48);
    moyo_operations_free(prim_operations);

    // Specific Hall number: hall_number 530 is no. 230
    MoyoOperations *hall_operations =
        moyo_operations_from_number(-1, MOYO_SETTING_HALL_NUMBER, 530, false);
    assert(hall_operations != NULL);
    assert(hall_operations->num_operations == 96);
    moyo_operations_free(hall_operations);

    // Invalid arguments return NULL
    assert(moyo_operations_from_number(231, MOYO_SETTING_SPGLIB, -1, false) == NULL);
    assert(moyo_operations_from_number(0, MOYO_SETTING_SPGLIB, -1, false) == NULL);
    assert(moyo_operations_from_number(INT32_MIN, MOYO_SETTING_SPGLIB, -1, false) == NULL);
    assert(moyo_operations_from_number(-1, MOYO_SETTING_HALL_NUMBER, -1, false) == NULL);
    assert(moyo_operations_from_number(-1, MOYO_SETTING_HALL_NUMBER, INT32_MIN, false) == NULL);
}

void test_operations_from_layer_number(void) {
    // LG 72 (p-3m1): order 12 with a primitive lattice
    MoyoOperations *operations =
        moyo_operations_from_layer_number(72, MOYO_LAYER_SETTING_STANDARD, -1, false);
    assert(operations != NULL);
    printf("num_operations (LG 72, conventional): %d\n", operations->num_operations);
    assert(operations->num_operations == 12);
    moyo_operations_free(operations);

    MoyoOperations *prim_operations =
        moyo_operations_from_layer_number(72, MOYO_LAYER_SETTING_STANDARD, -1, true);
    assert(prim_operations != NULL);
    printf("num_operations (LG 72, primitive): %d\n", prim_operations->num_operations);
    assert(prim_operations->num_operations == 12);
    moyo_operations_free(prim_operations);

    // Invalid arguments return NULL
    assert(moyo_operations_from_layer_number(81, MOYO_LAYER_SETTING_STANDARD, -1, false) == NULL);
    assert(moyo_operations_from_layer_number(0, MOYO_LAYER_SETTING_STANDARD, -1, false) == NULL);
    assert(moyo_operations_from_layer_number(INT32_MIN, MOYO_LAYER_SETTING_STANDARD, -1, false) ==
           NULL);
}

void test_magnetic_operations_from_uni_number(void) {
    // UNI 1242: R31'_c[R3] (BNS R_I3, 146.12)
    MoyoMagneticOperations *magnetic_operations =
        moyo_magnetic_operations_from_uni_number(1242, false);
    assert(magnetic_operations != NULL);
    printf("num_operations (UNI 1242, conventional): %d\n",
           magnetic_operations->num_operations);
    assert(magnetic_operations->num_operations == 18);
    int num_time_reversals = 0;
    for (int i = 0; i < magnetic_operations->num_operations; i++) {
        if (magnetic_operations->time_reversals[i]) {
            num_time_reversals++;
        }
    }
    printf("num_time_reversals (UNI 1242, conventional): %d\n", num_time_reversals);
    assert(num_time_reversals == 9);
    moyo_magnetic_operations_free(magnetic_operations);

    MoyoMagneticOperations *prim_magnetic_operations =
        moyo_magnetic_operations_from_uni_number(1242, true);
    assert(prim_magnetic_operations != NULL);
    printf("num_operations (UNI 1242, primitive): %d\n",
           prim_magnetic_operations->num_operations);
    assert(prim_magnetic_operations->num_operations == 6);
    moyo_magnetic_operations_free(prim_magnetic_operations);

    // Invalid arguments return NULL
    assert(moyo_magnetic_operations_from_uni_number(1652, false) == NULL);
    assert(moyo_magnetic_operations_from_uni_number(0, false) == NULL);
    assert(moyo_magnetic_operations_from_uni_number(INT32_MIN, false) == NULL);
}

int main(void) {
    test_operations_from_number();
    test_operations_from_layer_number();
    test_magnetic_operations_from_uni_number();
    return 0;
}
