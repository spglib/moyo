#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "moyoc.h"

void test_hall_symbol_entry(void) {
    MoyoHallSymbolEntry *entry = moyo_hall_symbol_entry_new(1);
    assert(entry != NULL);
    printf("hall_number=%d number=%d hall_symbol=%s hm_short=%s centering=%s\n",
           entry->hall_number, entry->number, entry->hall_symbol, entry->hm_short,
           entry->centering);
    assert(entry->hall_number == 1);
    assert(entry->number == 1);
    assert(entry->arithmetic_number == 1);
    assert(strcmp(entry->hall_symbol, "P 1") == 0);
    assert(strcmp(entry->hm_short, "P 1") == 0);
    assert(strcmp(entry->centering, "P") == 0);
    moyo_hall_symbol_entry_free(entry);

    // Ia-3d (no. 230)
    MoyoHallSymbolEntry *last_entry = moyo_hall_symbol_entry_new(530);
    assert(last_entry != NULL);
    assert(last_entry->number == 230);
    assert(strcmp(last_entry->centering, "I") == 0);
    moyo_hall_symbol_entry_free(last_entry);

    assert(moyo_hall_symbol_entry_new(0) == NULL);
    assert(moyo_hall_symbol_entry_new(INT32_MIN) == NULL);
    assert(moyo_hall_symbol_entry_new(531) == NULL);
}

void test_layer_hall_symbol_entry(void) {
    MoyoLayerHallSymbolEntry *entry = moyo_layer_hall_symbol_entry_new(116);
    assert(entry != NULL);
    printf("layer hall_number=%d number=%d hm_short=%s centering=%s\n",
           entry->hall_number, entry->number, entry->hm_short, entry->centering);
    assert(entry->hall_number == 116);
    assert(entry->number == 80);
    assert(strcmp(entry->centering, "P") == 0);
    moyo_layer_hall_symbol_entry_free(entry);

    assert(moyo_layer_hall_symbol_entry_new(0) == NULL);
    assert(moyo_layer_hall_symbol_entry_new(INT32_MIN) == NULL);
    assert(moyo_layer_hall_symbol_entry_new(117) == NULL);
}

void test_space_group_type(void) {
    // Ia-3d (no. 230)
    MoyoSpaceGroupType *space_group_type = moyo_space_group_type_new(230);
    assert(space_group_type != NULL);
    printf("number=%d hm_short=%s geometric_crystal_class=%s crystal_system=%s "
           "bravais_class=%s lattice_system=%s crystal_family=%s\n",
           space_group_type->number, space_group_type->hm_short,
           space_group_type->geometric_crystal_class, space_group_type->crystal_system,
           space_group_type->bravais_class, space_group_type->lattice_system,
           space_group_type->crystal_family);
    assert(space_group_type->number == 230);
    assert(strcmp(space_group_type->geometric_crystal_class, "m-3m") == 0);
    assert(strcmp(space_group_type->crystal_system, "Cubic") == 0);
    assert(strcmp(space_group_type->bravais_class, "cI") == 0);
    assert(strcmp(space_group_type->lattice_system, "Cubic") == 0);
    assert(strcmp(space_group_type->crystal_family, "Cubic") == 0);
    moyo_space_group_type_free(space_group_type);

    assert(moyo_space_group_type_new(0) == NULL);
    assert(moyo_space_group_type_new(INT32_MIN) == NULL);
    assert(moyo_space_group_type_new(231) == NULL);
}

void test_layer_group_type(void) {
    // LG 61 (p4/mmm)
    MoyoLayerGroupType *layer_group_type = moyo_layer_group_type_new(61);
    assert(layer_group_type != NULL);
    printf("number=%d hm_short=%s geometric_crystal_class=%s bravais_class=%s "
           "lattice_system=%s\n",
           layer_group_type->number, layer_group_type->hm_short,
           layer_group_type->geometric_crystal_class, layer_group_type->bravais_class,
           layer_group_type->lattice_system);
    assert(layer_group_type->number == 61);
    assert(strcmp(layer_group_type->geometric_crystal_class, "4/mmm") == 0);
    assert(strcmp(layer_group_type->bravais_class, "tp") == 0);
    assert(strcmp(layer_group_type->lattice_system, "Square") == 0);
    moyo_layer_group_type_free(layer_group_type);

    assert(moyo_layer_group_type_new(0) == NULL);
    assert(moyo_layer_group_type_new(INT32_MIN) == NULL);
    assert(moyo_layer_group_type_new(81) == NULL);
}

void test_magnetic_space_group_type(void) {
    // UNI 1242: R31'_c[R3] (BNS R_I3, 146.12)
    MoyoMagneticSpaceGroupType *magnetic_space_group_type =
        moyo_magnetic_space_group_type_new(1242);
    assert(magnetic_space_group_type != NULL);
    printf("uni_number=%d litvin_number=%d bns_number=%s og_number=%s number=%d "
           "construct_type=%d\n",
           magnetic_space_group_type->uni_number, magnetic_space_group_type->litvin_number,
           magnetic_space_group_type->bns_number, magnetic_space_group_type->og_number,
           magnetic_space_group_type->number, magnetic_space_group_type->construct_type);
    assert(magnetic_space_group_type->uni_number == 1242);
    assert(strcmp(magnetic_space_group_type->bns_number, "146.12") == 0);
    assert(magnetic_space_group_type->number == 146);
    assert(magnetic_space_group_type->construct_type == 4);
    moyo_magnetic_space_group_type_free(magnetic_space_group_type);

    assert(moyo_magnetic_space_group_type_new(0) == NULL);
    assert(moyo_magnetic_space_group_type_new(INT32_MIN) == NULL);
    assert(moyo_magnetic_space_group_type_new(1652) == NULL);
}

int main(void) {
    test_hall_symbol_entry();
    test_layer_hall_symbol_entry();
    test_space_group_type();
    test_layer_group_type();
    test_magnetic_space_group_type();
    return 0;
}
