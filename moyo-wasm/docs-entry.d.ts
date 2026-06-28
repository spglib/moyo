// Public API surface for the generated TypeDoc reference (see typedoc.json).
//
// wasm-pack emits `pkg/moyo_wasm.d.ts` with the full crystallographic API plus
// the low-level wasm-bindgen init plumbing (`InitInput`, `InitOutput`,
// `SyncInitInput` and their `__wbindgen_*` internals). Re-exporting an explicit
// allowlist here keeps the reference focused on the public API. `init` /
// `initSync` are kept because callers need them to load the module.
//
// When a new `#[wasm_bindgen]` function or `tsify` type becomes public, add it
// here so it appears in the reference.
export {
  default as init,
  initSync,
  // Analysis
  analyze_cell,
  // Space groups
  space_group_type,
  operations_from_number,
  hall_symbol_entry,
  hall_symbol_entries_from_number,
  arithmetic_crystal_class,
  wyckoff_positions,
  // Layer groups
  layer_group_type,
  operations_from_layer_number,
  layer_hall_symbol_entry,
  layer_hall_symbol_entries_from_number,
  layer_arithmetic_crystal_class,
  layer_wyckoff_positions,
  // Magnetic space groups
  magnetic_space_group_type,
  magnetic_operations_from_uni_number,
  magnetic_hall_symbol_entry,
  // Types
  type Lattice,
  type AngleTolerance,
  type MoyoSetting,
  type MoyoLayerSetting,
  type MoyoCentering,
  type MoyoLayerCentering,
  type MoyoCell,
  type MoyoOperation,
  type MoyoMagneticOperation,
  type MoyoDataset,
  type MoyoSpaceGroupType,
  type MoyoArithmeticCrystalClass,
  type MoyoHallSymbolEntry,
  type MoyoWyckoffPosition,
  type MoyoLayerGroupType,
  type MoyoLayerArithmeticCrystalClass,
  type MoyoLayerHallSymbolEntry,
  type MoyoLayerWyckoffPosition,
  type MoyoMagneticSpaceGroupType,
  type MoyoMagneticHallSymbolEntry,
} from "./pkg/moyo_wasm";
