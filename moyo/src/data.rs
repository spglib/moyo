mod arithmetic_crystal_class;
mod centering;
mod classification;
mod hall_symbol;
mod hall_symbol_database;
mod layer;
mod magnetic;
mod operations;
mod point_group;
mod setting;
mod wyckoff;

pub use arithmetic_crystal_class::{
    ArithmeticCrystalClassEntry, ArithmeticNumber, arithmetic_crystal_class_entry,
};
pub use centering::Centering;
pub use classification::{
    BravaisClass, CrystalFamily, CrystalSystem, GeometricCrystalClass, LatticeSystem,
};
pub use hall_symbol::{HallSymbol, LayerHallSymbol, MagneticHallSymbol};
pub use hall_symbol_database::{HallNumber, HallSymbolEntry, Number, hall_symbol_entry};
pub use layer::{
    LayerArithmeticCrystalClassEntry, LayerArithmeticNumber, LayerBravaisClass, LayerCentering,
    LayerCrystalSystem, LayerHallNumber, LayerHallSymbolEntry, LayerLatticeSystem, LayerNumber,
    LayerSetting, LayerWyckoffPosition, iter_layer_arithmetic_crystal_entry,
    iter_layer_hall_symbol_entry, iter_layer_wyckoff_positions_from_hall_number,
    layer_arithmetic_crystal_class_entry, layer_hall_symbol_entry,
};
pub use magnetic::{
    ConstructType, MagneticHallSymbolEntry, MagneticSpaceGroupType, NUM_MAGNETIC_SPACE_GROUP_TYPES,
    UNINumber, get_magnetic_space_group_type, magnetic_hall_symbol_entry,
};
pub use operations::{
    magnetic_operations_from_uni_number, operations_from_layer_number, operations_from_number,
};
pub use setting::Setting;

pub(super) use arithmetic_crystal_class::iter_arithmetic_crystal_entry;
pub(super) use layer::{LayerPointGroupRepresentative, iter_layer_wyckoff_positions};
pub(super) use magnetic::uni_number_range;
pub(super) use point_group::PointGroupRepresentative;
pub use wyckoff::{WyckoffPosition, iter_wyckoff_positions_from_hall_number};

pub(super) use wyckoff::{WyckoffPositionSpace, iter_wyckoff_positions};
