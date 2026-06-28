mod arithmetic_crystal_class;
mod centering;
mod hall_symbol;
mod layer_arithmetic_crystal_class;
mod layer_centering;
mod layer_group_type;
mod layer_hall_symbol;
mod layer_setting;
mod magnetic_space_group_type;
mod setting;
mod space_group_type;
mod wyckoff;

pub use arithmetic_crystal_class::PyArithmeticCrystalClass;
pub use centering::PyCentering;
pub use hall_symbol::PyHallSymbolEntry;
pub use layer_arithmetic_crystal_class::PyLayerArithmeticCrystalClass;
pub use layer_centering::PyLayerCentering;
pub use layer_group_type::PyLayerGroupType;
pub use layer_hall_symbol::PyLayerHallSymbolEntry;
pub use layer_setting::PyLayerSetting;
pub use magnetic_space_group_type::PyMagneticSpaceGroupType;
pub use setting::PySetting;
pub use space_group_type::PySpaceGroupType;
pub use wyckoff::PyWyckoffPosition;

use pyo3::prelude::*;

use super::base::{PyMagneticOperations, PyMoyoError, PyOperations};
use moyo::data::{
    LayerNumber, LayerSetting, Number, UNINumber,
    magnetic_operations_from_uni_number as core_magnetic_operations_from_uni_number,
    operations_from_layer_number as core_operations_from_layer_number,
    operations_from_number as core_operations_from_number,
};

/// Return symmetry operations for the given space-group ITA ``number``.
///
/// Parameters
/// ----------
/// number : int
///     ITA number of the space group (1 - 230).
/// setting : Setting, optional
///     Setting of the space group. If ``None``, the default setting is used.
/// primitive : bool, optional
///     If ``True``, return operations of the primitive cell; otherwise return operations
///     of the conventional cell.
#[pyfunction]
#[pyo3(signature = (number, *, setting=None, primitive=false))]
pub fn operations_from_number(
    number: Number,
    setting: Option<PySetting>,
    primitive: bool,
) -> Result<PyOperations, PyMoyoError> {
    let setting = setting.map(|s| s.0).unwrap_or_default();
    let operations = core_operations_from_number(number, setting, primitive)?;
    Ok(PyOperations::from(operations))
}

/// Return symmetry operations for the given layer-group ``number``.
///
/// Parameters
/// ----------
/// number : int
///     Layer-group number (1 - 80).
/// setting : LayerSetting, optional
///     Setting of the layer group. If ``None``, the default setting is used.
/// primitive : bool, optional
///     If ``True``, return operations of the primitive cell; otherwise return operations
///     of the conventional cell.
#[pyfunction]
#[pyo3(signature = (number, *, setting=None, primitive=false))]
pub fn operations_from_layer_number(
    number: LayerNumber,
    setting: Option<PyLayerSetting>,
    primitive: bool,
) -> Result<PyOperations, PyMoyoError> {
    let setting = setting.map(LayerSetting::from).unwrap_or_default();
    let operations = core_operations_from_layer_number(number, setting, primitive)?;
    Ok(PyOperations::from(operations))
}

/// Return magnetic symmetry operations for the given UNI number.
///
/// Parameters
/// ----------
/// uni_number : int
///     UNI (and BNS) serial number of the magnetic space group.
/// primitive : bool, optional
///     If ``True``, return magnetic operations of the primitive cell; otherwise return
///     operations of the conventional cell.
#[pyfunction]
#[pyo3(signature = (uni_number, *, primitive=false))]
pub fn magnetic_operations_from_uni_number(
    uni_number: UNINumber,
    primitive: bool,
) -> Result<PyMagneticOperations, PyMoyoError> {
    let magnetic_operations = core_magnetic_operations_from_uni_number(uni_number, primitive)?;
    Ok(PyMagneticOperations::from(magnetic_operations))
}
