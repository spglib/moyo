use std::collections::hash_map::Entry;
use std::collections::{HashMap, VecDeque};

use log::debug;
use nalgebra::{matrix, Vector3};

use super::centering::Centering;
use super::hall_symbol_database::{hall_symbol_entry, HallNumber};
use super::magnetic_hall_symbol_database::magnetic_hall_symbol_entry;
use super::magnetic_space_group::UNINumber;
use crate::base::{
    MagneticOperation, MagneticOperations, Operation, Operations, OriginShift, Rotation,
    TimeReversal, Transformation, Translation, EPS,
};

const MAX_DENOMINATOR: i32 = 12;

// Hall symbol. See A1.4.2.3 in ITB (2010).
//
// Extended Backus-Naur form (EBNF) for Hall symbols
// ----------------------------------------------------------
// <Hall symbol>    := <L> <N>+ <V>?
// <L>              := "-"? <lattice symbol>
// <lattice symbol> := [PABCIRHF]  # Table A1.4.2.2
// <N>              := <nfold> <A>? <T>?
// <nfold>          := "-"? ("1" | "2" | "3" | "4" | "6")
// <A>              := [xyz] | "^" | "=" | "*"  # Table A.1.4.2.4, Table A.1.4.2.5, Table A.1.4.2.6
//                                              # The original axis symbol "'" is substituted with "^"
//                                              # The original axis symbol """ is substituted with "="
// <T>              := [abcnuvwd] | [1-6]  # Table A.1.4.2.3
// <V>              := "(" [0-11] [0-11] [0-11] ")"
#[derive(Debug)]
pub struct HallSymbol {
    pub hall_symbol: String,
    pub centering: Centering,
    pub centering_translations: Vec<Translation>,
    /// Generators of the space group except pure translations
    pub generators: Operations,
}

impl HallSymbol {
    pub fn new(hall_symbol: &str) -> Option<Self> {
        let tokens = tokenize(hall_symbol);
        let (inversion_at_origin, lattice_symbol, ns, origin_shift) = parse(&tokens)?;

        // From translation subgroup
        let centering_translations = get_centering_translations(&lattice_symbol);

        // change basis by (I, v)
        // (R, tau) -> (R, tau + v - Rv)
        let mut generators = vec![];

        if inversion_at_origin {
            generators.push(Operation::new(-Rotation::identity(), 2.0 * origin_shift));
        }

        for (rotation, translation, time_reversal) in ns {
            if time_reversal {
                debug!("Use MagneticHallSymbol for magnetic space groups.");
                return None;
            }
            let translation_mod1 = (translation + origin_shift
                - rotation.map(|e| e as f64) * origin_shift)
                .map(|e| e.rem_euclid(1.0));
            generators.push(Operation::new(rotation, translation_mod1));
        }

        Some(Self {
            hall_symbol: hall_symbol.to_string(),
            centering: lattice_symbol,
            centering_translations,
            generators,
        })
    }

    /// Traverse all the symmetry operations up to translations by conventional cell.
    /// The order of operations is guaranteed to be fixed.
    pub fn traverse(&self) -> Operations {
        let mut queue = VecDeque::new();
        let mut hm_translations = HashMap::new();
        let mut operations = vec![];

        queue.push_back(Operation::identity());

        while !queue.is_empty() {
            let ops_lhs = queue.pop_front().unwrap();
            let entry = hm_translations.entry(ops_lhs.rotation);
            if let Entry::Occupied(_) = entry {
                continue;
            }
            entry.or_insert(ops_lhs.translation);
            operations.push(ops_lhs.clone());

            for rhs in self.generators.iter() {
                let new_ops = ops_lhs.clone() * rhs.clone();
                if !hm_translations.contains_key(&new_ops.rotation) {
                    let new_translation_mod1 = purify_translation_mod1(&new_ops.translation);
                    queue.push_back(Operation::new(new_ops.rotation, new_translation_mod1));
                }
            }
        }

        operations
    }

    pub fn primitive_traverse(&self) -> Operations {
        let (_, primitive_operations) = self.traverse_and_primitive_traverse();
        primitive_operations
    }

    pub fn traverse_and_primitive_traverse(&self) -> (Operations, Operations) {
        let operations = self.traverse();
        let primitive_operations = Transformation::from_linear(self.centering.linear())
            .inverse_transform_operations(&operations);
        (operations, primitive_operations)
    }

    pub fn from_hall_number(hall_number: HallNumber) -> Option<Self> {
        if let Some(entry) = hall_symbol_entry(hall_number) {
            Self::new(entry.hall_symbol)
        } else {
            None
        }
    }

    pub fn primitive_generators(&self) -> Operations {
        Transformation::from_linear(self.centering.linear())
            .inverse_transform_operations(&self.generators)
    }
}

// Magnetic Hall symbol [J. Appl. Cryst. (2021). 54, 338-342]
//
// Extended Backus-Naur form (EBNF) for Magnetic Hall symbols
// ----------------------------------------------------------
// <Magnetic Hall symbol> := <L> <N-MSG>+ <V>?
// <N-MSG>                := <nfold> <A>? <T>? "'"?  # "'" represents time reversal
#[derive(Debug)]
pub struct MagneticHallSymbol {
    pub magnetic_hall_symbol: String,
    pub centering: Centering,
    pub centering_translations: Vec<Translation>,
    /// Generators of the magnetic space group except pure translations, which may include an anti translation
    pub generators: MagneticOperations,
}

impl MagneticHallSymbol {
    pub fn new(magnetic_hall_symbol: &str) -> Option<Self> {
        let tokens = tokenize(magnetic_hall_symbol);
        let (inversion_at_origin, lattice_symbol, ns, origin_shift) = parse(&tokens)?;

        // From translation subgroup
        let centering_translations = get_centering_translations(&lattice_symbol);

        // change basis by (I, v)
        // (R, tau) -> (R, tau + v - Rv)
        let mut generators = vec![];

        if inversion_at_origin {
            generators.push(MagneticOperation::new(
                -Rotation::identity(),
                2.0 * origin_shift,
                false,
            ));
        }

        for (rotation, translation, time_reversal) in ns {
            let translation_mod1 = (translation + origin_shift
                - rotation.map(|e| e as f64) * origin_shift)
                .map(|e| e.rem_euclid(1.0));
            generators.push(MagneticOperation::new(
                rotation,
                translation_mod1,
                time_reversal,
            ));
        }

        Some(Self {
            magnetic_hall_symbol: magnetic_hall_symbol.to_string(),
            centering: lattice_symbol,
            centering_translations,
            generators,
        })
    }

    /// Traverse all the magnetic symmetry operations up to translations by conventional cell.
    /// The order of operations is guaranteed to be fixed.
    pub fn traverse(&self) -> MagneticOperations {
        let mut queue = VecDeque::new();
        let mut hm_translations = HashMap::new();
        let mut operations = vec![];

        queue.push_back(MagneticOperation::identity());

        while !queue.is_empty() {
            let ops_lhs = queue.pop_front().unwrap();
            let entry = hm_translations.entry((ops_lhs.operation.rotation, ops_lhs.time_reversal));
            if let Entry::Occupied(_) = entry {
                continue;
            }
            entry.or_insert(ops_lhs.operation.translation);
            operations.push(ops_lhs.clone());

            for rhs in self.generators.iter() {
                let new_ops = ops_lhs.clone() * rhs.clone();
                let new_translation_mod1 = purify_translation_mod1(&new_ops.operation.translation);

                if !hm_translations
                    .contains_key(&(new_ops.operation.rotation, new_ops.time_reversal))
                {
                    queue.push_back(MagneticOperation::new(
                        new_ops.operation.rotation,
                        new_translation_mod1,
                        new_ops.time_reversal,
                    ));
                }
            }
        }

        operations
    }

    pub fn primitive_traverse(&self) -> MagneticOperations {
        let operations = self.traverse();
        Transformation::from_linear(self.centering.linear())
            .inverse_transform_magnetic_operations(&operations)
    }

    pub fn from_uni_number(uni_number: UNINumber) -> Option<Self> {
        if let Some(entry) = magnetic_hall_symbol_entry(uni_number) {
            Self::new(entry.magnetic_hall_symbol)
        } else {
            None
        }
    }

    pub fn primitive_generators(&self) -> MagneticOperations {
        Transformation::from_linear(self.centering.linear())
            .inverse_transform_magnetic_operations(&self.generators)
    }
}

/// Tokenize string by whitespaces
fn tokenize(hall_symbol: &str) -> Vec<&str> {
    let tokens = hall_symbol
        .split_whitespace()
        .filter(|s| !s.is_empty())
        .collect();
    tokens
}

#[allow(clippy::type_complexity)]
fn parse(
    tokens: &Vec<&str>,
) -> Option<(
    bool,
    Centering,
    Vec<(Rotation, Translation, TimeReversal)>,
    OriginShift,
)> {
    let (inversion_at_origin, lattice_symbol) = parse_lattice(tokens[0])?;

    let mut ns = vec![];
    let mut rotation_count = 0;
    let mut origin_shift = Vector3::<f64>::zeros();
    let mut prev_nfold = String::new();
    let mut prev_axis = String::new();

    for cursor in 1..tokens.len() {
        if tokens[cursor].starts_with('(') {
            // token = (vx, vy, vz)
            origin_shift = parse_origin_shift(tokens.iter().skip(cursor).cloned().collect())?;
            break;
        } else {
            // We need to remember the location of "N" symbols and preceding "N" symbol
            // because its default axis depends on them!
            let (rotation, translation, time_reversal, nfold, axis) =
                parse_operation(tokens[cursor], rotation_count, prev_nfold, prev_axis)?;
            ns.push((rotation, translation, time_reversal));
            prev_axis = axis;
            prev_nfold = nfold;
            rotation_count += 1;
        }
    }

    Some((inversion_at_origin, lattice_symbol, ns, origin_shift))
}

fn parse_lattice(token: &str) -> Option<(bool, Centering)> {
    let mut pos = 0;
    let inversion_at_origin = match token.chars().nth(pos).unwrap() {
        '-' => {
            pos += 1;
            true
        }
        _ => false,
    };
    let lattice_symbol = match token.chars().nth(pos).unwrap() {
        'P' => Centering::P,
        'A' => Centering::A,
        'B' => Centering::B,
        'C' => Centering::C,
        'I' => Centering::I,
        'R' => Centering::R,
        'F' => Centering::F,
        _ => return None,
    };
    Some((inversion_at_origin, lattice_symbol))
}

fn parse_origin_shift(tokens: Vec<&str>) -> Option<Vector3<f64>> {
    // Trim '(' and ')'
    let tokens = tokens
        .iter()
        .map(|s| {
            if let Some(stripped) = s.strip_prefix('(') {
                stripped
            } else if s.chars().nth(s.len() - 1).unwrap() == ')' {
                &s[..s.len() - 1]
            } else {
                s
            }
        })
        .filter(|s| !s.is_empty())
        .collect::<Vec<&str>>();
    if tokens.len() != 3 {
        return None;
    }
    let origin_shift = Vector3::<f64>::new(
        tokens[0].parse::<f64>().unwrap() / MAX_DENOMINATOR as f64,
        tokens[1].parse::<f64>().unwrap() / MAX_DENOMINATOR as f64,
        tokens[2].parse::<f64>().unwrap() / MAX_DENOMINATOR as f64,
    );
    Some(origin_shift)
}

fn parse_operation(
    token: &str,
    count: usize,
    prev_nfold: String,
    prev_axis: String,
) -> Option<(Rotation, Translation, TimeReversal, String, String)> {
    let mut pos = 0;

    let improper = match token.chars().nth(pos).unwrap() {
        '-' => {
            pos += 1;
            true
        }
        _ => false,
    };

    let nfold = token.chars().nth(pos).unwrap().to_string();
    pos += 1;

    let mut axis = String::new();
    if pos < token.len() {
        if token.chars().nth(pos).unwrap() == '^' {
            axis += "p";
            pos += 1;
        } else if token.chars().nth(pos).unwrap() == '=' {
            axis += "pp";
            pos += 1;
        }
    }

    if pos < token.len() {
        let c = token.chars().nth(pos).unwrap();
        if (c == 'x') || (c == 'y') || (c == 'z') || (c == '*') {
            axis.push(c);
            pos += 1;
        }
    }
    if ((axis == "p") || (axis == "pp"))
        && ((prev_axis == "x") || (prev_axis == "y") || (prev_axis == "z"))
    {
        // See Table A1.4.2.5
        axis += prev_axis.as_str();
    }

    if nfold == "1" {
        axis.push('z');
    }

    // Default axes. See A1.4.2.3.1
    if axis.is_empty() || (axis == "p") || (axis == "pp") {
        match count {
            0 => {
                // axis direction of c
                axis.push('z');
            }
            1 => {
                if (prev_nfold == "2") || (prev_nfold == "4") {
                    // axis direction of a
                    axis.push('x');
                } else if (prev_nfold == "3") || (prev_nfold == "6") {
                    // axis direction of a-b
                    axis += "pz";
                } else {
                    return None;
                }
            }
            2 => {
                if nfold == "3" {
                    // axis direction of a+b+c
                    axis.push('*');
                } else {
                    return None;
                }
            }
            _ => {
                return None;
            }
        }
    }

    let mut rotation = parse_rotation_matrix(format!("{}{}", nfold, axis).as_str())?;
    if improper {
        rotation *= -1;
    }

    // translation and time reversal
    let mut translation = Translation::zeros();
    let mut time_reversal = false;

    while pos < token.len() {
        let c = token.chars().nth(pos).unwrap();
        // translations are applied additively
        if "123456".contains(c) {
            // always along z-axis!
            translation = Translation::new(
                0.0,
                0.0,
                c.to_string().parse::<f64>().unwrap() / nfold.parse::<f64>().unwrap(),
            );
        } else if "abcnuvwd".contains(c) {
            translation += parse_translation_vector(c)?;
        } else if c == '\'' {
            time_reversal = true;
        } else {
            break;
        }
        pos += 1;
    }

    assert_eq!(pos, token.len());
    Some((
        rotation,
        translation,
        time_reversal,
        nfold,
        axis.to_string(),
    ))
}

fn parse_rotation_matrix(axis: &str) -> Option<Rotation> {
    match axis {
        "1x" | "1y" | "1z" => Some(matrix![
            1, 0, 0;
            0, 1, 0;
            0, 0, 1;
        ]),
        "2x" => Some(matrix![
            1, 0, 0;
            0, -1, 0;
            0, 0, -1;
        ]),
        "2y" => Some(matrix![
            -1, 0, 0;
            0, 1, 0;
            0, 0, -1;
        ]),
        "2z" => Some(matrix![
            -1, 0, 0;
            0, -1, 0;
            0, 0, 1;
        ]),
        "3x" => Some(matrix![
            1, 0, 0;
            0, 0, -1;
            0, 1, -1;
        ]),
        "3y" => Some(matrix![
            -1, 0, 1;
            0, 1, 0;
            -1, 0, 0;
        ]),
        "3z" => Some(matrix![
            0, -1, 0;
            1, -1, 0;
            0, 0, 1;
        ]),
        "4x" => Some(matrix![
            1, 0, 0;
            0, 0, -1;
            0, 1, 0;
        ]),
        "4y" => Some(matrix![
            0, 0, 1;
            0, 1, 0;
            -1, 0, 0;
        ]),
        "4z" => Some(matrix![
            0, -1, 0;
            1, 0, 0;
            0, 0, 1;
        ]),
        "6x" => Some(matrix![
            1, 0, 0;
            0, 1, -1;
            0, 1, 0;
        ]),
        "6y" => Some(matrix![
            0, 0, 1;
            0, 1, 0;
            -1, 0, 1;
        ]),
        "6z" => Some(matrix![
            1, -1, 0;
            1, 0, 0;
            0, 0, 1;
        ]),
        "2px" => Some(matrix![
            -1, 0, 0;
            0, 0, -1;
            0, -1, 0;
        ]),
        "2ppx" => Some(matrix![
            -1, 0, 0;
            0, 0, 1;
            0, 1, 0;
        ]),
        "2py" => Some(matrix![
            0, 0, -1;
            0, -1, 0;
            -1, 0, 0;
        ]),
        "2ppy" => Some(matrix![
            0, 0, 1;
            0, -1, 0;
            1, 0, 0;
        ]),
        "2pz" => Some(matrix![
            0, -1, 0;
            -1, 0, 0;
            0, 0, -1;
        ]),
        "2ppz" => Some(matrix![
            0, 1, 0;
            1, 0, 0;
            0, 0, -1;
        ]),
        "3*" => Some(matrix![
            0, 0, 1;
            1, 0, 0;
            0, 1, 0;
        ]),
        _ => None,
    }
}

fn parse_translation_vector(symbol: char) -> Option<Translation> {
    match symbol {
        'a' => Some(Translation::new(1.0 / 2.0, 0.0, 0.0)),
        'b' => Some(Translation::new(0.0, 1.0 / 2.0, 0.0)),
        'c' => Some(Translation::new(0.0, 0.0, 1.0 / 2.0)),
        'n' => Some(Translation::new(1.0 / 2.0, 1.0 / 2.0, 1.0 / 2.0)),
        'u' => Some(Translation::new(1.0 / 4.0, 0.0, 0.0)),
        'v' => Some(Translation::new(0.0, 1.0 / 4.0, 0.0)),
        'w' => Some(Translation::new(0.0, 0.0, 1.0 / 4.0)),
        'd' => Some(Translation::new(1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0)),
        _ => None,
    }
}

fn get_centering_translations(lattice_symbol: &Centering) -> Vec<Translation> {
    lattice_symbol
        .lattice_points()
        .iter()
        .filter(|&&translation| relative_ne!(translation, Translation::zeros(), epsilon = EPS))
        .cloned()
        .collect::<Vec<_>>()
}

fn purify_translation_mod1(translation: &Translation) -> Translation {
    translation.map(|e| {
        let mut eint = (e * (MAX_DENOMINATOR as f64)).round() as i32;
        eint = eint.rem_euclid(MAX_DENOMINATOR);
        (eint as f64) / (MAX_DENOMINATOR as f64)
    })
}

#[cfg(test)]
mod tests {
    use nalgebra::{matrix, vector};
    use rstest::rstest;
    use test_log::test as test_with_log;

    use super::*;

    #[rstest]
    #[case("P 2 2ab -1ab", Centering::P, 0, 3, 8)] // No. 51
    #[case("P 31 2 (0 0 4)", Centering::P, 0, 2, 6)] // No. 151
    #[case("P 65", Centering::P, 0, 1, 6)] // No. 170
    #[case("P 61 2 (0 0 5)", Centering::P, 0, 2, 12)] // No. 178
    #[case("-P 6c 2c", Centering::P, 0, 3, 24)] // No. 194
    #[case("F 4d 2 3", Centering::F, 3, 3, 24)] // No. 210
    fn test_hall_symbol_small(
        #[case] hall_symbol: &str,
        #[case] lattice_symbol: Centering,
        #[case] num_centering_translations: usize,
        #[case] num_generators: usize,
        #[case] num_operations: usize, // without centerings
    ) {
        let hs = HallSymbol::new(hall_symbol).unwrap();
        assert_eq!(hs.centering, lattice_symbol);
        assert_eq!(hs.centering_translations.len(), num_centering_translations);
        assert_eq!(hs.generators.len(), num_generators);
        let operations = hs.traverse();
        assert_eq!(operations.len(), num_operations);
    }

    #[test_with_log]
    fn test_hall_symbol_generators() {
        // No. 178
        let hs = HallSymbol::new("P 61 2 (0 0 5)").unwrap();
        let generators = hs.generators;
        assert_eq!(generators.len(), 2);
        assert_eq!(
            generators[0].rotation,
            matrix![
                1, -1, 0;
                1, 0, 0;
                0, 0, 1;
            ]
        );
        assert_relative_eq!(generators[0].translation, vector![0.0, 0.0, 1.0 / 6.0]);
        assert_eq!(
            generators[1].rotation,
            matrix![
                0, -1, 0;
                -1, 0, 0;
                0, 0, -1;
            ]
        );
        assert_relative_eq!(generators[1].translation, vector![0.0, 0.0, 5.0 / 6.0]);
    }

    #[rstest]
    #[case("C 2c -2 1c'", Centering::C, 1, 3, 8)] // 36.177 (type-IV)
    #[case("P 31 2 1c' (0 0 4)", Centering::P, 0, 3, 12)] // 151.32 (type-IV)
    #[case("P 6c 2c' -1'", Centering::P, 0, 3, 24)] // 194.265 (type-III)
    #[case("F 4d 2 3 1'", Centering::F, 3, 4, 48)] // 210.53 (type-II)
    fn test_magnetic_hall_symbol_small(
        #[case] magnetic_hall_symbol: &str,
        #[case] lattice_symbol: Centering,
        #[case] num_centering_translations: usize,
        #[case] num_generators: usize,
        #[case] num_operations: usize, // without centerings
    ) {
        let mhs = MagneticHallSymbol::new(magnetic_hall_symbol).unwrap();
        assert_eq!(mhs.centering, lattice_symbol);
        assert_eq!(mhs.centering_translations.len(), num_centering_translations);
        assert_eq!(mhs.generators.len(), num_generators);
        let magnetic_operations = mhs.traverse();
        assert_eq!(magnetic_operations.len(), num_operations);
    }
}
