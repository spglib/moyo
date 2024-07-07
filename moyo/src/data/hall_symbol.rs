use std::collections::hash_map::Entry;
use std::collections::{HashMap, VecDeque};

use nalgebra::{matrix, Matrix3, Vector3};
use strum_macros::EnumIter;

use super::hall_symbol_database::{hall_symbol_entry, HallNumber};
use crate::base::{Linear, Operations, Rotation, Transformation, Translation, EPS};

const MAX_DENOMINATOR: i32 = 12;

/// Hall symbol. See A1.4.2.3 in ITB (2010).
///
/// Extended Backus-Naur form (EBNF) for Hall symbols
/// ----------------------------------------------------------
/// <Hall symbol>    := <L> <N>+ <V>?
/// <L>              := "-"? <lattice symbol>
/// <lattice symbol> := [PABCIRHF]  # Table A1.4.2.2
/// <N>              := <nfold> <A>? <T>?
/// <nfold>          := "-"? ("1" | "2" | "3" | "4" | "6")
/// <A>              := [xyz] | "'" | '"' | "=" | "*"  # Table A.1.4.2.4, Table A.1.4.2.5, Table A.1.4.2.6
/// <T>              := [abcnuvwd] | [1-6]  # Table A.1.4.2.3
/// <V>              := "(" [0-11] [0-11] [0-11] ")"
#[derive(Debug)]
pub struct HallSymbol {
    pub hall_symbol: String,
    pub centering: Centering,
    pub centering_translations: Vec<Translation>,
    /// Generators of the space group except pure translations
    pub generators: Operations,
}

#[derive(Debug, Copy, Clone, PartialEq, EnumIter)]
pub enum Centering {
    P, // Primitive
    A, // A-face centered
    B, // B-face centered
    C, // C-face centered
    I, // Body centered
    R, // Rhombohedral (obverse setting)
    F, // Face centered
}

impl Centering {
    pub fn order(&self) -> usize {
        match self {
            Centering::P => 1,
            Centering::A => 2,
            Centering::B => 2,
            Centering::C => 2,
            Centering::I => 2,
            Centering::R => 3,
            Centering::F => 4,
        }
    }

    /// Inverse matrices of https://github.com/spglib/spglib/blob/39a95560dd831c2d16f162126921ac1e519efa31/src/spacegroup.c#L373-L384
    /// Transformation matrix from primitive to conventional cell.
    pub fn linear(&self) -> Linear {
        match self {
            Centering::P => Linear::identity(),
            Centering::A => Linear::new(
                1, 0, 0, //
                0, 1, 1, //
                0, -1, 1, //
            ),
            Centering::B => Linear::new(
                1, 0, -1, //
                0, 1, 0, //
                1, 0, 1, //
            ),
            Centering::C => Linear::new(
                1, -1, 0, //
                1, 1, 0, //
                0, 0, 1, //
            ),
            Centering::R => Linear::new(
                1, 0, 1, //
                -1, 1, 1, //
                0, -1, 1, //
            ),
            Centering::I => Linear::new(
                0, 1, 1, //
                1, 0, 1, //
                1, 1, 0, //
            ),
            Centering::F => Linear::new(
                -1, 1, 1, //
                1, -1, 1, //
                1, 1, -1, //
            ),
        }
    }

    /// Transformation matrix from conventional to primitive cell.
    pub fn inverse(&self) -> Matrix3<f64> {
        self.linear().map(|e| e as f64).try_inverse().unwrap()
    }

    pub fn lattice_points(&self) -> Vec<Vector3<f64>> {
        match self {
            Centering::P => {
                vec![Translation::new(0.0, 0.0, 0.0)]
            }
            Centering::A => {
                vec![
                    Translation::new(0.0, 0.0, 0.0),
                    Translation::new(0.0, 0.5, 0.5),
                ]
            }
            Centering::B => {
                vec![
                    Translation::new(0.0, 0.0, 0.0),
                    Translation::new(0.5, 0.0, 0.5),
                ]
            }
            Centering::C => {
                vec![
                    Translation::new(0.0, 0.0, 0.0),
                    Translation::new(0.5, 0.5, 0.0),
                ]
            }
            Centering::I => {
                vec![
                    Translation::new(0.0, 0.0, 0.0),
                    Translation::new(0.5, 0.5, 0.5),
                ]
            }
            Centering::R => {
                // obverse setting
                vec![
                    Translation::new(0.0, 0.0, 0.0),
                    Translation::new(2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0),
                    Translation::new(1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0),
                ]
            }
            // Centering::H => {
            //     vec![
            //         Translation::new(0.0, 0.0, 0.0),
            //         Translation::new(2.0 / 3.0, 1.0 / 3.0, 0.0),
            //         Translation::new(1.0 / 3.0, 2.0 / 3.0, 0.0),
            //     ]
            // }
            Centering::F => {
                vec![
                    Translation::new(0.0, 0.0, 0.0),
                    Translation::new(0.0, 0.5, 0.5),
                    Translation::new(0.5, 0.0, 0.5),
                    Translation::new(0.5, 0.5, 0.0),
                ]
            }
        }
    }
}

impl HallSymbol {
    pub fn new(hall_symbol: &str) -> Option<Self> {
        let tokens = Self::tokenize(hall_symbol);

        let (inversion_at_origin, lattice_symbol) = Self::parse_lattice(tokens[0])?;

        let mut ns = vec![];
        let mut rotation_count = 0;
        let mut origin_shift = Vector3::<f64>::zeros();
        let mut prev_nfold = String::new();
        let mut prev_axis = String::new();

        for cursor in 1..tokens.len() {
            if tokens[cursor].chars().nth(0).unwrap() == '(' {
                // token = (vx, vy, vz)
                origin_shift =
                    Self::parse_origin_shift(tokens.iter().skip(cursor).cloned().collect())?;
                break;
            } else {
                // We need to remember the location of "N" symbols and preceding "N" symbol
                // because its default axis depends on them!
                let (rotation, translation, nfold, axis) =
                    Self::parse_operation(tokens[cursor], rotation_count, prev_nfold, prev_axis)?;
                ns.push((rotation, translation));
                prev_axis = axis;
                prev_nfold = nfold;
                rotation_count += 1;
            }
        }

        // From translation subgroup
        let centering_translations = lattice_symbol
            .lattice_points()
            .iter()
            .filter(|&&translation| relative_ne!(translation, Translation::zeros(), epsilon = EPS))
            .cloned()
            .collect::<Vec<_>>();

        // change basis by (I, v)
        // (R, tau) -> (R, tau + v - Rv)
        let mut rotations = vec![];
        let mut translations = vec![];

        if inversion_at_origin {
            rotations.push(-Rotation::identity());
            translations.push(2.0 * origin_shift);
        }

        for (rotation, translation) in ns {
            rotations.push(rotation);
            let translation_mod1 = (translation + origin_shift
                - rotation.map(|e| e as f64) * origin_shift)
                .map(|e| e.rem_euclid(1.0));
            translations.push(translation_mod1);
        }

        Some(Self {
            hall_symbol: hall_symbol.to_string(),
            centering: lattice_symbol,
            centering_translations,
            generators: Operations::new(rotations, translations),
        })
    }

    /// Traverse all the symmetry operations up to translations by conventional cell.
    /// The order of operations is guaranteed to be fixed.
    pub fn traverse(&self) -> Operations {
        let mut queue = VecDeque::new();
        let mut hm_translations = HashMap::new();
        let mut rotations = vec![];
        let mut translations = vec![];

        queue.push_back((Rotation::identity(), Translation::zeros()));

        while !queue.is_empty() {
            let (rotation_lhs, translation_lhs) = queue.pop_front().unwrap();
            let entry = hm_translations.entry(rotation_lhs);
            if let Entry::Occupied(_) = entry {
                continue;
            }
            entry.or_insert(translation_lhs);
            rotations.push(rotation_lhs);
            translations.push(translation_lhs);

            for (rotation_rhs, translation_rhs) in self
                .generators
                .rotations
                .iter()
                .zip(self.generators.translations.iter())
            {
                let new_rotation = rotation_lhs * rotation_rhs;
                let new_translation =
                    rotation_lhs.map(|e| e as f64) * translation_rhs + translation_lhs;
                let new_translation_mod1 = new_translation.map(|e| {
                    let mut eint = (e * (MAX_DENOMINATOR as f64)).round() as i32;
                    eint = eint.rem_euclid(MAX_DENOMINATOR);
                    (eint as f64) / (MAX_DENOMINATOR as f64)
                });

                if !hm_translations.contains_key(&new_rotation) {
                    queue.push_back((new_rotation, new_translation_mod1));
                }
            }
        }

        Operations::new(rotations, translations)
    }

    pub fn from_hall_number(hall_number: HallNumber) -> Self {
        let entry = hall_symbol_entry(hall_number);
        Self::new(entry.hall_symbol).unwrap()
    }

    pub fn primitive_generators(&self) -> Operations {
        Transformation::from_linear(self.centering.linear())
            .inverse_transform_operations(&self.generators)
    }

    fn tokenize(hall_symbol: &str) -> Vec<&str> {
        let tokens = hall_symbol
            .split_whitespace()
            .filter(|s| !s.is_empty())
            .collect();
        tokens
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
                if s.chars().nth(0).unwrap() == '(' {
                    &s[1..]
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
    ) -> Option<(Rotation, Translation, String, String)> {
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

        let prime_axis_symbol = '\'';
        let mut axis = String::new();
        if pos < token.len() {
            if token.chars().nth(pos).unwrap() == prime_axis_symbol {
                axis += "p";
                pos += 1;
            } else if (token.chars().nth(pos).unwrap() == '"')
                || (token.chars().nth(pos).unwrap() == '=')
            {
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

        let mut rotation = Self::rotation_matrix(format!("{}{}", nfold, axis).as_str())?;
        if improper {
            rotation *= -1;
        }

        // translation
        let mut translation = Translation::zeros();

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
                translation += Self::translation_vector(c)?;
            } else {
                break;
            }
            pos += 1;
        }

        assert_eq!(pos, token.len());
        Some((rotation, translation, nfold, axis.to_string()))
    }

    fn rotation_matrix(axis: &str) -> Option<Rotation> {
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

    fn translation_vector(symbol: char) -> Option<Translation> {
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
}

#[cfg(test)]
mod tests {
    use nalgebra::{matrix, vector};
    use rstest::rstest;
    use strum::IntoEnumIterator;

    use super::{Centering, HallSymbol};
    use crate::base::Transformation;
    use crate::data::hall_symbol_database::iter_hall_symbol_entry;

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
        assert_eq!(hs.generators.num_operations(), num_generators);
        let operations = hs.traverse();
        assert_eq!(operations.num_operations(), num_operations);
    }

    #[test]
    fn test_hall_symbol_generators() {
        // No. 178
        let hs = HallSymbol::new("P 61 2 (0 0 5)").unwrap();
        let generators = hs.generators;
        assert_eq!(generators.num_operations(), 2);
        assert_eq!(
            generators.rotations[0],
            matrix![
                1, -1, 0;
                1, 0, 0;
                0, 0, 1;
            ]
        );
        assert_relative_eq!(generators.translations[0], vector![0.0, 0.0, 1.0 / 6.0]);
        assert_eq!(
            generators.rotations[1],
            matrix![
                0, -1, 0;
                -1, 0, 0;
                0, 0, -1;
            ]
        );
        assert_relative_eq!(generators.translations[1], vector![0.0, 0.0, 5.0 / 6.0]);
    }

    #[test]
    fn test_hall_symbol_whole() {
        for entry in iter_hall_symbol_entry() {
            let hs = HallSymbol::new(entry.hall_symbol).unwrap();
            assert_eq!(48 % hs.traverse().num_operations(), 0);
        }
    }

    #[test]
    fn test_conventional_transformation_matrix() {
        for centering in Centering::iter() {
            assert_eq!(
                Transformation::from_linear(centering.linear()).size,
                centering.order()
            );
        }
    }
}
