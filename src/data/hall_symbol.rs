use itertools::Itertools;
use nalgebra::{matrix, Matrix3, Vector3};

use crate::base::lattice::Lattice;
use crate::base::operation::{Operations, Rotation, Translation};

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
    pub inversion_at_origin: bool,
    pub lattice_symbol: LatticeSymbol,
    pub generators: Operations,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum LatticeSymbol {
    P, // Primitive
    A, // A-face centered
    B, // B-face centered
    C, // C-face centered
    I, // Body centered
    R, // Rhombohedral (obverse setting)
    H, // Hexagonal
    F, // Face centered
}

impl HallSymbol {
    pub fn new(hall_symbol: &str) -> Option<Self> {
        let tokens = Self::tokenize(&hall_symbol);

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

        let mut rotations = vec![];
        let mut translations = vec![];
        // change basis by (I, v)
        // (R, tau) -> (R, tau + v - Rv)

        // From translation subgroup
        for translation in Self::lattice_points(lattice_symbol) {
            if relative_ne!(translation, Translation::zeros()) {
                rotations.push(Rotation::identity());
                translations.push(translation);
            }
        }

        if inversion_at_origin {
            rotations.push(-Rotation::identity());
            translations.push(2.0 * origin_shift);
        }

        for (rotation, translation) in ns {
            rotations.push(rotation);
            translations
                .push(translation + origin_shift - rotation.map(|e| e as f64) * translation);
        }

        Some(Self {
            hall_symbol: hall_symbol.to_string(),
            inversion_at_origin,
            lattice_symbol,
            generators: Operations::new(
                Lattice::new(Matrix3::<f64>::identity()),
                rotations,
                translations,
            ),
        })
    }

    fn tokenize(hall_symbol: &str) -> Vec<&str> {
        let tokens = hall_symbol
            .split_whitespace()
            .filter(|s| !s.is_empty())
            .collect();
        tokens
    }

    fn parse_lattice(token: &str) -> Option<(bool, LatticeSymbol)> {
        let mut pos = 0;
        let inversion_at_origin = match token.chars().nth(pos).unwrap() {
            '-' => {
                pos += 1;
                true
            }
            _ => false,
        };
        let lattice_symbol = match token.chars().nth(pos).unwrap() {
            'P' => LatticeSymbol::P,
            'A' => LatticeSymbol::A,
            'B' => LatticeSymbol::B,
            'C' => LatticeSymbol::C,
            'I' => LatticeSymbol::I,
            'R' => LatticeSymbol::R,
            'H' => LatticeSymbol::H,
            'F' => LatticeSymbol::F,
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
            tokens[0].parse::<i32>().unwrap() as f64 / 12.0,
            tokens[1].parse::<i32>().unwrap() as f64 / 12.0,
            tokens[2].parse::<i32>().unwrap() as f64 / 12.0,
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
        if pos <= token.len() - 1 {
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

        if pos <= token.len() - 1 {
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
        if (axis == "") || (axis == "p") || (axis == "pp") {
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

        while pos <= token.len() - 1 {
            let c = token.chars().nth(pos).unwrap();
            // translations are applied additively
            if "12346".contains(c) {
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

    fn lattice_points(lattice_symbol: LatticeSymbol) -> Vec<Vector3<f64>> {
        match lattice_symbol {
            LatticeSymbol::P => {
                vec![Translation::new(0.0, 0.0, 0.0)]
            }
            LatticeSymbol::A => {
                vec![
                    Translation::new(0.0, 0.0, 0.0),
                    Translation::new(0.0, 0.5, 0.5),
                ]
            }
            LatticeSymbol::B => {
                vec![
                    Translation::new(0.0, 0.0, 0.0),
                    Translation::new(0.5, 0.0, 0.5),
                ]
            }
            LatticeSymbol::C => {
                vec![
                    Translation::new(0.0, 0.0, 0.0),
                    Translation::new(0.5, 0.5, 0.0),
                ]
            }
            LatticeSymbol::I => {
                vec![
                    Translation::new(0.0, 0.0, 0.0),
                    Translation::new(0.5, 0.5, 0.5),
                ]
            }
            LatticeSymbol::R => {
                // obverse setting
                vec![
                    Translation::new(0.0, 0.0, 0.0),
                    Translation::new(2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0),
                    Translation::new(1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0),
                ]
            }
            LatticeSymbol::H => {
                vec![
                    Translation::new(0.0, 0.0, 0.0),
                    Translation::new(2.0 / 3.0, 1.0 / 3.0, 0.0),
                    Translation::new(1.0 / 3.0, 2.0 / 3.0, 0.0),
                ]
            }
            LatticeSymbol::F => {
                vec![
                    Translation::new(0.0, 0.0, 0.0),
                    Translation::new(0.0, 0.5, 0.5),
                    Translation::new(0.5, 0.0, 0.5),
                    Translation::new(0.5, 0.5, 0.0),
                ]
            }
        }
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
    use rstest::rstest;

    use super::{HallSymbol, LatticeSymbol};

    #[rstest]
    #[case("P 2 2ab -1ab", false, LatticeSymbol::P, 3)] // No. 51
    #[case("P 31 2 (0 0 4)", false, LatticeSymbol::P, 2)] // No. 151
    #[case("-P 6c 2c", true, LatticeSymbol::P, 3)] // No. 194
    #[case("F 4d 2 3", false, LatticeSymbol::F, 6)] // No. 210
    fn test_hall_symbol(
        #[case] hall_symbol: &str,
        #[case] inversion_at_origin: bool,
        #[case] lattice_symbol: LatticeSymbol,
        #[case] num_generators: usize,
    ) {
        let hs = HallSymbol::new(hall_symbol).unwrap();
        assert_eq!(hs.inversion_at_origin, inversion_at_origin);
        assert_eq!(hs.lattice_symbol, lattice_symbol);
        assert_eq!(hs.generators.num_operations(), num_generators);
        dbg!(&hs);
    }
}
