use nalgebra::{Matrix3, Vector3};

struct WyckoffPosition {
    /// Wyckoff letter
    pub letter: char,
    /// Multiplicity in conventional cell
    pub multiplicity: usize,
    pub site_symmetry: &'static str,
    pub coordinates: &'static str,
}

struct WyckoffPositionSpace {
    pub linear: Matrix3<i32>,
    pub origin: Vector3<f64>,
}

impl WyckoffPositionSpace {
    /// Parse short-hand notation
    /// EBNF for Wyckoff coordinates is as follows (we always ignore space):
    ///     <shorthand>   ::= <term>, <term>, <term>
    ///     <term>        ::= "-"?<factor> ([+-]<factor>)* ([+-]<translation>)?
    ///     <factor>      ::= <integer>? <variable>
    ///     <variable>    ::= "x" || "y" || "z"
    ///     <translation> ::= <integer> ("/" <integer>)?
    ///     <integer>     ::= ("0" || "1" || "2" || "3" || "4" || "5" || "6" || "7" || "8" || "9")+

    pub fn new(coordinates: &str) -> Self {
        let coordinates = coordinates.replace(' ', "");
        let terms = coordinates.split(",").collect::<Vec<_>>();
        assert_eq!(terms.len(), 3);

        let mut linear = Matrix3::zeros();
        let mut origin = Vector3::zeros();
        let variables = ['x', 'y', 'z'];
        for (i, term) in terms.iter().enumerate() {
            let mut tokens_with_sign = vec![]; // like [(-1, "y"), (1, "2x"), (-1, "1/2")]
            let mut sign = 1;
            let mut token = vec![];
            for c in term.chars() {
                if c == '+' {
                    assert!(!token.is_empty());
                    tokens_with_sign.push((sign, token.clone().into_iter().collect::<String>()));
                    sign = 1;
                    token.clear();
                } else if c == '-' {
                    if !token.is_empty() {
                        tokens_with_sign
                            .push((sign, token.clone().into_iter().collect::<String>()));
                        token.clear();
                    }
                    sign = -1;
                } else {
                    token.push(c);
                }
            }
            if !token.is_empty() {
                tokens_with_sign.push((sign, token.clone().into_iter().collect::<String>()));
            }

            for (sign, token) in tokens_with_sign {
                if token.chars().last().unwrap().is_digit(10) {
                    // translation case
                    let nums = token.split('/').collect::<Vec<_>>();
                    if nums.len() == 1 {
                        origin[i] += (sign as f64) * token.parse::<f64>().unwrap();
                    } else {
                        let numerator = nums[0].parse::<f64>().unwrap();
                        let denominator = nums[1].parse::<f64>().unwrap();
                        origin[i] += (sign as f64) * numerator / denominator;
                    }
                } else {
                    // variable case
                    for j in 0..3 {
                        if token.chars().last().unwrap() != variables[j] {
                            continue;
                        }
                        let coeff = if token.chars().count() - 1 == 0 {
                            1
                        } else {
                            token[..token.len() - 1].parse::<i32>().unwrap()
                        };
                        linear[(i, j)] += sign * coeff;
                    }
                }
            }
        }

        Self { linear, origin }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::WyckoffPositionSpace;
    use nalgebra::{matrix, vector};

    #[rstest]
    #[case("-y, x, z+1/2", matrix![0, -1, 0; 1, 0, 0; 0, 0, 1], vector![0.0, 0.0, 0.5])]
    #[case("x,x-y+1/4,z+1/4", matrix![1, 0, 0; 1, -1, 0; 0, 0, 1], vector![0.0, 0.25, 0.25])]
    #[case("-x+2z,y,z", matrix![-1, 0, 2; 0, 1, 0; 0, 0, 1], vector![0.0, 0.0, 0.0])]
    #[case("1/4,1/4,1/4", matrix![0, 0, 0; 0, 0, 0; 0, 0, 0], vector![0.25, 0.25, 0.25])]
    fn test_wyckoff_position_space(
        #[case] coordinates: &str,
        #[case] linear: nalgebra::Matrix3<i32>,
        #[case] origin: nalgebra::Vector3<f64>,
    ) {
        let space = WyckoffPositionSpace::new(coordinates);
        assert_eq!(space.linear, linear);
        assert_relative_eq!(space.origin, origin);
    }
}
