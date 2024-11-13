//! Get rules on simplices.

use crate::types::NumericalQuadratureDefinition;
use crate::types::QuadratureError;
use itertools::Itertools;
use paste::paste;

macro_rules! simplex_rule {
    ($cell:ident, $dim:tt) => {
        paste! {
use crate::simplex_rule_definitions::[<SIMPLEX_RULE_DEFINITIONS_ $cell:upper>];

/// Return a simplex rule on a [<$cell:lower>] for a given number of points.
///
/// If the rule does not exist `Err(())` is returned.
pub fn [<simplex_rule_ $cell:lower>](
    npoints: usize,
) -> Result<NumericalQuadratureDefinition, QuadratureError> {
    if let Some((order, points, weights)) = [<SIMPLEX_RULE_DEFINITIONS_ $cell:upper>]
        .get(&npoints)
    {
        Ok(NumericalQuadratureDefinition {
            dim: [<$dim>],
            order: *order,
            npoints,
            weights: weights.to_vec(),
            points: points.to_vec(),
        })
    } else {
        Err(QuadratureError::RuleNotFound)
    }
}


/// Return a vector with the numbers of points for which simplex rules are available on a [<$cell:lower>].
pub fn [<available_rules_ $cell:lower>]() -> Vec<usize> {
    [<SIMPLEX_RULE_DEFINITIONS_ $cell:upper>].iter().map(|(npoints, _)| *npoints).collect_vec()
}

        }
    };
}

simplex_rule!(interval, 1);
simplex_rule!(triangle, 2);
simplex_rule!(quadrilateral, 2);
simplex_rule!(tetrahedron, 3);
simplex_rule!(hexahedron, 3);
simplex_rule!(prism, 3);
simplex_rule!(pyramid, 3);

#[cfg(test)]
mod test {
    use super::*;
    use approx::*;

    fn volume_interval() -> f64 {
        1.0
    }
    fn volume_triangle() -> f64 {
        0.5
    }
    fn volume_quadrilateral() -> f64 {
        1.0
    }
    fn volume_tetrahedron() -> f64 {
        1.0 / 6.0
    }
    fn volume_hexahedron() -> f64 {
        1.0
    }
    fn volume_prism() -> f64 {
        0.5
    }
    fn volume_pyramid() -> f64 {
        1.0 / 3.0
    }

    macro_rules! test_cell {
        ($cell:ident) => {
            paste! {

                #[test]
                fn [<test_volume_ $cell:lower>]() {
                    let rules = [<available_rules_ $cell:lower>]();
                    for npoints in rules {
                        let rule = [<simplex_rule_ $cell:lower>](npoints).unwrap();
                        let volume_actual: f64 = rule.weights.iter().sum();
                        assert_relative_eq!(volume_actual, [<volume_ $cell:lower>](), max_relative=1E-14);
                        }

                }

            }
        };
    }

    test_cell!(interval);
    test_cell!(triangle);
    test_cell!(quadrilateral);
    test_cell!(tetrahedron);
    test_cell!(hexahedron);
    test_cell!(prism);
    test_cell!(pyramid);
}
