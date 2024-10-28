//! Print out all available quadrature rules

use bempp_quadrature::simplex_rules::available_rules;
use strum::IntoEnumIterator;

pub fn main() {
    for cell_type in ndelement::types::ReferenceCellType::iter() {
        let mut rules = available_rules(cell_type);
        rules.as_mut().unwrap().sort();
        println!("Available rules for {:?}: {:?}", cell_type, rules);
    }
}
