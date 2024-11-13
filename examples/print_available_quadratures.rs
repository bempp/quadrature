//! Print out all available quadrature rules

use bempp_quadrature::simplex_rules::{
    available_rules_hexahedron, available_rules_interval, available_rules_prism,
    available_rules_pyramid, available_rules_quadrilateral, available_rules_tetrahedron,
    available_rules_triangle,
};

pub fn main() {
    let mut rules = available_rules_interval();
    rules.sort();
    println!("Available rules for interval: {:?}", rules);

    let mut rules = available_rules_triangle();
    rules.sort();
    println!("Available rules for triangle: {:?}", rules);

    let mut rules = available_rules_quadrilateral();
    rules.sort();
    println!("Available rules for quadrilateral: {:?}", rules);

    let mut rules = available_rules_tetrahedron();
    rules.sort();
    println!("Available rules for tetrahedron: {:?}", rules);

    let mut rules = available_rules_hexahedron();
    rules.sort();
    println!("Available rules for hexahedron: {:?}", rules);

    let mut rules = available_rules_prism();
    rules.sort();
    println!("Available rules for prism: {:?}", rules);

    let mut rules = available_rules_pyramid();
    rules.sort();
    println!("Available rules for pyramid: {:?}", rules);
}
