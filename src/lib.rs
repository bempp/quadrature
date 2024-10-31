//! Quadrature
#![cfg_attr(feature = "strict", deny(warnings), deny(unused_crate_dependencies))]
#![warn(missing_docs)]

pub mod duffy;
pub mod simplex_rule_definitions;
pub mod simplex_rules;
pub mod types;

#[cfg(test)]
mod test {
    use strum as _; // Hack to show that strum is used, as cargo test is not picking up dependencies in examples
}

