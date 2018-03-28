//!
//! Contains core types and traits for geometry.
//!
//! # Algebra
//! For basic linear algebra such as vectors and matrics, the crate re-exports *cgmath*
//! types and provides type aliases for commonly used types, e.g. `Vec2` is an alias
//! for `cgmath::Vector2<f32>`.
//!
//! ```
//! use aitios_geom::prelude::*;
//! use aitios_geom::Vec2;
//! use std::f32::consts::SQRT_2;
//!
//! assert_eq!(
//!     Vec2::new(1.0, 1.0).magnitude(),
//!     SQRT_2,
//!     "The diagonal of a unit squere is sqrt(2)"
//! );
//! ```
//!
//! # Vertices
//! There are some traits that are implemented by types representing 3D vertices or
//! point samples:
//! * `Position`,
//! * `Normal`,
//! * `Texcoords`.
//!
//! Additionally, the `Vertex` type provides a standard type for owned vertices
//! implementing all of the above traits.
//!
//! # Triangles
//! The [`Triangle`](tri/trait.Triangle.html) trait can be implemented by types that can expose three
//! vertices which at least provide the `Position` trait. Based on this functionality,
//! implementing types get access to a host of provided methods of the `Triangle`
//! trait.
//!
//!

#[cfg_attr(test, macro_use)]
extern crate cgmath;

pub mod tri;
pub mod linalg;
pub mod vtx;
pub mod aabb;
pub mod bounds;
pub mod intersect;

pub use linalg::*;
pub use vtx::*;
pub use tri::*;
pub use aabb::*;
pub use intersect::*;
pub use bounds::*;

pub mod prelude {
    pub use bounds::Bounds;
    pub use cgmath::prelude::*;
    pub use intersect::IntersectRay;
    pub use tri::{Triangle, FromVertices, Interpolation, InterpolateVertex};
}
