//!
//! Core types and traits for geometry, including re-exports for vectors and matrices
//! from [`cgmath`](../cgmath/index.html) and functionality for triangles, AABB and
//! intersection.
//!
//! This crate is used in the weathering engine *aitios*, which uses it to procedurally
//! generate weathered materials in scenes. You might find another use for it, though.
//!
//! # Algebra
//! For basic linear algebra such as vectors and matrics, the crate re-exports
//! [`cgmath`](../cgmath/index.html) types and provides type aliases for commonly used types,
//! e.g. `Vec2` is an alias for `cgmath::Vector2<f32>`.
//!
//! ## Example
//! This is a minimal example using this crate to test if the diagonal of a unit square
//! is really √2:
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
//! * [`Position`](vtx/trait.Position.html),
//! * [`Normal`](vtx/trait.Normal.html),
//! * [`Texcoords`](vtx/trait.Texcoords.html).
//!
//! Additionally, the [`Vertex`](vtx/struct.Vertex.html) type provides a standard type for
//! owned vertices implementing all of the above traits.
//!
//! # Triangles
//! The [`Triangle`](tri/trait.Triangle.html) trait can be implemented by types that can expose three
//! vertices which at least provide the `Position` trait. Based on this functionality,
//! implementing types get access to a host of provided methods of the `Triangle`
//! trait.
//!
//! ## Example
//!
//! ```
//! use aitios_geom::{
//!     Vec3,
//!     InnerSpace,
//!     Triangle,
//!     TupleTriangle
//! };
//!
//! /// Calculate tha area of a triangle by evaluating:
//! /// ||(V1 − V0) × (V2 − V0)|| / 2.
//! fn area<T>(triangle: T) -> f32
//!     where T: Triangle
//! {
//!     let (v0, v1, v2) = triangle.positions();
//!     0.5 * ((v1 - v0).cross(v2 - v0)).magnitude()
//! }
//!
//! # fn main() {
//! let some_triangle = TupleTriangle(
//!     // Instead of Vec3, use anything that implements `Position`,
//!     // even your own types.
//!     Vec3::new(0.0, 0.0, 42.0),
//!     Vec3::new(1.0, 0.0, 42.0),
//!     Vec3::new(0.0, 1.0, 42.0),
//! );
//! assert_eq!(
//!     area(some_triangle),
//!     0.5,
//!     "The triangle A(0,0) B(1,0) C(0,1) has exactly half the area of a unit square"
//! )
//! # }
//! ```
//!

#[cfg_attr(test, macro_use)]
extern crate cgmath;

// Export modules as public so qualified imports are possible
// Example:
// use aitios_geom::tri::{Triangle, Iter};
pub mod aabb;
pub mod bounds;
pub mod intersect;
pub mod linalg;
pub mod tri;
pub mod vtx;

// Make qualified imports an opt-in by bringing the most important
// types into top-level position so you can do
// Example:
// use aitios_geom::Triangle;
pub use aabb::Aabb;
pub use bounds::Bounds;
pub use intersect::IntersectRay;
pub use linalg::*;
pub use tri::{
    FromVertices, InterpolateVertex, Interpolation, TangentSpace, Triangle, TupleTriangle,
};
pub use vtx::{Normal, Position, Texcoords, Vertex};

pub mod prelude {
    //! Exposes the cgmath prelude, as well as `Bounds`, `IntersectRay`
    //! and some triangle traits.
    pub use bounds::Bounds;
    pub use cgmath::prelude::*;
    pub use intersect::IntersectRay;
    pub use tri::{FromVertices, InterpolateVertex, Interpolation, Triangle};
}
