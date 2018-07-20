//! Contains the trait [`Triangle`](trait.Triangle.html), containing the core
//! functionality around triangles.
//!
//! Some more advanced functionality is available in other traits exposed by this
//! module. They need to be brought into scope first to be used.
//!
//! # Examples
//! Use [`TangentSpace`](trait.TangentSpace.html) to obtain the geometric normal
//! of a triangle:
//!
//! ```
//! use aitios_geom::{
//!     TupleTriangle, Vec3,
//!     TangentSpace, FromVertices
//! };
//!
//! let triangle = TupleTriangle(
//!     Vec3::new(0.0, 0.0, -42.0),
//!     Vec3::new(1.0, 0.0, -42.0),
//!     Vec3::new(0.0, 1.0, -42.0)
//! );
//!
//! assert_eq!(
//!     // `normal` is a function exposed by `TangentSpace`.
//!     triangle.normal(),
//!     Vec3::new(0.0, 0.0, 1.0),
//!     "The geometric normal of a 3D triangle points normal to counter-clockwise rotation order"
//! );
//!
//! assert_eq!(
//!     // `flip` is a function exposed by `FromVertices` and inverts the vertex ordering.
//!     triangle.flip()
//!         // Which in turn flips the geometric normal.
//!         .normal(),
//!     Vec3::new(0.0, 0.0, -1.0),
//!     "The geometric normal of a 3D triangle points normal to counter-clockwise rotation order"
//! );
//! ```

mod from_verts;
mod interpolation;
mod intersect;
mod iter;
mod split;
mod tangent;
mod tri;
mod tuple;

pub use self::from_verts::FromVertices;
pub use self::interpolation::{InterpolateVertex, Interpolation};
pub use self::intersect::*;
pub use self::iter::*;
pub use self::split::{split_at_edge_midpoints, Split};
pub use self::tangent::*;
pub use self::tri::*;
pub use self::tuple::TupleTriangle;
