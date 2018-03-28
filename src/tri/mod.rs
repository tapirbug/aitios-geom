pub mod iter;

mod interpolation;
mod intersect;
mod split;
mod tangent;
mod tri;
mod tuple;

pub use self::interpolation::{Interpolation, InterpolateVertex};
pub use self::intersect::*;
pub use self::tri::*;
pub use self::tangent::*;
pub use self::tuple::TupleTriangle;
pub use self::split::split_at_edge_midpoints;
