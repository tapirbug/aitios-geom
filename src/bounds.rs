//! Exposes the `Bounds` type for objects that have spatial dimensions.

use super::aabb::Aabb;
use linalg::Vec3;
use tri::Triangle;

/// Implemented by types that have spatial dimenions in three-dimensional
/// euclidean space.
///
/// The bounds can be accessed as an axis-aligned bounding box.
pub trait Bounds {
    fn bounds(&self) -> Aabb;
}

impl<T: Triangle> Bounds for T {
    fn bounds(&self) -> Aabb {
        let (a, b, c) = self.positions();
        Aabb {
            min: Vec3::new(
                a.x.min(b.x).min(c.x),
                a.y.min(b.y).min(c.y),
                a.z.min(b.z).min(c.z),
            ),
            max: Vec3::new(
                a.x.max(b.x).max(c.x),
                a.y.max(b.y).max(c.y),
                a.z.max(b.z).max(c.z),
            ),
        }
    }
}

impl Bounds for Aabb {
    fn bounds(&self) -> Aabb {
        *self
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use linalg::Vec3;
    use tri::FromVertices;
    use tri::TupleTriangle as Tri;

    #[test]
    fn test_bounds() {
        // cw (backside) on X/Y-Plane when looking in negative Z direction, normal should point in negative Z direction
        let tri = Tri::new(
            Vec3::new(-1.0, 1.0, 0.0),
            Vec3::new(1.0, 1.0, 0.0),
            Vec3::new(0.0, -1.0, 0.0),
        );

        assert_ulps_eq!(Vec3::new(-1.0, -1.0, 0.0), tri.bounds().min);
        assert_ulps_eq!(Vec3::new(1.0, 1.0, 0.0), tri.bounds().max);
    }
}
