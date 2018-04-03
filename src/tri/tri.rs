use super::iter::TriangleVertexIter;

use ::linalg::Vec3;
use ::vtx::Position;
use ::cgmath::prelude::*;
use std::f32::EPSILON;

/// Represents a triangle in three-dimensional space.
///
/// Its individual vertices can be obtained with the required
/// `vertices` method, which is the only method required for implementors.
/// As a minimum, returned vertices are expected to implement [`Position`](../vtx/trait.Position.html).
///
/// If a custom implementation is not desired, *aitios_geom* provides the [`TupleTriangle`](struct.TupleTriangle.html)
/// type for a triangle consisting of three vertices of the same type.
///
/// The trait provides some common triangle functionality, such as calculating the area.
/// For further functionality, the trait can be used as a bound for templated functions.
///
/// # Examples
///
/// For basic usage, have a look at calculating the area of a `TupleTriangle<Vec3>`:
/// ```
/// # #[macro_use]
/// # extern crate aitios_geom;
/// use aitios_geom::prelude::*;
/// use aitios_geom::{Vec3, TupleTriangle};
///
/// // A right triangle with short sides of length 5 and 8.
/// # fn main() {
/// let side_a = 5.0;
/// let side_b = 8.0;
/// let expected_area = 0.5 * side_a * side_b;
///
/// let triangle = TupleTriangle(
///     Vec3::new(0.0, 0.0, 0.0), // Vec3 implements `Position` by returning itself
///     Vec3::new(side_a, 0.0, 0.0),
///     Vec3::new(0.0, side_b, 0.0)
/// );
///
/// assert_ulps_eq!(expected_area, triangle.area());
/// # }
/// ```
///
/// Since `Triangle` is a trait, it plays nicely as a trait bound for templated functions
/// or types. This way, you can use it to build additional triangle functionality
/// not provided by the crate.
///
/// For instance, there is no `circumference` method provided in the triangle trait.
/// Should you need such functionality, you could easily build a function that works with any
/// type implementing `Triangle`:
///
/// ```
/// # #[macro_use]
/// # extern crate aitios_geom;
/// use aitios_geom::{
///     Vec3,
///     Triangle,
///     TupleTriangle
/// };
/// use aitios_geom::prelude::*;
///
/// /// Calculates the circumference of any triangle.
/// fn circumference<T>(triangle: &T) -> f32
///     where T : Triangle
/// {
///     let (v0, v1, v2) = triangle.positions();
///     let (a, b, c) = (
///         (v2 - v1).magnitude(),
///         (v2 - v0).magnitude(),
///         (v1 - v0).magnitude()
///     );
///     a + b + c
/// }
///
/// # fn main() {
/// // Make a triangle with
/// // a = √2, b = 1, c = 1
/// let triangle = TupleTriangle(
///     Vec3::new(0.0, 0.0, 0.0),
///     Vec3::new(1.0, 0.0, 0.0),
///     Vec3::new(0.0, 1.0, 0.0)
/// );
///
/// assert_ulps_eq!(
///     f32::sqrt(2.0) + 1.0 + 1.0,
///     circumference(&triangle)
/// );
/// # }
/// ```
pub trait Triangle {
    type Vertex : Position;

    fn vertices(&self) -> (&Self::Vertex, &Self::Vertex, &Self::Vertex);

    fn to_array(&self) -> [&Self::Vertex; 3] {
        let (a, b, c) = self.vertices();
        [a, b, c]
    }

    fn positions(&self) -> (Vec3, Vec3, Vec3) {
        let (a, b, c) = self.vertices();
        (a.position(), b.position(), c.position())
    }

    /// ||(V1 −V0)×(V2 −V0)||/2
    fn area(&self) -> f32 {
        let (v0, v1, v2) = self.positions();

        0.5 * ((v1 - v0).cross(v2 - v0)).magnitude()
    }

    fn centroid(&self) -> Vec3 {
        let one_over_three =  1.0 / 3.0;
        let (a, b, c) = self.positions();

        one_over_three * a +
        one_over_three * b +
        one_over_three * c
    }

    /// Gets the center of a sphere that runs through all of the three
    /// triangle vertices.
    fn circumcenter(&self) -> Vec3 {
        let (a, b, c) = self.positions();

        let ac = c - a;
        let ab = b - a;
        let ab_cross_ac = ab.cross(ac);

        // this is the vector from a vertex A to the circumsphere center
        let to_circumsphere_center = (ab_cross_ac.cross(ab) * ac.magnitude2() + ac.cross(ab_cross_ac) * ab.magnitude2()) /
            (2.0 * ab_cross_ac.magnitude2());

        a +  to_circumsphere_center
    }

    /// Checks if the triangle is completely inside the given sphere
    fn is_inside_sphere(&self, center: Vec3, radius: f32) -> bool {
        let (a, b, c) = self.vertices();
        let radius_sqr = radius * radius;
        [a, b, c].iter()
            .all(|v| center.distance2(v.position()) < radius_sqr)
    }

    fn iter<'a>(&self) -> TriangleVertexIter<Self> {
        TriangleVertexIter::new(self)
    }

    /// Returns the minimum bounding sphere center and squared radius
    /// of the triangle.
    ///
    /// See: http://realtimecollisiondetection.net/blog/?p=20
    fn minimum_bounding_sphere_sqr(&self) -> (Vec3, f32) {
        //void MinimumBoundingCircle(Circle &circle, Point a, Point b, Point c) {
        let (a, b, c) = self.positions();
        let dot_abab = (b - a).dot(b - a);
        let dot_abac = (b - a).dot(c - a);
        let dot_acac = (c - a).dot(c - a);
        let d = 2.0 * (dot_abab * dot_acac - dot_abac * dot_abac);
        let mut reference_point = a;

        let center = if d.abs() <= EPSILON {
            let min = Vec3::new(
                a.x.min(b.x).min(c.x),
                a.y.min(b.y).min(c.y),
                a.z.min(b.z).min(c.z)
            );
            let max = Vec3::new(
                a.x.max(b.x).max(c.x),
                a.y.max(b.y).max(c.y),
                a.z.max(b.z).max(c.z)
            );

            // a, b, and c lie on a line. Circle center is center of AABB of the
            // points, and radius is distance from circle center to AABB corner
            reference_point = min;
            0.5 * (min + max)
        } else {
            let s = (dot_abab * dot_acac - dot_acac * dot_abac) / d;
            let t = (dot_acac * dot_abab - dot_abab * dot_abac) / d;
            // s controls height over AC, t over AB, (1-s-t) over BC
            if s <= 0.0 {
                0.5 * (a + c)
            } else if t <= 0.0 {
                0.5 * (a + b)
            } else if (s + t) >= 1.0 {
                reference_point = b;
                0.5 * (b + c)
            } else {
                a + s*(b - a) + t*(c - a)
            }
        };

        let radius_sqr = center.distance2(reference_point);

        (center, radius_sqr)
    }
}

pub trait FromVertices : Triangle {
    fn new(v0: Self::Vertex, v1: Self::Vertex, v2: Self::Vertex) -> Self;
}

#[cfg(test)]
mod test {
    use super::*;
    use ::tri::TupleTriangle;

    #[test]
    fn test_centroid() {
        let tri = TupleTriangle::new(Vec3::new(-10.0, 0.0, 0.0), Vec3::new(10.0, 0.0, 0.0), Vec3::new(0.0, 10.0, 0.0));
        assert_ulps_eq!(Vec3::new(0.0, 10.0/3.0, 0.0),  tri.centroid());
    }

    #[test]
    fn test_circumcenter() {
        let tri = TupleTriangle::new(
            Vec3::new(3.0, 2.0, 0.0),
            Vec3::new(1.0, 4.0, 0.0),
            Vec3::new(5.0, 4.0, 0.0)
        );
        let circumcenter = tri.circumcenter();

        assert_eq!(Vec3::new(3.0, 4.0, 0.0), circumcenter);
    }

    #[test]
    fn test_is_inside_sphere() {
        let tri = TupleTriangle::new(Vec3::new(-10.0, 0.0, 0.0), Vec3::new(10.0, 0.0, 0.0), Vec3::new(0.0, 10.0, 0.0));

        assert!(tri.is_inside_sphere(Vec3::zero(), 10.000001));
        assert!(!tri.is_inside_sphere(Vec3::zero(), 0.0));
        assert!(!tri.is_inside_sphere(Vec3::zero(), 9.9999999));
    }

    #[test]
    fn test_area() {
        let tri = TupleTriangle::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(10.0, 0.0, 0.0), Vec3::new(10.0, 10.0, 0.0));
        let area = tri.area();
        assert_eq!(50.0, area);
    }

    #[test]
    fn test_area_colinear() {
        let tri = TupleTriangle::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(10.0, 0.0, 0.0), Vec3::new(100.0, 0.0, 0.0));
        let area = tri.area();
        assert_eq!(0.0, area);
    }

    #[test]
    fn test_iter() {
        let vertex0 = Vec3::new(-1.0, -1.0, 0.0);
        let vertex1 = Vec3::new(1.0, -1.0, 0.0);
        let vertex2 = Vec3::new(0.0, 1.0, 0.0);
        let tri = TupleTriangle::new(vertex0, vertex1, vertex2);

        let mut iter = tri.iter();
        assert_eq!(vertex0, *iter.next().unwrap());
        assert_eq!(vertex1, *iter.next().unwrap());
        assert_eq!(vertex2, *iter.next().unwrap());
    }

    #[test]
    fn test_minimum_bounding_sphere() {
        let origin_offset = Vec3::new(100.0, -1000.35, 3.124145);
        let vertex0 = origin_offset + Vec3::new(3.0, -4.0, 0.0);
        let vertex1 = origin_offset + Vec3::new(3.0, -4.0, 0.0);
        let vertex2 = origin_offset - Vec3::new(3.0, -4.0, 0.0);
        let tri = TupleTriangle::new(vertex0, vertex1, vertex2);
        let (center, radius_sqr) = tri.minimum_bounding_sphere_sqr();
        assert_ulps_eq!(radius_sqr, 25.0);
        assert_ulps_eq!(origin_offset, center);
    }
}
