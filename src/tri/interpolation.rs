use tri::{Triangle, FromVertices};
use vtx::Position;
use ::linalg::Vec3;
use std::ops::{Mul, Add};
use ::cgmath::prelude::*;

pub trait Interpolation : Triangle {
    /// Compute barycentric coordinates [u, v, w] for
    /// the closest point to p on the triangle.
    fn barycentric_at(&self, p: Vec3) -> [f32; 3] {
        let (a, b, c) = self.positions();

        let v0 = b - a;
        let v1 = c - a;
        let v2 = p - a;

        let d00 = v0.dot(v0);
        let d01 = v0.dot(v1);
        let d11 = v1.dot(v1);
        let d20 = v2.dot(v0);
        let d21 = v2.dot(v1);
        let denom = d00 * d11 - d01 * d01;

        let v = (d11 * d20 - d01 * d21) / denom;
        let w = (d00 * d21 - d01 * d20) / denom;
        let u = 1.0 - v - w;

        [u, v, w]
    }

    /// Uses the given function to extract a single value from each vertex and
    /// then interpolates a new value for the given position based on the vertex
    /// values.
    fn interpolate_at<F, T, M>(&self, position: Vec3, vertex_to_val_fn: F) -> M
        where F: Fn(&Self::Vertex) -> T,
            T: Mul<f32, Output = M>,
            M : Add<M, Output = M>
    {
        self.interpolate_bary(self.barycentric_at(position), vertex_to_val_fn)
    }

    /// Uses the given function to extract a single value from each vertex and
    /// then interpolates a new value for the given barycentric coordinates
    /// based on the vertex values.
    fn interpolate_bary<F, T, M>(&self, bary: [f32; 3], vertex_to_val_fn: F) -> M
        where F: Fn(&Self::Vertex) -> T,
            T: Mul<f32, Output = M>,
            M : Add<M, Output = M>
    {
        let (a, b, c) = self.vertices();

        vertex_to_val_fn(a) * bary[0] +
        vertex_to_val_fn(b) * bary[1] +
        vertex_to_val_fn(c) * bary[2]
    }
}

// REVIEW could move this into main Interpolation type but with
// fn interpolate… -> Self::Vertex
// where Self : FromVertices,
//       Self::Vertex : Position + Clone + Mul<f32, Output = V> + Add<V, Output = V>
// would save the need for a InterpolateVertex type
// on the other hand InterpolateVertex is useful for trait bounds and should at least
// remain as a tag trait with no methods
pub trait InterpolateVertex : Interpolation {
    fn interpolate_vertex_bary(&self, bary: [f32; 3]) -> Self::Vertex;

    fn interpolate_vertex_at(&self, position: Vec3) -> Self::Vertex {
        self.interpolate_vertex_bary(self.barycentric_at(position))
    }
}

impl<T : Triangle> Interpolation for T { }

impl<T, V> InterpolateVertex for T
    where T : FromVertices<Vertex = V>,
        V : Position + Clone + Mul<f32, Output = V> + Add<V, Output = V>
    {

    fn interpolate_vertex_bary(&self, bary: [f32; 3]) -> Self::Vertex {
        self.interpolate_bary(bary, |v| v.clone())
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ::tri::TupleTriangle;
    use ::vtx::Vertex;
    use ::linalg::{Vec2, Vec3};

    #[test]
    fn test_barycentric_at() {
        let tri = TupleTriangle::new(
            Vec3::new(3.0, 2.0, 0.0),
            Vec3::new(1.0, 4.0, 0.0),
            Vec3::new(5.0, 4.0, 0.0)
        );

        let (a, b, c) = tri.positions();
        let barys_at_a = tri.barycentric_at(a);
        let barys_at_b = tri.barycentric_at(b);
        let barys_at_c = tri.barycentric_at(c);
        let barys_at_centroid = tri.barycentric_at(tri.centroid());

        assert_eq!([1.0, 0.0, 0.0], barys_at_a);
        assert_eq!([0.0, 1.0, 0.0], barys_at_b);
        assert_eq!([0.0, 0.0, 1.0], barys_at_c);

        assert_ulps_eq!((1.0 / 3.0), barys_at_centroid[0]);
        assert_ulps_eq!((1.0 / 3.0), barys_at_centroid[1]);
        assert_ulps_eq!((1.0 / 3.0), barys_at_centroid[2]);
    }

    #[test]
    fn interpolate_position() {
        let vertex0 = Vec3::new(-1.0, -1.0, 0.0);
        let vertex1 = Vec3::new(1.0, -1.0, 0.0);
        let vertex2 = Vec3::new(0.0, 1.0, 0.0);
        let tri = TupleTriangle::new(vertex0, vertex1, vertex2);

        let point_on_there = Vec3::new(0.0, 0.5, 0.0);

        assert_eq!(
            point_on_there,
            tri.interpolate_at(point_on_there, |v| *v),
            "Interpolating the position value should yield the same point"
        );
    }

    #[test]
    fn interpolate_vertex() {
        let vertex0 = Vec3::new(-1.0, -1.0, 0.0);
        let vertex1 = Vec3::new(1.0, -1.0, 0.0);
        let vertex2 = Vec3::new(0.0, 1.0, 0.0);
        let tri = TupleTriangle::new(vertex0, vertex1, vertex2);

        let point_on_there = Vec3::new(0.0, 0.5, 0.0);

        // For a vertex type of Vec3 interpolating a vertex at a position value should yield the same point
        assert_ulps_eq!(
            point_on_there,
            tri.interpolate_vertex_at(point_on_there)
        );
    }

    #[test]
    fn interpolate_owned_vertex() {
        let vertex0 = Vertex {
            position: Vec3::new(-1.0, -1.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            texcoords: Vec2::new(0.0, 0.0),
        };
        let vertex1 = Vertex {
            position: Vec3::new(1.0, -1.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            texcoords: Vec2::new(1.0, 0.0),
        };
        let vertex2 = Vertex {
            position: Vec3::new(0.0, 1.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            texcoords: Vec2::new(0.5, 1.0),
        };
        let tri = TupleTriangle::new(vertex0, vertex1, vertex2);

        let on_there = Vec3::new(0.0, 0.0, 0.0);
        //assert_ulps_eq!(center, Vec3::new(0.0, -0.3333333, 0.0));

        let middle_vertex = tri.interpolate_vertex_at(on_there);

        // Expected vertex position interpolated at centroid to be equal to centroid
        assert_ulps_eq!(
            middle_vertex.position,
            on_there
        );

        // Expected texcoords interpolated at centroid to be 0.5/0.5 in the example
        assert_ulps_eq!(
            middle_vertex.texcoords,
            Vec2::new(0.5, 0.5)
        );

        // Expected normal to be unchanged since they are uniform on the triangle
        assert_ulps_eq!(
            middle_vertex.normal,
            Vec3::new(0.0, 0.0, 1.0)
        );
    }
}
