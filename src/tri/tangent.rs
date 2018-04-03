use super::Triangle;
use ::linalg::{Mat3, Vec3};
use ::cgmath::prelude::*;

/// Implementors of this trait either have defined or can calculate:
///
/// * tangent,
/// * binormal,
/// * normal.
///
/// Implements some tangent space operations.
pub trait TangentSpace {
    /// Calculates a tangent space based on vertex positions and returns it as three
    /// vectors that form an orthonormal basis. The first vector will be a tangent,
    /// the second a binormal and the third the face normal.
    ///
    /// Note that there are infinite possible tangent spaces. The resulting tangent
    /// is parallel to the edge C, that is from vertex 0 to vertex 1. The
    /// tangent space is not guaranteed to be aligned with texture space, which is
    /// normally a common way to align it.
    fn tangent_space(&self) -> (Vec3, Vec3, Vec3);

    fn tangent(&self) -> Vec3 {
        let (tangent, _, _) = self.tangent_space();
        tangent
    }

    fn binormal(&self) -> Vec3 {
        let (_, binormal, _) = self.tangent_space();
        binormal
    }

    /// Calculates a face normal for the triangle based on the vertex positions
    /// and the cross product.
    ///
    /// Panics for empty triangles, since the resulting cross product is always zero.
    fn normal(&self) -> Vec3 {
        let (_, _, normal) = self.tangent_space();
        normal
    }

    fn world_to_tangent_matrix(&self) -> Mat3 {
        let (tangent, binormal, normal) = self.tangent_space();
        Mat3::from_cols(tangent, binormal, normal).transpose()
    }

    fn tangent_to_world_matrix(&self) -> Mat3 {
        let (tangent, binormal, normal) = self.tangent_space();
        Mat3::from_cols(tangent, binormal, normal)
    }

    /// Transforms the given direction vector into tangent space, setting the height component to zero and then
    /// transforming back into world space. The resulting direction should be parallel to the tangential plane
    /// in world space. If the given direction happens to be parallel to the normal, a zero vector is returned.
    fn project_onto_tangential_plane(&self, incoming_direction_world: Vec3) -> Vec3 {
        let world_to_tangent = self.world_to_tangent_matrix();
        let tangent_to_world = world_to_tangent.invert()
            .expect("Expected tangent space matrix to be invertible");

        let tangent_space_direction = world_to_tangent * incoming_direction_world;
        let tangent_space_direction_flat = tangent_space_direction
            .truncate() // drop Z
            .extend(0.0); // and set to zero

        let scaled_projected = tangent_to_world * tangent_space_direction_flat;

        if scaled_projected.is_zero() {
            scaled_projected
        } else {
            scaled_projected.normalize()
        }
    }
}

impl<T : Triangle> TangentSpace for T {
    /// Calculates a tangent space based on vertex positions and returns it as three
    /// vectors that form an orthonormal basis. The first vector will be a tangent,
    /// the second a binormal and the third the face normal.
    ///
    /// Note that there are infinite possible tangent spaces. The resulting tangent
    /// is parallel to the edge C, that is from vertex 0 to vertex 1. The
    /// tangent space is not guaranteed to be aligned with texture space, which is
    /// normally a common way to align it.
    fn tangent_space(&self) -> (Vec3, Vec3, Vec3) {
        let (a, b, c) = self.positions();

        let a_to_b = b - a;
        let a_to_c = c - a;

        let scaled_normal = a_to_b.cross(a_to_c);

        // If all three points are colinear, result is always zero vector: v.cross(v) = 0, v.cross(Vector3::zero()) = 0
        assert!(!scaled_normal.is_zero(), "Face normal is undefined for triangle with zero area: [{:?}, {:?}, {:?}]", a, b, c);

        let normal = scaled_normal.normalize();
        let tangent = a_to_b.normalize();
        let binormal = (tangent.cross(normal)).normalize();

        (tangent, binormal, normal)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ::tri::TupleTriangle as Tri;
    use ::tri::FromVertices;

    #[test]
    fn test_calculate_face_normal_from_positions() {
        // cw (backside) on X/Y-Plane when looking in negative Z direction, normal should point in negative Z direction
        let tri = Tri::new(
            Vec3::new(-1.0, 1.0, 0.0),
            Vec3::new(1.0, 1.0, 0.0),
            Vec3::new(0.0, -1.0, 0.0)
        );

        let (tangent, binormal, normal) = tri.tangent_space();
        assert_eq!(Vec3::new(0.0, 0.0, -1.0), normal);
        assert_eq!(Vec3::new(1.0, 0.0, 0.0), tangent);
        assert_eq!(Vec3::new(0.0, 1.0, 0.0), binormal);

        // ccw on X/Y-Plane when looking in negative Z direction, normal should point in positive z direction
        let tri = Tri::new(
            Vec3::new(-1.0, 1.0, 0.0),
            Vec3::new(0.0, -1.0, 0.0),
            Vec3::new(1.0, 1.0, 0.0)
        );


        assert_eq!(Vec3::new(0.0, 0.0, 1.0), tri.normal());
        assert_eq!((Vec3::new(0.0, -1.0, 0.0) - Vec3::new(-1.0, 1.0, 0.0)).normalize(), tri.tangent());
        assert_eq!(tri.binormal().dot(tri.normal()), 0.0);
    }

    #[test]
    #[should_panic]
    fn test_zero_area_triangle_normal_calculation_panics() {
        let tri = Tri::new(
            Vec3::new(-1.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 1.0)
        );

        tri.normal();
    }

    #[test]
    fn test_project_direction_on_tangential_plane_on_floor() {
        // triangle flat on the floor with normal facing toward y
        // projecting should just drop the Z value in this case
        let floor_tri = Tri::new(
            Vec3::new(-1.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 1.0),
            Vec3::new(0.0, 0.0, -1.0)
        );

        let up_right_positive_z = Vec3::new(1.0, 1.0, 1.0).normalize();

        assert_eq!(
            Vec3::new(1.0, 0.0, 1.0).normalize(),
            floor_tri.project_onto_tangential_plane(up_right_positive_z),
            "Projecting a vector onto a triangle flat on the floor should yield the same vector with the z value dropped"
        );
    }

    #[test]
    fn test_project_direction_on_slope_tri() {
        // Triangle with 45° upward slope in X direction
        let slope_tri = Tri::new(
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, -1.0)
        );

        let down = Vec3::new(0.0, -1.0, 0.0);
        let projected = slope_tri.project_onto_tangential_plane(down);

        assert_eq!(
            Vec3::new(-1.0, -1.0, 0.0).normalize(),
            projected
        );
    }

    #[test]
    fn test_project_direction_parallel_to_normal() {
        // triangle flat on the floor with normal facing toward y
        // projecting should just drop the Z value in this case
        let floor_tri = Tri::new(
            Vec3::new(-1.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 1.0),
            Vec3::new(0.0, 0.0, -1.0)
        );

        let up = Vec3::new(0.0, 1.0, 0.0).normalize();

        assert_eq!(
            Vec3::new(0.0, 0.0, 0.0),
            floor_tri.project_onto_tangential_plane(up),
            "Projecting a vector onto a triangle flat on the floor should yield the same vector with the z value dropped"
        );
    }
}