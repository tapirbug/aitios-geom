//! Defines a trait for ray intersections.
use cgmath::Vector3;

/// Implemented by types that can look for the closest intersection of a given
/// ray with the types geometry.
pub trait IntersectRay {
    /// Finds Some(t), where t is the closest value to zero in the ray equation `ray(t) = ray_origin + t * ray_direction`
    /// that intersects the called entity. Returns `None` if no intersection.
    fn ray_intersection_parameter(
        &self,
        ray_origin: Vector3<f32>,
        ray_direction: Vector3<f32>,
    ) -> Option<f32>;

    /// Finds the closest point on the given ray that intersects the called entity.
    fn ray_intersection_point(
        &self,
        ray_origin: Vector3<f32>,
        ray_direction: Vector3<f32>,
    ) -> Option<Vector3<f32>> {
        self.ray_intersection_parameter(ray_origin, ray_direction)
            .map(move |t| ray_origin + t * ray_direction)
    }
}
