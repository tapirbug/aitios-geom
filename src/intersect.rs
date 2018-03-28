
use ::cgmath::Vector3;

pub trait IntersectRay {
    /// Finds the closest value for t in the ray equation ray(t) = ray_origin + t * ray_direction
    /// for the called object or None if no intersection
    fn ray_intersection_parameter(&self, ray_origin: Vector3<f32>, ray_direction: Vector3<f32>) -> Option<f32>;

    fn ray_intersection_point(&self, ray_origin: Vector3<f32>, ray_direction: Vector3<f32>) -> Option<Vector3<f32>> {
        self.ray_intersection_parameter(ray_origin, ray_direction)
            .map(move |t| ray_origin + t * ray_direction)
    }
}
