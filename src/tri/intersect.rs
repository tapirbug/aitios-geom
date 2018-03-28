use ::intersect::IntersectRay;
use super::Triangle;
use ::linalg::Vec3;
use ::cgmath::prelude::*;

/// Implements the moller trombore ray-triangle-intersection algorithm
impl<T : Triangle> IntersectRay for T
{
    fn ray_intersection_parameter(&self, ray_origin: Vec3, ray_direction: Vec3) -> Option<f32> {
        let (vertex0, vertex1, vertex2) = self.positions();

        let epsilon = 0.0000001;

        let edge1 = vertex1 - vertex0;
        let edge2 = vertex2 - vertex0;

        let h = ray_direction.cross(edge2);
        let a = edge1.dot(h);

        if a > -epsilon && a < epsilon {
            return None;
        }

        let f = 1.0 / a;
        let s = ray_origin - vertex0;
        let u = f * (s.dot(h));

        if u < 0.0 || u > 1.0 {
            return None;
        }

        let q = s.cross(edge1);
        let v = f * ray_direction.dot(q);

        if v < 0.0 || (u + v) > 1.0 {
            return None;
        }

        let t = f * edge2.dot(q);

        if t < epsilon {
            return None;
        }

        Some(t)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ::tri::TupleTriangle as Tri;
    use ::tri::FromVertices;

    #[test]
    fn intersect_ray_with_tri() {
        let ray_origin = Vec3::new(0.0, 0.0, 0.0);
        let ray_direction = Vec3::new(0.0, 0.0, 1.0);

        let vertex0 = Vec3::new(-1.0, -1.0, 100.0);
        let vertex1 = Vec3::new(1.0, -1.0, 100.0);
        let vertex2 = Vec3::new(0.0, 1.0, 200.0);

        let tri = Tri::new(vertex0, vertex1, vertex2);

        assert!(tri.ray_intersection_parameter(ray_origin, ray_direction).unwrap() > 0.0);
        assert_ulps_eq!(tri.ray_intersection_point(ray_origin, ray_direction).unwrap(), Vec3::new(0.0, 0.0, 150.0));
    }

    #[test]
    fn intersect_ray_with_tri_and_miss() {
        let ray_origin = Vec3::new(0.0, 0.0, 0.0);
        let ray_direction = Vec3::new(0.0, 0.0, -1.0);

        let vertex0 = Vec3::new(-1.0, -1.0, 100.0);
        let vertex1 = Vec3::new(1.0, -1.0, 100.0);
        let vertex2 = Vec3::new(0.0, 1.0, 200.0);

        let tri = Tri::new(vertex0, vertex1, vertex2);

        assert_eq!(tri.ray_intersection_parameter(ray_origin, ray_direction), None);
        assert_eq!(tri.ray_intersection_point(ray_origin, ray_direction), None);
    }
}
