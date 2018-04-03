use std::f32::{INFINITY, NEG_INFINITY};

use ::cgmath::Vector3;

use super::intersect::IntersectRay;

/// An axis-aligned bounding box in 3D
#[derive(Debug, Copy, Clone)]
pub struct Aabb {
    pub min: Vector3<f32>,
    pub max: Vector3<f32>
}

impl Aabb {
    /// Creates the smallest aabb that encloses all of the points returned
    /// by the given iterator.
    /// Returns an aabb with max at negative infinity and min at positive infinity if
    /// the given iterator was empty.
    pub fn from_points<P>(points: P) -> Aabb
        where P: IntoIterator<Item = Vector3<f32>>
    {
        points.into_iter()
            .fold(
                Aabb {
                    min: Vector3::new(INFINITY, INFINITY, INFINITY),
                    max: Vector3::new(NEG_INFINITY, NEG_INFINITY, NEG_INFINITY)
                },
                |Aabb { min, max }, p| {
                    let min_x = if p.x < min.x { p.x } else { min.x };
                    let min_y = if p.y < min.y { p.y } else { min.y };
                    let min_z = if p.z < min.z { p.z } else { min.z };

                    let max_x = if p.x > max.x { p.x } else { max.x };
                    let max_y = if p.y > max.y { p.y } else { max.y };
                    let max_z = if p.z > max.z { p.z } else { max.z };

                    Aabb {
                        min: Vector3::new(min_x, min_y, min_z),
                        max: Vector3::new(max_x, max_y, max_z)
                    }
                }
            )
    }

    /// Returns the smallest aabb that encloses all of the aabb in the given iterator.
    /// Returns an aabb with max at negative infinity and min at positive infinity if
    /// the given iterator was empty.
    pub fn union<A>(aabbs: A) -> Aabb
        where A: IntoIterator<Item = Aabb>
    {
        aabbs.into_iter()
            .fold(
                Aabb {
                    min: Vector3::new(INFINITY, INFINITY, INFINITY),
                    max: Vector3::new(NEG_INFINITY, NEG_INFINITY, NEG_INFINITY)
                },
                |Aabb { min: acc_min, max: acc_max }, Aabb { min: aabb_min, max: aabb_max }| {
                    let min_x = if aabb_min.x < acc_min.x { aabb_min.x } else { acc_min.x };
                    let min_y = if aabb_min.y < acc_min.y { aabb_min.y } else { acc_min.y };
                    let min_z = if aabb_min.z < acc_min.z { aabb_min.z } else { acc_min.z };

                    let max_x = if aabb_max.x > acc_max.x { aabb_max.x } else { acc_max.x };
                    let max_y = if aabb_max.y > acc_max.y { aabb_max.y } else { acc_max.y };
                    let max_z = if aabb_max.z > acc_max.z { aabb_max.z } else { acc_max.z };

                    Aabb {
                        min: Vector3::new(min_x, min_y, min_z),
                        max: Vector3::new(max_x, max_y, max_z)
                    }
                }
            )
    }

    fn is_point_outside(&self, point: Vector3<f32>) -> bool {
        point.x < self.min.x || point.x > self.max.x ||
            point.y < self.min.y || point.y > self.max.y ||
            point.z < self.min.z || point.z > self.max.z
    }

    pub fn is_point_inside(&self, point: Vector3<f32>) -> bool {
        !self.is_point_outside(point)
    }

    pub fn is_aabb_inside(&self, other: &Aabb) -> bool {
        self.is_point_inside(other.min) && self.is_point_inside(other.max)
    }

    pub fn intersects_ray(&self, ray_origin: Vector3<f32>, ray_direction: Vector3<f32>) -> bool {
        match self.line_intersection_min_max_parameters(ray_origin, ray_direction) {
            // If one is > 0, ray originates inside, if two are > 0 ray intersects from the outside
            Some((min_t, max_t)) => min_t > 0.0 || max_t > 0.0,
            None => false
        }
    }

    pub fn volume(&self) -> f32 {
        let dims = self.max - self.min;
        dims.x * dims.y * dims.z
    }

    /// Finds the smallest and the highest value for t in the line equation `line_origin = t * line_dir` where the given line
    /// intersects with the given axis-aligned bounding box. Since a line and not a ray is given, negative values for t are
    /// also considered an intersection.
    ///
    /// The line is specified with the given origin and direction.
    ///
    /// If the line does not intersect the aabb, None is returned.
    ///
    /// See: https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection
    fn line_intersection_min_max_parameters(&self, line_origin: Vector3<f32>, line_dir: Vector3<f32>) -> Option<(f32, f32)> {
        let line_dir_inv = 1.0 / line_dir;
        let line_sign = [
            if line_dir_inv.x < 0.0 { 1 } else { 0 },
            if line_dir_inv.y < 0.0 { 1 } else { 0 },
            if line_dir_inv.z < 0.0 { 1 } else { 0 }
        ];
        let bounds = [
            self.min,
            self.max
        ];

        let mut tmin = (bounds[line_sign[0]].x - line_origin.x) * line_dir_inv.x;
        let mut tmax = (bounds[1-line_sign[0]].x - line_origin.x) * line_dir_inv.x;

        let tymin = (bounds[line_sign[1]].y - line_origin.y) * line_dir_inv.y;
        let tymax = (bounds[1-line_sign[1]].y - line_origin.y) * line_dir_inv.y;

        if (tmin > tymax) || (tymin > tmax) {
            return None;
        }
        if tymin > tmin {
            tmin = tymin;
        }
        if tymax < tmax {
            tmax = tymax;
        }

        let tzmin = (bounds[line_sign[2]].z - line_origin.z) * line_dir_inv.z;
        let tzmax = (bounds[1-line_sign[2]].z - line_origin.z) * line_dir_inv.z;

        if (tmin > tzmax) || (tzmin > tmax) {
            return None;
        }
        if tzmin > tmin {
            tmin = tzmin;
        }
        if tzmax < tmax {
            tmax = tzmax;
        }

        Some((tmin, tmax))
    }
}

impl IntersectRay for Aabb {
    /// Finds the closest value for t in the ray equation ray(t) = ray_origin + t * ray_direction
    /// for the called object or None if no intersection
    fn ray_intersection_parameter(&self, ray_origin: Vector3<f32>, ray_direction: Vector3<f32>) -> Option<f32> {
        if let Some((tmin, tmax)) = self.line_intersection_min_max_parameters(ray_origin, ray_direction) {
            if tmin >= 0.0 {
                Some(tmin)
            } else if tmax >= 0.0 {
                // tmin negative, tmax positive, origin is inside AABB, closest t is zero
                Some(0.0)
            } else {
                // tmin and tmax are negative, intersection is in opposite direction and does not count
                None
            }
        } else {
            None
        }
    }
}

#[cfg(test)]
mod test {

    use super::*;
    use std::iter;

    #[test]
    fn test_aabb_from_points_empty() {
        let aabb = Aabb::from_points(iter::empty());

        assert_eq!(aabb.min.x, INFINITY, "Expected infinite AABB from empty points");
        assert_eq!(aabb.min.y, INFINITY, "Expected infinite AABB from empty points");
        assert_eq!(aabb.min.z, INFINITY, "Expected infinite AABB from empty points");

        assert_eq!(aabb.max.x, NEG_INFINITY, "Expected infinite AABB from empty points");
        assert_eq!(aabb.max.y, NEG_INFINITY, "Expected infinite AABB from empty points");
        assert_eq!(aabb.max.z, NEG_INFINITY, "Expected infinite AABB from empty points");
    }

    #[test]
    fn test_aabb_from_single_point() {
        let point = Vector3::new(1.0, 2.0, 3.0);
        let aabb = Aabb::from_points(iter::once(point));

        assert_eq!(aabb.min, point, "Built AABB from single point {:?} and expected min to be equal, but was {:?}", point, aabb.min);
        assert_eq!(aabb.max, point, "Built AABB from single point {:?} and expected max to be equal, but was {:?}", point, aabb.max);
    }

    #[test]
    fn test_aabb_from_points_triangle() {
        let aabb = Aabb::from_points(vec![
            Vector3::new(-0.5, -0.5, 1.0),
            Vector3::new(0.5, -0.5, 1.0),
            Vector3::new(0.0, 0.5, -1.0)
        ]);

        assert_eq!(aabb.min, Vector3::new(-0.5, -0.5, -1.0));
        assert_eq!(aabb.max, Vector3::new(0.5, 0.5, 1.0));
    }

    #[test]
    fn test_inside() {
        let aabb = Aabb::from_points(vec![
            Vector3::new(-0.5, -0.5, 1.0),
            Vector3::new(0.5, -0.5, 1.0),
            Vector3::new(0.0, 0.5, -1.0)
        ]);

        let point = Vector3::new(0.0, 0.0, 0.0);
        assert!(aabb.is_point_inside(point));

        let point = Vector3::new(10.0, 0.0, 0.0);
        assert!(!aabb.is_point_inside(point));

        let other_aabb = Aabb {
            min: Vector3::new(-0.1, -0.1, -0.1),
            max: Vector3::new(0.1, 0.1, 0.1),
        };
        assert!(aabb.is_aabb_inside(&other_aabb));

        let other_aabb = Aabb {
            min: Vector3::new(99.9, -0.1, -0.1),
            max: Vector3::new(100.1, 0.1, 0.1),
        };
        assert!(!aabb.is_aabb_inside(&other_aabb));
    }

    #[test]
    fn test_aabb_line_intersection() {
        // Given one aabb around the origin
        let aabb = Aabb::from_points(vec![
            Vector3::new(-0.5, -0.5, 1.0),
            Vector3::new(0.5, -0.5, 1.0),
            Vector3::new(0.0, 0.5, -1.0)
        ]);

        // Line above the origin, with direction facing up, hit
        let line_origin = Vector3::new(0.0, 10.0, 0.0);
        let line_direction = Vector3::new(0.0, -1.0, 0.0);
        let intersection_params = aabb.line_intersection_min_max_parameters(line_origin, line_direction);
        assert!(intersection_params.is_some(), "Line above the origin, with direction facing up, should hit an aabb below");
        let (min_t, max_t) = intersection_params.unwrap();
        assert!(min_t > 0.0 && max_t > 0.0 && max_t > min_t);

        // Line above the origin, with direction facing down should also hit an aabb below, this is not a ray
        let line_origin = Vector3::new(0.0, 10.0, 0.0);
        let line_direction = Vector3::new(0.0, 1.0, 0.0);
        let intersection_params = aabb.line_intersection_min_max_parameters(line_origin, line_direction);
        assert!(intersection_params.is_some(), "Line above the origin, with direction facing down should also hit an aabb below, this is not a ray");
        let (min_t, max_t) = intersection_params.unwrap();
        assert!(min_t < 0.0 && max_t < 0.0 && max_t > min_t);

        // Line above the origin, facing right, should miss
        let line_origin = Vector3::new(0.0, 10.0, 0.0);
        let line_direction = Vector3::new(1.0, 0.0, 0.0);
        let intersection_params = aabb.line_intersection_min_max_parameters(line_origin, line_direction);
        assert!(intersection_params.is_none(), "Line above the origin, facing right, should miss");

        // Line above the origin, facing towards Z, should miss
        let line_origin = Vector3::new(0.0, 10.0, 0.0);
        let line_direction = Vector3::new(0.0, 0.0, 1.0);
        let intersection_params = aabb.line_intersection_min_max_parameters(line_origin, line_direction);
        assert!(intersection_params.is_none(), "Line above the origin, facing right, should miss");

    }

    #[test]
    fn test_aabb_intersects_ray_inside() {
        // Given one aabb around the origin
        let aabb = Aabb::from_points(vec![
            Vector3::new(-0.5, -0.5, -1.0),
            Vector3::new(0.5, 0.5, 1.0)
        ]);

        let ray_origin = Vector3::new(0.0, 0.0, 0.0);
        let ray_direction = Vector3::new(1.0, 1.0, 0.0);
        assert!(aabb.is_point_inside(ray_origin));
        assert!(aabb.intersects_ray(ray_origin, ray_direction), "No matter the direction, a ray originating inside an aabb should always hit it");
        let ray_direction = Vector3::new(-1.0, 0.0, 1.0);
        assert!(aabb.is_point_inside(ray_origin));
        assert!(aabb.intersects_ray(ray_origin, ray_direction), "No matter the direction, a ray originating inside an aabb should always hit it");
        let ray_direction = Vector3::new(0.0, 1.0, 1.0);
        assert!(aabb.is_point_inside(ray_origin));
        assert!(aabb.intersects_ray(ray_origin, ray_direction), "No matter the direction, a ray originating inside an aabb should always hit it");
    }

    #[test]
    fn test_aabb_intersects_ray_hit_from_outside() {
        // Given one aabb around the origin
        let aabb = Aabb::from_points(vec![
            Vector3::new(-0.5, -0.5, -1.0),
            Vector3::new(0.5, 0.5, 1.0)
        ]);

        let ray_origin = Vector3::new(10.0, 0.0, 0.0);
        let ray_direction = Vector3::new(-1.0, 0.0, 0.0);
        assert!(aabb.intersects_ray(ray_origin, ray_direction), "Ray shot from the right to an abb on the left should hit");
    }

    #[test]
    fn test_aabb_intersects_miss() {
        // Given one aabb around the origin
        let aabb = Aabb::from_points(vec![
            Vector3::new(-0.5, -0.5, -1.0),
            Vector3::new(0.5, 0.5, 1.0)
        ]);

        let ray_origin = Vector3::new(10.0, 0.0, 0.0);
        let ray_direction = Vector3::new(1.0, 0.0, 0.0);
        assert!(!aabb.intersects_ray(ray_origin, ray_direction), "Ray shot from the right to the right should miss an aabb on the left");
    }
}
