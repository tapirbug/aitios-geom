use std::f32::{INFINITY, NEG_INFINITY};
use std::iter::FromIterator;
use std::borrow::Borrow;
use linalg::Vec3;
use intersect::IntersectRay;

/// An axis-aligned bounding box in 3D
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Aabb {
    pub min: Vec3,
    pub max: Vec3
}

impl Aabb {
    /// Creates the smallest AABB encompassing all of the points
    /// in the given iterator over points or point references.
    ///
    /// If the iterator only has a single point, the AABB will have
    /// min and max at that point. An empty iterator yields `Aabb::empty()`.
    ///
    /// ```
    /// use aitios_geom::{Aabb, Vec3};
    ///
    /// let p1 = Vec3::new(-0.5, -0.5, 1.0);
    /// let p2 = Vec3::new(0.5, -0.5, 1.0);
    /// let p3 = Vec3::new(0.0, 0.5, -1.0);
    ///
    /// // Consumes the points in the vector
    /// let aabb = Aabb::from_points(vec![p1, p2, p3]);
    ///
    /// // All points are contained
    /// assert!(aabb.is_point_inside(p1) && aabb.is_point_inside(p2) && aabb.is_point_inside(p3));
    ///
    /// // Any kind of reference is okay too
    /// let aabb = Aabb::from_points(vec![p1, p2, p3].iter());
    /// let aabb = Aabb::from_points(vec![&p1, &p2, &p3]);
    /// let aabb = Aabb::from_points(vec![Box::new(p1), Box::new(p2), Box::new(p3)]);
    /// ```
    pub fn from_points<I, P>(iter: I) -> Self
        where I : IntoIterator<Item=P>,
            P : Borrow<Vec3>
    {
        iter.into_iter()
            .fold(
                Aabb::empty(),
                |Aabb { min, max }, p| {
                    let &Vec3 { x, y, z} = p.borrow();
                    let min = Vec3::new(
                        min.x.min(x),
                        min.y.min(y),
                        min.z.min(z)
                    );
                    let max = Vec3::new(
                        max.x.max(x),
                        max.y.max(y),
                        max.z.max(z)
                    );
                    Aabb { min, max }
                }
            )
    }

    /// Constructs an AABB that covers no possible point. Its min point
    /// is set to positive infinity, and max to negative infinity.
    /// Hence, it passes every outside check and fails every inside
    /// check.
    ///
    /// There are many possible representations of the empty AABB, but
    /// this one is often the nicest for fold operations.
    pub fn empty() -> Self {
        Aabb {
            min: Vec3::new(INFINITY, INFINITY, INFINITY),
            max: Vec3::new(NEG_INFINITY, NEG_INFINITY, NEG_INFINITY)
        }
    }

    /// Checks whether the Aabb is empty, i.e. has a zero volume.
    pub fn is_empty(&self) -> bool {
        self.min.x >= self.max.x ||
        self.min.y >= self.max.y ||
        self.min.z >= self.max.z
    }

    /// Creates the largest possible AABB with the min point set at negative
    /// infinity and max at positive infinity. It will pass every inside check
    /// and fail every outside check.
    pub fn infinite() -> Self {
        Aabb {
            min: Vec3::new(NEG_INFINITY, NEG_INFINITY, NEG_INFINITY),
            max: Vec3::new(INFINITY, INFINITY, INFINITY)
        }
    }

    pub fn is_infinite(&self) -> bool {
        self == &Self::infinite()
    }

    /// Returns the smallest aabb that encloses self and the given aabb.
    /// Returns an aabb with max at negative infinity and min at positive infinity if
    /// the given iterator was empty.
    pub fn union(&self, other: &Self) -> Self
    {
        Aabb {
            min: Vec3::new(
                self.min.x.min(other.min.x),
                self.min.y.min(other.min.y),
                self.min.z.min(other.min.z)
            ),
            max: Vec3::new(
                self.max.x.max(other.max.x),
                self.max.y.max(other.max.y),
                self.max.z.max(other.max.z)
            )
        }
    }

    fn is_point_outside(&self, point: Vec3) -> bool {
        point.x < self.min.x || point.x > self.max.x ||
            point.y < self.min.y || point.y > self.max.y ||
            point.z < self.min.z || point.z > self.max.z
    }

    pub fn is_point_inside(&self, point: Vec3) -> bool {
        !self.is_point_outside(point)
    }

    pub fn is_aabb_inside(&self, other: &Aabb) -> bool {
        self.is_point_inside(other.min) && self.is_point_inside(other.max)
    }

    pub fn intersects_ray(&self, ray_origin: Vec3, ray_direction: Vec3) -> bool {
        match self.line_intersection_min_max_parameters(ray_origin, ray_direction) {
            // If one is > 0, ray originates inside, if two are > 0 ray intersects from the outside
            Some((min_t, max_t)) => min_t > 0.0 || max_t > 0.0,
            None => false
        }
    }

    /// Returns the volume of the AABB.
    ///
    /// If the min point is larger than the max point in any dimension,
    /// a volume of zero will be returned.
    pub fn volume(&self) -> f32 {
        if self.is_empty() {
            0.0
        } else {
            let dims = self.max - self.min;
            dims.x * dims.y * dims.z
        }
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
    fn line_intersection_min_max_parameters(&self, line_origin: Vec3, line_dir: Vec3) -> Option<(f32, f32)> {
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

/// When created from an iterator over owned or borrowed Aabb, return the smallest possible Aabb
/// encompassing all the given Aabb.
impl<T> FromIterator<T> for Aabb
    where T : Borrow<Aabb>
{
    fn from_iter<I>(iter: I) -> Self
        where I : IntoIterator<Item=T>
    {
        iter.into_iter()
            .fold(
                Aabb::empty(),
                |acc, a| acc.union(a.borrow())
            )
    }
}

// FIXME why cannot implement multiple times?
/*
/// When created from an iterator over owned or borrowed points, return the smallest possible Aabb
/// encompassing all the given points.
///
/// ```
/// use aitios_geom::{Vec3, Aabb};
///
/// let points = vec![
///     Vec3::new(0.0, 0.0, 0.0),
///     Vec3::new(0.5, 0.5, 0.5),
///     Vec3::new(1.0, 1.0, 1.0)
/// ];
/// let aabb : Aabb = points.iter().collect();
///
/// assert!(
///     aabb.is_point_inside(points[0]) &&
///     aabb.is_point_inside(points[1]) &&
///     aabb.is_point_inside(points[2]))
/// );
/// ```
impl<T> FromIterator<T> for Aabb
    where T : Borrow<Vec3>
{
    fn from_iter<I>(iter: I) -> Self
        where I : IntoIterator<Item=T>
    {
        iter.into_iter()
            .fold(
                Aabb::empty(),
                |Aabb { min, max }, p| {
                    let &Vec3 { x, y, z} = p.borrow();
                    let min = Vec3::new(
                        min.x.min(x),
                        min.y.min(y),
                        min.z.min(z)
                    );
                    let max = Vec3::new(
                        max.x.max(x),
                        max.y.max(y),
                        max.z.max(z)
                    );
                    Aabb { min, max }
                }
            )
    }
}
*/

impl IntersectRay for Aabb {
    /// Finds the closest value for t in the ray equation ray(t) = ray_origin + t * ray_direction
    /// for the called object or None if no intersection
    fn ray_intersection_parameter(&self, ray_origin: Vec3, ray_direction: Vec3) -> Option<f32> {
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
    use std::f32::{MIN as LOWEST, MAX as HIGHEST};

    #[test]
    fn empty() {
        let empty = Aabb::empty();

        let normal_point = Vec3::new(0.0, 0.0, 1.0);
        let infinite = Vec3::new(INFINITY, INFINITY, INFINITY);
        let neg_infinite = Vec3::new(INFINITY, INFINITY, INFINITY);
        let faraway_point = Vec3::new(HIGHEST, HIGHEST, HIGHEST);
        let faraway_point2 = Vec3::new(HIGHEST, LOWEST, HIGHEST);

        assert!(empty.is_point_outside(normal_point));
        assert!(empty.is_point_outside(infinite));
        assert!(empty.is_point_outside(neg_infinite));
        assert!(empty.is_point_outside(faraway_point));
        assert!(empty.is_point_outside(faraway_point2));
        assert!(empty.is_empty());
        assert_eq!(empty.volume(), 0.0);

        let empty_on_y_axis = Aabb { min: Vec3::new(1.0, 1.0, 1.0), max: Vec3::new(10.0, 1.0, 100.0) };
        assert!(empty_on_y_axis.is_empty());

        let empty_same_point = Aabb { min: Vec3::new(1.0, 1.0, 1.0), max: Vec3::new(1.0, 1.0, 1.0) };
        assert!(empty_same_point.is_empty());

        let normal_aabb = Aabb { min: Vec3::new(0.0, 0.0, 0.0), max: Vec3::new(1.0, 1.0, 1.0) };
        assert_eq!(empty.union(&normal_aabb), normal_aabb);
    }

    #[test]
    fn infinities() {
        let infinite = Aabb::infinite();

        let normal_point = Vec3::new(0.0, 0.0, 1.0);
        let infinite_point = Vec3::new(INFINITY, INFINITY, INFINITY);
        let neg_infinite_point = Vec3::new(INFINITY, INFINITY, INFINITY);
        let faraway_point = Vec3::new(HIGHEST, HIGHEST, HIGHEST);
        let faraway_point2 = Vec3::new(HIGHEST, LOWEST, HIGHEST);

        assert!(infinite.is_point_inside(normal_point));
        assert!(infinite.is_point_inside(infinite_point));
        assert!(infinite.is_point_inside(neg_infinite_point));
        assert!(infinite.is_point_inside(faraway_point));
        assert!(infinite.is_point_inside(faraway_point2));
        assert!(infinite.is_infinite());
        assert_eq!(infinite.volume(), INFINITY);

        let normal_aabb = Aabb { min: Vec3::new(0.0, 0.0, 0.0), max: Vec3::new(1.0, 1.0, 1.0) };
        assert_eq!(infinite.union(&normal_aabb), infinite);
    }

    #[test]
    fn test_union_from_iter() {
        let aabbs = vec![
            Aabb {
                min: Vec3::new(-10.0, -10.0, -10.0),
                max: Vec3::new(-9.0, -9.0, -9.0),
            },
            Aabb {
                min: Vec3::new(1.0, 1.0, 1.0),
                max: Vec3::new(100.0, 100.0, 100.0),
            }
        ];

        let expected_combined = Aabb {
            min: Vec3::new(-10.0, -10.0, -10.0),
            max: Vec3::new(100.0, 100.0, 100.0),
        };

        // first, create from iterator over references
        let combined_from_refs : Aabb = aabbs.iter()
            .collect();

        assert_eq!(expected_combined, combined_from_refs);

        let combined_from_owned : Aabb = aabbs.into_iter()
            .map(|aabb : Aabb| aabb) // this wouldn't compile unless we are really dealing with owned values
            .collect();

        assert_eq!(expected_combined, combined_from_owned);

        // Adding an empty one should not change the outcome
        let mut aabbs = vec![
            Aabb {
                min: Vec3::new(-10.0, -10.0, -10.0),
                max: Vec3::new(-9.0, -9.0, -9.0),
            },
            Aabb {
                min: Vec3::new(1.0, 1.0, 1.0),
                max: Vec3::new(100.0, 100.0, 100.0),
            },
            Aabb::empty()
        ];
        // first, create from iterator over references
        let combined : Aabb = aabbs.iter()
            .collect();

        assert_eq!(expected_combined, combined);

        // But adding infinity should change a lot
        aabbs.push(Aabb::infinite());
        let expected_combined = Aabb::infinite();

        let combined : Aabb = aabbs.iter()
            .collect();

        assert_eq!(expected_combined, combined);
    }

    #[test]
    fn aabb_from_points_empty() {
        let aabb = Aabb::from_points(iter::empty::<Vec3>());
        assert!(aabb.is_empty());

        let aabb = Aabb::from_points(iter::empty::<&Vec3>());
        assert!(aabb.is_empty());
    }

    #[test]
    fn aabb_from_single_point() {
        let point = Vec3::new(1.0, 2.0, 3.0);
        let aabb = Aabb::from_points(iter::once(point));

        assert_eq!(aabb.min, point, "Built AABB from single point {:?} and expected min to be equal, but was {:?}", point, aabb.min);
        assert_eq!(aabb.max, point, "Built AABB from single point {:?} and expected max to be equal, but was {:?}", point, aabb.max);
        assert!(aabb.is_empty());
    }

    #[test]
    fn aabb_from_points_triangle() {
        let aabb = Aabb::from_points(vec![
            Vec3::new(-0.5, -0.5, 1.0),
            Vec3::new(0.5, -0.5, 1.0),
            Vec3::new(0.0, 0.5, -1.0)
        ]);

        assert_eq!(aabb.min, Vec3::new(-0.5, -0.5, -1.0));
        assert_eq!(aabb.max, Vec3::new(0.5, 0.5, 1.0));
    }

    #[test]
    fn volume() {
        let aabb = Aabb {
            min: Vec3::new(0.0, 0.0, 0.0),
            max: Vec3::new(1.0, 1.0, 1.0),
        };

        assert_eq!(aabb.volume(), 1.0);
    }

    #[test]
    fn empty_volumes() {
        // Different versions of empty Aabb should all report
        // zero.
        // If min is greater or equal to max in any dimension,
        // the volume should be zero.
        let empty = Aabb::empty();

        // Also test that multiple such cases do not
        // cancel each other out e.g. x*y*z when
        // x = -10, y = -10, z = 10 is positive
        let one_dimension_min_gt_max = Aabb {
            min: Vec3::new(1.1, 0.0, 0.0),
            max: Vec3::new(1.0, 1.0, 1.0)
        };
        let two_dimension_min_gt_max = Aabb {
            min: Vec3::new(1.1, 1.1, 0.0),
            max: Vec3::new(1.0, 1.0, 1.0)
        };
        let three_dimension_min_gt_max = Aabb {
            min: Vec3::new(1.1, 1.1, 1.1),
            max: Vec3::new(1.0, 1.0, 1.0)
        };

        assert_eq!(empty.volume(), 0.0);
        assert_eq!(one_dimension_min_gt_max.volume(), 0.0);
        assert_eq!(two_dimension_min_gt_max.volume(), 0.0);
        assert_eq!(three_dimension_min_gt_max.volume(), 0.0);
    }

    #[test]
    fn inside() {
        let aabb = Aabb::from_points(vec![
            Vec3::new(-0.5, -0.5, 1.0),
            Vec3::new(0.5, -0.5, 1.0),
            Vec3::new(0.0, 0.5, -1.0)
        ]);

        let point = Vec3::new(0.0, 0.0, 0.0);
        assert!(aabb.is_point_inside(point));

        let point = Vec3::new(10.0, 0.0, 0.0);
        assert!(!aabb.is_point_inside(point));

        let other_aabb = Aabb {
            min: Vec3::new(-0.1, -0.1, -0.1),
            max: Vec3::new(0.1, 0.1, 0.1),
        };
        assert!(aabb.is_aabb_inside(&other_aabb));

        let other_aabb = Aabb {
            min: Vec3::new(99.9, -0.1, -0.1),
            max: Vec3::new(100.1, 0.1, 0.1),
        };
        assert!(!aabb.is_aabb_inside(&other_aabb));
    }

    /// Tests the special case with only three points that lie on the same axis.
    /// The aabb should still contain all points, despite having zero volume
    #[test]
    fn axis_aligned_triangle() {
        let p1 = Vec3::new(-0.5, -0.5, 1.0);
        let p2 = Vec3::new(0.5, -0.5, 1.0);
        let p3 = Vec3::new(0.0, 0.5, 1.0);

        let aabb = Aabb::from_points([p1, p2, p3].iter());

        assert!(aabb.is_point_inside(p1));
        assert!(aabb.is_point_inside(p2));
        assert!(aabb.is_point_inside(p3));
        assert!(aabb.is_empty());
        assert_eq!(aabb.volume(), 0.0);
    }

    #[test]
    fn aabb_line_intersection() {
        // Given one aabb around the origin
        let aabb = Aabb::from_points(vec![
            Vec3::new(-0.5, -0.5, 1.0),
            Vec3::new(0.5, -0.5, 1.0),
            Vec3::new(0.0, 0.5, -1.0)
        ]);

        // Line above the origin, with direction facing up, hit
        let line_origin = Vec3::new(0.0, 10.0, 0.0);
        let line_direction = Vec3::new(0.0, -1.0, 0.0);
        let intersection_params = aabb.line_intersection_min_max_parameters(line_origin, line_direction);
        assert!(intersection_params.is_some(), "Line above the origin, with direction facing up, should hit an aabb below");
        let (min_t, max_t) = intersection_params.unwrap();
        assert!(min_t > 0.0 && max_t > 0.0 && max_t > min_t);

        // Line above the origin, with direction facing down should also hit an aabb below, this is not a ray
        let line_origin = Vec3::new(0.0, 10.0, 0.0);
        let line_direction = Vec3::new(0.0, 1.0, 0.0);
        let intersection_params = aabb.line_intersection_min_max_parameters(line_origin, line_direction);
        assert!(intersection_params.is_some(), "Line above the origin, with direction facing down should also hit an aabb below, this is not a ray");
        let (min_t, max_t) = intersection_params.unwrap();
        assert!(min_t < 0.0 && max_t < 0.0 && max_t > min_t);

        // Line above the origin, facing right, should miss
        let line_origin = Vec3::new(0.0, 10.0, 0.0);
        let line_direction = Vec3::new(1.0, 0.0, 0.0);
        let intersection_params = aabb.line_intersection_min_max_parameters(line_origin, line_direction);
        assert!(intersection_params.is_none(), "Line above the origin, facing right, should miss");

        // Line above the origin, facing towards Z, should miss
        let line_origin = Vec3::new(0.0, 10.0, 0.0);
        let line_direction = Vec3::new(0.0, 0.0, 1.0);
        let intersection_params = aabb.line_intersection_min_max_parameters(line_origin, line_direction);
        assert!(intersection_params.is_none(), "Line above the origin, facing right, should miss");

    }

    #[test]
    fn aabb_intersects_ray_inside() {
        // Given one aabb around the origin
        let aabb = Aabb::from_points(vec![
            Vec3::new(-0.5, -0.5, -1.0),
            Vec3::new(0.5, 0.5, 1.0)
        ]);

        let ray_origin = Vec3::new(0.0, 0.0, 0.0);
        let ray_direction = Vec3::new(1.0, 1.0, 0.0);
        assert!(aabb.is_point_inside(ray_origin));
        assert!(aabb.intersects_ray(ray_origin, ray_direction), "No matter the direction, a ray originating inside an aabb should always hit it");
        let ray_direction = Vec3::new(-1.0, 0.0, 1.0);
        assert!(aabb.is_point_inside(ray_origin));
        assert!(aabb.intersects_ray(ray_origin, ray_direction), "No matter the direction, a ray originating inside an aabb should always hit it");
        let ray_direction = Vec3::new(0.0, 1.0, 1.0);
        assert!(aabb.is_point_inside(ray_origin));
        assert!(aabb.intersects_ray(ray_origin, ray_direction), "No matter the direction, a ray originating inside an aabb should always hit it");
    }

    #[test]
    fn aabb_intersects_ray_hit_from_outside() {
        // Given one aabb around the origin
        let aabb = Aabb::from_points(vec![
            Vec3::new(-0.5, -0.5, -1.0),
            Vec3::new(0.5, 0.5, 1.0)
        ]);

        let ray_origin = Vec3::new(10.0, 0.0, 0.0);
        let ray_direction = Vec3::new(-1.0, 0.0, 0.0);
        assert!(aabb.intersects_ray(ray_origin, ray_direction), "Ray shot from the right to an abb on the left should hit");
    }

    #[test]
    fn aabb_intersects_miss() {
        // Given one aabb around the origin
        let aabb = Aabb::from_points(vec![
            Vec3::new(-0.5, -0.5, -1.0),
            Vec3::new(0.5, 0.5, 1.0)
        ]);

        let ray_origin = Vec3::new(10.0, 0.0, 0.0);
        let ray_direction = Vec3::new(1.0, 0.0, 0.0);
        assert!(!aabb.intersects_ray(ray_origin, ray_direction), "Ray shot from the right to the right should miss an aabb on the left");
    }
}
