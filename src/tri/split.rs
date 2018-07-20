use tri::{FromVertices, InterpolateVertex};
use vtx::Position;

/// Takes a triangle and returns four newly created triangles that are
/// the resulting of splitting the input triangle into four by introducing
/// new vertices at edge midpoints
pub fn split_at_edge_midpoints<T, V>(tri: &T) -> Split<T>
where
    T: FromVertices<Vertex = V> + InterpolateVertex<Vertex = V> + Sized,
    V: Position + Clone,
{
    Split::new(tri)
}

/// Iterates over the resulting triangles when a triangle was split
/// at its edge midpoints.
///
/// Call [split_at_edge_midpoints](fn.split_at_edge_midpoints.html) to
/// obtain an iterator over sub-triangles.
pub struct Split<'a, T: 'a + InterpolateVertex> {
    next_idx: usize,
    src_tri: &'a T,
    mids: [T::Vertex; 3],
}

impl<'a, T> Split<'a, T>
where
    T: InterpolateVertex,
{
    fn new(tri: &'a T) -> Self {
        let mids = [
            tri.interpolate_vertex_bary([0.5, 0.5, 0.0]),
            tri.interpolate_vertex_bary([0.0, 0.5, 0.5]),
            tri.interpolate_vertex_bary([0.5, 0.0, 0.5]),
        ];

        let src_tri = tri;

        Split {
            next_idx: 0,
            src_tri,
            mids,
        }
    }
}

impl<'a, T, V> Iterator for Split<'a, T>
where
    T: InterpolateVertex<Vertex = V> + FromVertices<Vertex = V> + Sized,
    V: Position + Clone,
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        match self.next_idx {
            // The three triangles that have one vertex of the outer triangle and two edge midpoints
            // as a vertex
            mid_idx0 @ 0...2 => {
                self.next_idx += 1;
                let v0 = self.mids[mid_idx0].clone();
                let v1 = self
                    .src_tri
                    .iter()
                    .skip((mid_idx0 + 1) % 3)
                    .next()
                    .unwrap()
                    .clone();
                let v2 = self.mids[(mid_idx0 + 1) % 3].clone();
                Some(T::new(v0, v1, v2))
            }
            // The last triangle in the center, having all the midpoints as vertices
            3 => {
                self.next_idx += 1;
                Some(T::new(
                    self.mids[0].clone(),
                    self.mids[1].clone(),
                    self.mids[2].clone(),
                ))
            }
            // Then we are done with all four
            _ => None,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use linalg::Vec3;
    use tri::{Triangle, TupleTriangle};

    #[test]
    fn test_splitting_at_edge_midpoints() {
        let tri = TupleTriangle::new(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(100.0, 0.0, 0.0),
            Vec3::new(0.0, 100.0, 0.0),
        );

        let triangles: Vec<_> = split_at_edge_midpoints(&tri).collect();

        assert_eq!(triangles.len(), 4);

        assert!(
            triangles
                .iter()
                .any(|t| {
                    t.iter().map(|v| v.position()).all(|p| {
                        ulps_eq!(p, Vec3::new(50.0, 0.0, 0.0))
                            || ulps_eq!(p, Vec3::new(50.0, 50.0, 0.0))
                            || ulps_eq!(p, Vec3::new(0.0, 50.0, 0.0))
                    })
                }),
            "Expected a triangle in the lower left corner"
        );

        assert!(
            triangles
                .iter()
                .any(|t| {
                    t.iter().map(|v| v.position()).all(|p| {
                        ulps_eq!(p, Vec3::new(50.0, 0.0, 0.0))
                            || ulps_eq!(p, Vec3::new(100.0, 0.0, 0.0))
                            || ulps_eq!(p, Vec3::new(50.0, 50.0, 0.0))
                    })
                }),
            "Expected a triangle in the right corner"
        );

        assert!(
            triangles
                .iter()
                .any(|t| {
                    t.iter().map(|v| v.position()).all(|p| {
                        ulps_eq!(p, Vec3::new(50.0, 50.0, 0.0))
                            || ulps_eq!(p, Vec3::new(0.0, 100.0, 0.0))
                            || ulps_eq!(p, Vec3::new(0.0, 50.0, 0.0))
                    })
                }),
            "Expected a triangle in the top corner"
        );

        assert!(
            triangles
                .iter()
                .any(|t| {
                    t.iter().map(|v| v.position()).all(|p| {
                        ulps_eq!(p, Vec3::new(50.0, 50.0, 0.0))
                            || ulps_eq!(p, Vec3::new(50.0, 0.0, 0.0))
                            || ulps_eq!(p, Vec3::new(0.0, 50.0, 0.0))
                    })
                }),
            "Expected a triangle in the center"
        );
    }
}
