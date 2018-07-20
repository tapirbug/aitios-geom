use super::tri::Triangle;

/// Iterates over shared references to triangle vertices.
#[derive(Debug, Copy, Clone)]
pub struct Iter<'a, T: Triangle + ?Sized + 'a> {
    tri: &'a T,
    next_idx: usize,
}

impl<'a, T: Triangle + ?Sized> Iter<'a, T> {
    pub fn new(triangle: &'a T) -> Self {
        Iter {
            tri: triangle,
            next_idx: 0,
        }
    }
}

impl<'a, T: Triangle + 'a> Iterator for Iter<'a, T> {
    type Item = &'a T::Vertex;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        match self.tri.vertices() {
            (a, _, _) if self.next_idx == 0 => {
                self.next_idx += 1;
                Some(a)
            }
            (_, b, _) if self.next_idx == 1 => {
                self.next_idx += 1;
                Some(b)
            }
            (_, _, c) if self.next_idx == 2 => {
                self.next_idx += 1;
                Some(c)
            }
            _ => None,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use linalg::Vec3;
    use tri::FromVertices;
    use tri::TupleTriangle as Tri;

    #[test]
    fn test_iter() {
        let vertex0 = Vec3::new(-1.0, -1.0, 0.0);
        let vertex1 = Vec3::new(1.0, -1.0, 0.0);
        let vertex2 = Vec3::new(0.0, 1.0, 0.0);
        let tri = Tri::new(vertex0, vertex1, vertex2);

        let mut iter = tri.iter();
        assert_eq!(vertex0, *iter.next().unwrap());
        assert_eq!(vertex1, *iter.next().unwrap());
        assert_eq!(vertex2, *iter.next().unwrap());
    }
}
