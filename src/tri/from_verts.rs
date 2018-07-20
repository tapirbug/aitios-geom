use super::Triangle;

/// Implemented by triangles that can be built from owned vertices.
pub trait FromVertices: Triangle {
    fn new(v0: Self::Vertex, v1: Self::Vertex, v2: Self::Vertex) -> Self;

    /// Flips the first two vertices. This operation inverts the winding order
    /// and flips the geometric normal.
    fn flip(&self) -> Self
    where
        Self: Sized,
        Self::Vertex: Clone,
    {
        let (a, b, c) = self.vertices();
        Self::new(b.clone(), a.clone(), c.clone())
    }
}
