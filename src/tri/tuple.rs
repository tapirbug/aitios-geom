use super::tri::{Triangle, FromVertices};
use ::vtx::Position;

/// Implementation of [`Triangle`](trait.Triangle.html) that owns its vertices and stores them
/// in a tuple.
#[derive(Debug, Copy, Clone)]
pub struct TupleTriangle<V>(pub V, pub V, pub V);

impl<V : Position> FromVertices for TupleTriangle<V> {
    fn new(v0: V, v1: V, v2: V) -> Self {
        TupleTriangle(v0, v1, v2)
    }
}

impl<V : Position> Triangle for TupleTriangle<V> {
    type Vertex = V;

    fn vertices(&self) -> (&Self::Vertex, &Self::Vertex, &Self::Vertex) {
        let &TupleTriangle(ref v0, ref v1, ref v2) = self;
        (v0, v1, v2)
    }
}
