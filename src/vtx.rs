use ::linalg::{Vec2, Vec3};
use std::ops::{Add, Mul};

pub trait Position {
    fn position(&self) -> Vec3;
}

/// Vectors can be directly used as triangle vertices
impl Position for Vec3 {
    fn position(&self) -> Vec3 {
        *self
    }
}

pub trait Normal {
    fn normal(&self) -> Vec3;
}

pub trait Texcoords {
    fn texcoords(&self) -> Vec2;
}

/// An owned vertex.
#[derive(Debug, Clone)]
pub struct Vertex {
    pub position: Vec3,
    pub normal: Vec3,
    pub texcoords: Vec2,
}

impl Position for Vertex {
    fn position(&self) -> Vec3 {
        self.position
    }
}

impl Normal for Vertex {
    fn normal(&self) -> Vec3 {
        self.normal
    }
}

impl Texcoords for Vertex {
    fn texcoords(&self) -> Vec2 {
        self.texcoords
    }
}

impl Mul<f32> for Vertex {
    type Output = Vertex;

    fn mul(self, scalar: f32) -> Self::Output {
        Vertex {
            position: self.position * scalar,
            normal: self.normal * scalar,
            texcoords: self.texcoords * scalar,
        }
    }
}

impl Add for Vertex {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Vertex {
            position: self.position + rhs.position,
            normal: self.normal + rhs.normal,
            texcoords: self.texcoords + rhs.texcoords,
        }
    }
}

impl<'a> Mul<f32> for &'a Vertex {
    type Output = Vertex;

    fn mul(self, scalar: f32) -> Self::Output {
        self.clone() * scalar
    }
}

impl<'a> Add for &'a Vertex {
    type Output = Vertex;

    fn add(self, rhs: Self) -> Self::Output {
        self.clone() + rhs.clone()
    }
}

impl<'a> From<VertexRef<'a>> for Vertex {
    fn from(from: VertexRef<'a>) -> Self {
        Vertex {
            position: Vec3::new(from.position[0], from.position[1], from.position[2]),
            normal: Vec3::new(from.normal[0], from.normal[1], from.normal[2]),
            texcoords: Vec2::new(from.texcoords[0], from.texcoords[1])
        }
    }
}

impl<'a, 'b> From<&'b VertexRef<'a>> for Vertex {
    fn from(from: &'b VertexRef<'a>) -> Self {
        Vertex {
            position: Vec3::new(from.position[0], from.position[1], from.position[2]),
            normal: Vec3::new(from.normal[0], from.normal[1], from.normal[2]),
            texcoords: Vec2::new(from.texcoords[0], from.texcoords[1])
        }
    }
}

/// A vertex that references into an array
///
/// FIXME it turned out that VertexRef (48 bytes) is larger than Vertex (32 bytes)
///       so it would only make sense to have VertexRef if it was mutable
///       maybe make a triangles_mut() for mesh some time
#[derive(Debug, Copy, Clone)]
pub struct VertexRef<'a> {
    pub position: &'a [f32],
    pub normal: &'a [f32],
    pub texcoords: &'a [f32]
}

impl<'a> Position for VertexRef<'a> {
    fn position(&self) -> Vec3 {
        Vec3::new(self.position[0], self.position[1], self.position[2])
    }
}

impl<'a> Normal for VertexRef<'a> {
    fn normal(&self) -> Vec3 {
        Vec3::new(self.normal[0], self.normal[1], self.normal[2])
    }
}

impl<'a> Texcoords for VertexRef<'a> {
    fn texcoords(&self) -> Vec2 {
        Vec2::new(self.texcoords[0], self.texcoords[1])
    }
}

/*#[cfg(test)]
mod test {
    use super::*;
    use std::mem;

    #[test]
    fn size_test() {
        println!("VertexRef: {}", mem::size_of::<VertexRef>());
        println!("Vertex: {}", mem::size_of::<Vertex>());
        assert!(mem::size_of::<VertexRef>() < mem::size_of::<Vertex>());
    }
}*/
