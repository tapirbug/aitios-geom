//! Exposes the public interface of `cgmath` as well as some shorthands
//! for commonly used `cgmath` types.

pub use cgmath::*;

pub type Vec2 = Vector2<f32>;
pub type Vec3 = Vector3<f32>;
pub type Vec4 = Vector4<f32>;
pub type Mat3 = Matrix3<f32>;
pub type Mat4 = Matrix4<f32>;
