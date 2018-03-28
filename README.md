# aitios-geom
Encapsulates low-level geometry code and traits for use in *aitios*.

Contains traits and structs for triangles, aabb, and basic ray intersections.

Re-exports types from the brilliant [*cgmath*](https://crates.io/crates/cgmath) crate for vectors, matrices and other linear algebra.

See the crate documentation for additional details.

## Examples
If the `Triangle` type would not providea the `area` method, here is how you could build one that works with any type implementing `Triangle`:

    use aitios_geom::Triangle;

    /// Calculate tha area of a triangle by evaluating:
    /// ||(V1 − V0) × (V2 − V0)|| / 2
    fn area<T>(triangle: T) -> f32
        where T : Triangle
    {
        let (v0, v1, v2) = triangle.positions();
        0.5 * ((v1 - v0).cross(v2 - v0)).magnitude()
    }
