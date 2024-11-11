use nalgebra::{SimdPartialOrd, Vector3};

pub struct Triangle {
    pub vertex0: Vector3<f32>,
    pub vertex1: Vector3<f32>,
    pub vertex2: Vector3<f32>,
}
impl Triangle {
    /// ```
    /// use nalgebra::vector;
    /// use cvoxel::{Aabb, Triangle};
    /// let triangle = Triangle {
    ///     vertex0: vector![0.0, -3.0, 2.0],
    ///     vertex1: vector![-1.0, 2.0, 1.0],
    ///     vertex2: vector![1.0, -1.0, 3.0],
    /// };
    /// let aabb = triangle.aabb();
    /// assert_eq!(aabb, Aabb::new(vector![-1.0, -3.0, 1.0], vector![1.0, 2.0, 3.0]));
    /// ```
    pub fn aabb(&self) -> Aabb {
        let min = self.vertex0.simd_min(self.vertex1).simd_min(self.vertex2);
        let max = self.vertex0.simd_max(self.vertex1).simd_max(self.vertex2);
        Aabb {
            max, min
        }
    }
}

#[derive(PartialEq, Debug)]
pub struct Aabb {
    pub min: Vector3<f32>,
    pub max: Vector3<f32>,
}

impl Aabb {
    pub fn new(min: Vector3<f32>, max: Vector3<f32>) -> Self {
        Self {
            min, max
        }
    }
}

// pub fn voxelize(mesh: &[Triangle]) {

// }