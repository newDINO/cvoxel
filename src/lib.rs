use nalgebra::Point3;

pub struct Triangle {
    pub vertex0: Point3<f32>,
    pub vertex1: Point3<f32>,
    pub vertex2: Point3<f32>,
}

impl Triangle {
    /// ```
    /// use nalgebra::point;
    /// use cvoxel::{Aabb, Triangle};
    /// let triangle = Triangle {
    ///     vertex0: point![0.0, -3.0, 2.0],
    ///     vertex1: point![-1.0, 2.0, 1.0],
    ///     vertex2: point![1.0, -1.0, 3.0],
    /// };
    /// let aabb = triangle.aabb();
    /// assert_eq!(aabb, Aabb::new(point![-1.0, -3.0, 1.0], point![1.0, 2.0, 3.0]));
    /// ```
    pub fn aabb(&self) -> Aabb {
        let min = self.vertex0.inf(&self.vertex1).inf(&self.vertex2);
        let max = self.vertex0.sup(&self.vertex1).sup(&self.vertex2);
        Aabb { max, min }
    }
}

#[derive(PartialEq, Debug)]
pub struct Aabb {
    pub min: Point3<f32>,
    pub max: Point3<f32>,
}

impl Aabb {
    pub fn new(min: Point3<f32>, max: Point3<f32>) -> Self {
        Self { min, max }
    }
}

// pub fn voxelize(mesh: &[Triangle]) {

// }
