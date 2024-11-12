use std::ops::BitOr;
use bitvec::prelude::*;
use nalgebra::{point, vector, Isometry3, Point3, Translation, UnitQuaternion, Vector3};

pub struct Triangle {
    pub vertex0: Point3<f32>,
    pub vertex1: Point3<f32>,
    pub vertex2: Point3<f32>,
}
pub struct Ray {
    pub origin: Point3<f32>,
    pub direction: Vector3<f32>,
}

impl Triangle {
    /// ```
    /// use cvoxel::Triangle;
    /// use nalgebra::{vector, point};
    /// let triangle = Triangle {
    ///     vertex0: point![0.0, 0.0, 0.0],
    ///     vertex1: point![0.0, 0.0, 1.0],
    ///     vertex2: point![1.0, 0.0, 0.0],
    /// };
    /// let norm = triangle.cal_norm_cw();
    /// assert_eq!(norm, vector![0.0, 1.0, 0.0]);
    /// 
    /// ```
    pub fn cal_norm_cw(&self) -> Vector3<f32> {
        let v1 = self.vertex1 - self.vertex0;
        let v2 = self.vertex2 - self.vertex0;
        v1.cross(&v2) 
    }
    pub fn intersect(&self, ray: &Ray) -> Option<Point3<f32>> {
        let epsilon = 1e-8;

        // Triangle vertices
        let v0 = self.vertex0;
        let v1 = self.vertex1;
        let v2 = self.vertex2;

        // Edges of the triangle
        let edge1 = v1 - v0;
        let edge2 = v2 - v0;

        // Begin calculating the determinant - also used to calculate u parameter
        let h = ray.direction.cross(&edge2);
        let a = edge1.dot(&h);

        // If the determinant is near zero, the ray is parallel to the triangle
        if a.abs() < epsilon {
            return None;
        }

        let f = 1.0 / a;
        let s = ray.origin - v0;
        let u = f * s.dot(&h);

        // Check if the intersection point is outside the triangle
        if u < 0.0 || u > 1.0 {
            return None;
        }

        let q = s.cross(&edge1);
        let v = f * ray.direction.dot(&q);

        // The intersection point is outside the triangle
        if v < 0.0 || u + v > 1.0 {
            return None;
        }

        // Calculate the distance along the ray to the intersection point
        let t = f * edge2.dot(&q);

        // Ray intersection
        if t > epsilon {
            // Calculate the intersection point
            let intersection_point = ray.origin + ray.direction * t;
            Some(intersection_point)
        } else {
            // This means that there is a line intersection but not a ray intersection
            None
        }
    }
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
    pub fn merge(&self, other: &Aabb) -> Aabb {
        Aabb {
            min: self.min.inf(&other.min),
            max: self.max.sup(&other.max)
        }
    }
    pub fn size(&self) -> Vector3<f32> {
        self.max - self.min
    }
    pub fn middle(&self) -> Point3<f32> {
        0.5 * (self.min + self.max.coords)
    }
}
impl BitOr for Aabb {
    type Output = Aabb;
    fn bitor(self, rhs: Self) -> Self::Output {
        self.merge(&rhs)
    }
}

pub struct Voxels<T> {
    pub shape: Vector3<usize>,
    pub transform: Isometry3<f32>,
    pub dx: f32,
    pub data: T,
}

/// Returns None if the length of the mesh is zero.
pub fn voxelize(mesh: &[Triangle], dx: f32) -> Option<Voxels<BitVec>> {
    let total_aabb = mesh
        .iter()
        .map(|tri| tri.aabb())
        .reduce(|acc, cur| acc | cur)?;
    let size = total_aabb.size();
    let shape = vector![
        (size.x / dx).ceil() as usize,
        (size.y / dx).ceil() as usize,
        (size.z / dx).ceil() as usize,
    ];
    let mut data: BitVec = BitVec::with_capacity(shape.product());
    for k in 0..shape.z {
        let z = k as f32 * dx + total_aabb.min.z;
        for j in 0..shape.y {
            let y = j as f32 * dx + total_aabb.min.y;
            for i in 0..shape.x {
                let x = i as f32 * dx + total_aabb.min.x;
            }
        }
    }
    let transform = Isometry3 {
        translation: Translation::from(total_aabb.middle().coords),
        rotation: UnitQuaternion::identity(),
    };
    Some(Voxels {
        shape, transform, dx, data
    })
}
