use nalgebra::{vector, Isometry3, Point3, Translation, UnitQuaternion, Vector3};
use num_enum::TryFromPrimitive;
use std::ops::{BitAnd, BitOr};

mod debug;

pub struct Triangle {
    pub vertex0: Point3<f32>,
    pub vertex1: Point3<f32>,
    pub vertex2: Point3<f32>,
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

    /// Returns the traveling distance of the ray if there is an intersection,
    /// including traveling backwards resulting in negative distance.
    pub fn intersect(&self, ray: &Ray) -> Option<f32> {
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
        Some(f * edge2.dot(&q))
    }
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

pub struct Ray {
    pub origin: Point3<f32>,
    pub direction: Vector3<f32>,
}
impl Ray {
    pub fn new(origin: Point3<f32>, direction: Vector3<f32>) -> Self {
        Self { origin, direction }
    }
}

/// Axis-Aligned Bounding Box.
#[derive(PartialEq, Debug, Clone, Copy)]
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
            max: self.max.sup(&other.max),
        }
    }
    pub fn intersection(&self, other: &Aabb) -> Aabb {
        Aabb {
            min: self.min.sup(&other.min),
            max: self.max.inf(&other.max),
        }
    }
    pub fn intersect(&self, other: &Aabb) -> bool {
        self.intersection(other).is_empty()
    }
    pub fn is_empty(&self) -> bool {
        self.min.x > self.max.x || self.min.y > self.max.y
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
impl BitAnd for Aabb {
    type Output = Aabb;
    fn bitand(self, rhs: Self) -> Self::Output {
        self.intersection(&rhs)
    }
}

/// Oriented Bounding Box. Untested.
pub struct Obb {
    pub transform: Isometry3<f32>,
    pub half_size: Vector3<f32>,
}
impl Obb {
    pub fn aabb(&self) -> Aabb {
        let ws_half_extents = self
            .transform
            .rotation
            .to_rotation_matrix()
            .into_inner()
            .abs()
            * self.half_size;
        let center = Point3::from(self.transform.translation.vector);
        let min = center - ws_half_extents;
        let max = center + ws_half_extents;
        Aabb { min, max }
    }
    pub fn vertices(&self) -> [Point3<f32>; 8] {
        let t = self.transform;
        let e = self.half_size;
        [
            t * Point3::new(e.x, e.y, e.z),
            t * Point3::new(-e.x, e.y, e.z),
            t * Point3::new(e.x, -e.y, e.z),
            t * Point3::new(-e.x, -e.y, e.z),
            t * Point3::new(e.x, e.y, -e.z),
            t * Point3::new(-e.x, e.y, -e.z),
            t * Point3::new(e.x, -e.y, -e.z),
            t * Point3::new(-e.x, -e.y, -e.z),
        ]
    }
    pub fn intersect(&self, other: &Obb) -> bool {
        let r1 = self.transform.rotation.to_rotation_matrix();
        let m1 = r1.matrix();
        let x1 = m1.column(0).clone_owned();
        let y1 = m1.column(1).clone_owned();
        let z1 = m1.column(2).clone_owned();
        let v1 = self.vertices();

        let r2 = self.transform.rotation.to_rotation_matrix();
        let m2 = r2.matrix();
        let x2 = m2.column(0).clone_owned();
        let y2 = m2.column(1).clone_owned();
        let z2 = m2.column(2).clone_owned();
        let v2 = other.vertices();

        sat_intersect(&v1, &v2, x1)
            && sat_intersect(&v1, &v2, x2)
            && sat_intersect(&v1, &v2, y1)
            && sat_intersect(&v1, &v2, y2)
            && sat_intersect(&v1, &v2, z1)
            && sat_intersect(&v1, &v2, z2)
            && sat_intersect(&v1, &v2, x1.cross(&x2))
            && sat_intersect(&v1, &v2, x1.cross(&y2))
            && sat_intersect(&v1, &v2, x1.cross(&z2))
            && sat_intersect(&v1, &v2, y1.cross(&x2))
            && sat_intersect(&v1, &v2, y1.cross(&y2))
            && sat_intersect(&v1, &v2, y1.cross(&z2))
            && sat_intersect(&v1, &v2, z1.cross(&x2))
            && sat_intersect(&v1, &v2, z1.cross(&y2))
            && sat_intersect(&v1, &v2, z1.cross(&z2))
    }
}

/// Separating Axis Theorem.
pub fn sat_intersect(
    vertices1: &[Point3<f32>],
    vertices2: &[Point3<f32>],
    axis: Vector3<f32>,
) -> bool {
    if vertices1.is_empty() || vertices2.is_empty() {
        return false;
    }
    let fold_min_max = |(min, max): (f32, f32), acc: f32| (acc.min(min), acc.max(max));
    let (min1, max1) = vertices1
        .iter()
        .map(|v| v.coords.dot(&axis))
        .fold((f32::MAX, f32::MIN), fold_min_max);
    let (min2, max2) = vertices2
        .iter()
        .map(|v| v.coords.dot(&axis))
        .fold((f32::MAX, f32::MIN), fold_min_max);
    min1 < max2 && max1 > min2
}

#[derive(PartialEq, Eq, num_enum::TryFromPrimitive, Clone, Copy)]
#[repr(u8)]
pub enum CVoxelType {
    Body,
    Face,
    Edge,
    Corner,
    Air,
}

pub struct CVoxels {
    pub shape: Vector3<usize>,
    /// The position part is the center of the voxel object.
    pub transform: Isometry3<f32>,
    pub dx: f32,
    pub data: Vec<CVoxelType>,
}
impl CVoxels {
    /// Regenerate VoxelCType.
    pub fn regenerate_type(&mut self) {
        let area = self.area();
        for k in 0..self.shape.z {
            let k_part = k * area;
            for j in 0..self.shape.y {
                let j_part = j * self.shape.x;
                for i in 0..self.shape.x {
                    let coord = k_part + j_part + i;

                    if self.data[coord] == CVoxelType::Air {
                        continue;
                    }

                    let mut air_axis_cound: u8 = 0;

                    if k == 0
                        || k == self.shape.z - 1
                        || self.data[coord - area] == CVoxelType::Air
                        || self.data[coord + area] == CVoxelType::Air
                    {
                        air_axis_cound += 1;
                    };
                    if j == 0
                        || j == self.shape.y - 1
                        || self.data[coord - self.shape.x] == CVoxelType::Air
                        || self.data[coord + self.shape.x] == CVoxelType::Air
                    {
                        air_axis_cound += 1;
                    };
                    if i == 0
                        || i == self.shape.x - 1
                        || self.data[coord - 1] == CVoxelType::Air
                        || self.data[coord + 1] == CVoxelType::Air
                    {
                        air_axis_cound += 1;
                    };

                    self.data[coord] = CVoxelType::try_from_primitive(air_axis_cound).unwrap();
                }
            }
        }
    }

    pub fn obb(&self) -> Obb {
        Obb {
            half_size: self.size() * 0.5,
            transform: self.transform,
        }
    }

    pub fn aabb(&self) -> Aabb {
        self.obb().aabb()
    }

    /// self.shape.x * self.shape.y
    pub fn area(&self) -> usize {
        self.shape.x * self.shape.y
    }

    /// self.shape.cast::<f32>() * self.dx
    pub fn size(&self) -> Vector3<f32> {
        self.shape.cast::<f32>() * self.dx
    }

    pub fn from_indexed_mesh<T: num_traits::AsPrimitive<usize>>(
        vertices: &[[f32; 3]],
        indices: &[T],
        dx: f32,
    ) -> Option<CVoxels> {
        let trimesh = indices
            .iter()
            .map(|i| vertices[i.as_()])
            .collect::<Vec<[f32; 3]>>();
        Self::from_trimesh(&trimesh, dx)
    }

    /// Voxelize a mesh. Returns None if the length of the vertices is zero or can not be devided by 3.
    pub fn from_trimesh(vertices: &[[f32; 3]], dx: f32) -> Option<CVoxels> {
        if vertices.len() % 3 != 0 {
            return None;
        }
        let triangles = vertices
            .chunks(3)
            .map(|chunk| Triangle {
                vertex0: Point3::from(chunk[0]),
                vertex1: Point3::from(chunk[1]),
                vertex2: Point3::from(chunk[2]),
            })
            .collect::<Vec<Triangle>>();
        Self::from_triangles(&triangles, dx)
    }

    /// Voxelize a mesh. Returns None if the length of the triangles is zero.
    pub fn from_triangles(triangles: &[Triangle], dx: f32) -> Option<CVoxels> {
        // object attributes
        let total_aabb = triangles
            .iter()
            .map(|tri| tri.aabb())
            .reduce(|acc, cur| acc | cur)?;
        let size = total_aabb.size();
        let shape = vector![
            (size.x / dx).ceil() as usize,
            (size.y / dx).ceil() as usize,
            (size.z / dx).ceil() as usize,
        ];
        let transform = Isometry3 {
            translation: Translation::from(total_aabb.middle().coords),
            rotation: UnitQuaternion::identity(),
        };

        // fill voxels
        let mut data = Vec::with_capacity(shape.product());
        for k in 0..shape.z {
            let z = (k as f32 + 0.5) * dx + total_aabb.min.z;
            for j in 0..shape.y {
                let y = (j as f32 + 0.5) * dx + total_aabb.min.y;
                let mut inside = false;
                for i in 0..shape.x {
                    let xo = (i as f32 - 0.5) * dx + total_aabb.min.x;
                    let ray = Ray {
                        origin: Point3::new(xo, y, z),
                        direction: Vector3::new(1.0, 0.0, 0.0),
                    };
                    let mut max_t = 0.0;
                    for i in 0..triangles.len() {
                        if let Some(t) = triangles[i].intersect(&ray) {
                            if t > max_t && t <= dx {
                                max_t = t;
                                inside = triangles[i].cal_norm_cw().x < 0.0;
                            }
                        }
                    }
                    if inside {
                        data.push(CVoxelType::Body);
                    } else {
                        data.push(CVoxelType::Air)
                    }
                }
            }
        }

        // generate types
        let mut cvoxels = CVoxels {
            shape,
            transform,
            dx,
            data,
        };
        cvoxels.regenerate_type();

        Some(cvoxels)
    }
}
