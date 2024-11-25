use nalgebra::{vector, Isometry3, Point2, Point3, Translation, UnitQuaternion, Vector2, Vector3};
use num_enum::TryFromPrimitive;
use std::ops::{Add, BitAnd, BitOr, Div, Mul, Sub};

pub mod debug;

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
impl Add<Vector3<f32>> for Aabb {
    type Output = Aabb;
    fn add(self, rhs: Vector3<f32>) -> Self::Output {
        Aabb {
            min: self.min + rhs,
            max: self.max + rhs
        }
    }
}
impl Mul<f32> for Aabb {
    type Output = Aabb;
    fn mul(self, rhs: f32) -> Self::Output {
        Aabb {
            min: self.min * rhs,
            max: self.max * rhs
        }
    }
}
impl Div<f32> for Aabb {
    type Output = Aabb;
    fn div(self, rhs: f32) -> Self::Output {
        let inv = 1.0 / rhs;
        self * inv
    }
}
impl Sub<Vector3<f32>> for Aabb {
    type Output = Aabb;
    fn sub(self, rhs: Vector3<f32>) -> Self::Output {
        Aabb {
            min: self.min - rhs,
            max: self.max - rhs
        }
    }
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

#[derive(Default)]
pub struct CVoxels {
    pub shape: Vector3<usize>,
    /// The position part is the center of the voxel object.
    pub transform: Isometry3<f32>,
    pub dx: f32,
    pub data: Vec<CVoxelType>,
}
impl CVoxels {
    /// Aabb in world space
    pub fn aabb(&self) -> Aabb {
        let ws_half_extents =
            self.transform.rotation.to_rotation_matrix().matrix().abs() * 0.5 * self.size();
        Aabb {
            min: (self.transform.translation.vector - ws_half_extents).into(),
            max: (self.transform.translation.vector + ws_half_extents).into(),
        }
    }

    /// Aabb in local space
    pub fn local_aabb(&self) -> Aabb {
        let half_extents = 0.5 * self.size();
        Aabb {
            min: Point3::from(-half_extents),
            max: Point3::from(half_extents)
        }
    }

    /// self.shape.x * self.shape.y
    pub fn area(&self) -> usize {
        self.shape.x * self.shape.y
    }

    /// self.shape.cast::<f32>() * self.dx
    pub fn size(&self) -> Vector3<f32> {
        self.shape.cast::<f32>() * self.dx
    }

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

    /// Voxelize a mesh. Returns None if the length of the indices can not be devided by 3.
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

    /// Voxelize a mesh. Returns None if the length of the vertices can not be devided by 3.
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
        Some(Self::from_triangles(&triangles, dx))
    }

    pub fn from_triangles(triangles: &[Triangle], dx: f32) -> CVoxels {
        // object attributes
        let total_aabb = if let Some(aabb) = triangles
            .iter()
            .map(|tri| tri.aabb())
            .reduce(|acc, cur| acc | cur)
        {
            aabb
        } else {
            return Self::default();
        };
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

        cvoxels
    }

    pub fn intersection_aabb(&self, rhs: &CVoxels) -> Option<Aabb> {
        let mut vertices: [Point3<f32>; 8] = [
            Point3::new(-0.5, -0.5, -0.5),
            Point3::new(0.5, -0.5, -0.5),
            Point3::new(0.5, 0.5, -0.5),
            Point3::new(-0.5, 0.5, -0.5),
            Point3::new(-0.5, -0.5, 0.5),
            Point3::new(0.5, -0.5, 0.5),
            Point3::new(0.5, 0.5, 0.5),
            Point3::new(-0.5, 0.5, 0.5),
        ];
        let size = rhs.size();
        let transform = self.transform.inverse() * rhs.transform;
        for vertex in &mut vertices {
            *vertex = transform * Point3::from(vertex.coords.component_mul(&size));
        }
        const LINE_INDICES: [(u8, u8); 12] = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),

            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),

            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),
        ];

        let self_aabb = self.local_aabb();

        let mut lines = Vec::with_capacity(12);
        for i in 0..LINE_INDICES.len() {
            let (i1, i2) = LINE_INDICES[i];
            lines.push([vertices[i1 as usize], vertices[i2 as usize]]);
        }

        fn arg_min_max(x1: f32, x2: f32) -> ((usize, f32), (usize, f32)) {
            if x1 < x2 {
                ((0, x1), (1, x2))
            } else {
                ((1, x2), (0, x1))
            }
        }
        fn cross_2d(v1: Vector2<f32>, v2: Vector2<f32>) -> f32 {
            v1.x * v2.y - v1.y * v2.x
        }
        fn convex_sort(points: &mut [Point2<f32>]) {
            if points.len() <= 3 {
                return;
            }
            let ref_point = points[0];
            points[1..].sort_unstable_by(|a, b| {
                let cw = cross_2d(a - ref_point, b - ref_point) > 0.0;
                if cw {
                    std::cmp::Ordering::Less
                } else {
                    std::cmp::Ordering::Greater
                }
            });
        }

        // clipping
        // x max
        let mut crossings = Vec::new();
        let mut i = 0;
        while i < lines.len() {
            let ((min_i, min), (max_i, max)) = arg_min_max(lines[i][0].x, lines[i][1].x);
            if min > self_aabb.max.x {
                lines.swap_remove(i);
                continue;
            }
            if max > self_aabb.max.x {
                let crossing = lines[i][min_i].lerp(&lines[i][max_i], (self_aabb.max.x - min) / (max - min));
                lines[i][max_i] = crossing;
                crossings.push(Point2::new(crossing.y, crossing.z));
            }
            i += 1;
        }
        convex_sort(&mut crossings);
        for i in 1..crossings.len() {
            let p1 = crossings[i - 1];
            let p2 = crossings[i];
            lines.push([
                Point3::new(self_aabb.max.x, p1.x, p1.y),
                Point3::new(self_aabb.max.x, p2.x, p2.y),
            ])
        }
        // x min
        crossings.clear();
        let mut i = 0;
        while i < lines.len() {
            let ((min_i, min), (max_i, max)) = arg_min_max(lines[i][0].x, lines[i][1].x);
            if max < self_aabb.min.x {
                lines.swap_remove(i);
                continue;
            }
            if min < self_aabb.min.x {
                let crossing = lines[i][min_i].lerp(&lines[i][max_i], (self_aabb.min.x - min) / (max - min));
                lines[i][min_i] = crossing;
                crossings.push(Point2::new(crossing.y, crossing.z));
            }
            i += 1;
        }
        convex_sort(&mut crossings);
        for i in 1..crossings.len() {
            let p1 = crossings[i - 1];
            let p2 = crossings[i];
            lines.push([
                Point3::new(self_aabb.min.x, p1.x, p1.y),
                Point3::new(self_aabb.min.x, p2.x, p2.y),
            ])
        }
        // y max
        crossings.clear();
        let mut i = 0;
        while i < lines.len() {
            let ((min_i, min), (max_i, max)) = arg_min_max(lines[i][0].y, lines[i][1].y);
            if min > self_aabb.max.y {
                lines.swap_remove(i);
                continue;
            }
            if max > self_aabb.max.y {
                let crossing = lines[i][min_i].lerp(&lines[i][max_i], (self_aabb.max.y - min) / (max - min));
                lines[i][max_i] = crossing;
                crossings.push(Point2::new(crossing.x, crossing.z));
            }
            i += 1;
        }
        convex_sort(&mut crossings);
        for i in 1..crossings.len() {
            let p1 = crossings[i - 1];
            let p2 = crossings[i];
            lines.push([
                Point3::new(p1.x, self_aabb.max.y, p1.y),
                Point3::new(p2.x, self_aabb.max.y, p2.y),
            ])
        }
        // y min
        crossings.clear();
        let mut i = 0;
        while i < lines.len() {
            let ((min_i, min), (max_i, max)) = arg_min_max(lines[i][0].y, lines[i][1].y);
            if max < self_aabb.min.y {
                lines.swap_remove(i);
                continue;
            }
            if min < self_aabb.min.y {
                let crossing = lines[i][min_i].lerp(&lines[i][max_i], (self_aabb.min.y - min) / (max - min));
                lines[i][min_i] = crossing;
                crossings.push(Point2::new(crossing.x, crossing.z));
            }
            i += 1;
        }
        convex_sort(&mut crossings);
        for i in 1..crossings.len() {
            let p1 = crossings[i - 1];
            let p2 = crossings[i];
            lines.push([
                Point3::new(p1.x, self_aabb.min.y, p1.y),
                Point3::new(p2.x, self_aabb.min.y, p2.y),
            ])
        }

        // z max
        crossings.clear();
        let mut i = 0;
        while i < lines.len() {
            let ((min_i, min), (max_i, max)) = arg_min_max(lines[i][0].z, lines[i][1].z);
            if min > self_aabb.max.z {
                lines.swap_remove(i);
                continue;
            }
            if max > self_aabb.max.z {
                let crossing = lines[i][min_i].lerp(&lines[i][max_i], (self_aabb.max.z - min) / (max - min));
                lines[i][max_i] = crossing;
                crossings.push(Point2::new(crossing.x, crossing.y));
            }
            i += 1;
        }
        convex_sort(&mut crossings);
        for i in 1..crossings.len() {
            let p1 = crossings[i - 1];
            let p2 = crossings[i];
            lines.push([
                Point3::new(p1.x, p1.y, self_aabb.max.z),
                Point3::new(p2.x, p2.y, self_aabb.max.z),
            ])
        }
        // z min
        crossings.clear();
        let mut i = 0;
        while i < lines.len() {
            let ((min_i, min), (max_i, max)) = arg_min_max(lines[i][0].z, lines[i][1].z);
            if max < self_aabb.min.z {
                lines.swap_remove(i);
                continue;
            }
            if min < self_aabb.min.z {
                let crossing = lines[i][min_i].lerp(&lines[i][max_i], (self_aabb.min.z - min) / (max - min));
                lines[i][min_i] = crossing;
                crossings.push(Point2::new(crossing.x, crossing.y));
            }
            i += 1;
        }
        convex_sort(&mut crossings);
        for i in 1..crossings.len() {
            let p1 = crossings[i - 1];
            let p2 = crossings[i];
            lines.push([
                Point3::new(p1.x, p1.y, self_aabb.min.z),
                Point3::new(p2.x, p2.y, self_aabb.min.z),
            ])
        }

        const MAX_POINT: Point3<f32> = Point3::new(f32::MAX, f32::MAX, f32::MAX);
        const MIN_POINT: Point3<f32> = Point3::new(f32::MIN, f32::MIN, f32::MIN);
        let mut min = MAX_POINT;
        let mut max = MIN_POINT;

        for line in lines {
            for point in line {
                min = point.inf(&min);
                max = point.sup(&max);
            }
        }

        if min != MAX_POINT {
            Some(Aabb {
                min, max
            })
        } else {
            None
        }
    }

    /// This function does not use world space aabb for acceleration, but it uses `intersection_aabb()` to only consider the intersection part.
    /// Returns the first voxel pair found intersecting.
    pub fn intersected(&self, rhs: &CVoxels) -> Option<(usize, usize)> {
        let intersection_aabb = if let Some(aabb) = self.intersection_aabb(rhs) {
            aabb
        } else {
            return None;
        };

        let half_size = self.size() * 0.5;

        let intersection_part = (intersection_aabb + half_size) / self.dx;

        #[inline]
        fn component_floor(p: &Point3<f32>) -> Point3<usize> {
            Point3::new(p.x as usize, p.y as usize, p.z as usize)
        }
        #[inline]
        fn component_ceil(p: &Point3<f32>) -> Point3<usize> {
            Point3::new(p.x.ceil() as usize, p.y.ceil() as usize, p.z.ceil() as usize)
        }

        let start = component_floor(&intersection_part.min);
        let end = component_ceil(&intersection_part.max).inf(&self.shape.into());

        use std::ops::Range;
        #[inline]
        fn axis_range(dimensionless: f32, len: usize) -> Range<usize> {
            if dimensionless < -0.5 {
                Range { start: 0, end: 0 }
            } else if dimensionless < 0.5 {
                Range { start: 0, end: 1}
            } else if dimensionless < len as f32 - 0.5 {
                let start = dimensionless.round() as usize - 1;
                Range { start, end: start + 2 }
            } else if dimensionless < len as f32 + 0.5 {
                Range { start: len - 1, end: len}
            } else {
                Range { start: len, end: len }
            }
        }
        #[inline]
        fn axis_ranges(dimensionless: &Point3<f32>, shape: &Vector3<usize>) -> [Range<usize>; 3] {
            [
                axis_range(dimensionless.x, shape.x),
                axis_range(dimensionless.y, shape.y),
                axis_range(dimensionless.z, shape.z),
            ]
        }

        let to_rhs_transform = rhs.transform.inverse() * self.transform;
        let rhs_half_size = rhs.size() * 0.5;
        let rhs_inv_dx = 1.0 / rhs.dx;
        let rhs_area = rhs.area();

        let area = self.area();
        for k in start.z..end.z {
            let z_base = k * area;
            let z = (k as f32 + 0.5) * self.dx - half_size.z;
            for j in start.y..end.y {
                let yz_base = z_base + j * self.shape.x;
                let y = (j as f32 + 0.5) * self.dx - half_size.y;
                for i in start.x..end.x {
                    let index = yz_base + i;
                    let x = (i as f32 + 0.5) * self.dx - half_size.x;
                    if self.data[index] == CVoxelType::Air {
                        continue;
                    }
                    let coord = Point3::new(x, y, z);
                    let coord_in_rhs = to_rhs_transform * coord;
                    let unit_coord = (coord_in_rhs + rhs_half_size) * rhs_inv_dx;
                    let [range_x, range_y, range_z] = axis_ranges(&unit_coord, &rhs.shape);
                    
                    for rk in range_z {
                        let rz_base = rk * rhs_area;
                        for rj in range_y.clone() {
                            let rzy_base = rz_base + rj * rhs.shape.x;
                            for ri in range_x.clone() {
                                let rindex = rzy_base + ri;
                                if rhs.data[rindex] != CVoxelType::Air {
                                    return Some((index, rindex));
                                }
                            }
                        }
                    }
                }
            }
        }
        None
    }
}