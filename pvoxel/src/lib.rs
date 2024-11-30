use std::ops::{Deref, DerefMut};

use cvoxel::{CVoxelType, CVoxels};
use nalgebra::{Matrix3, Point3, UnitQuaternion, Vector3};

#[derive(PartialEq, Eq, Clone, Copy)]
pub enum RigidType {
    Dynamic,
    Fixed,
}

pub struct PVoxels {
    pub cvoxels: CVoxels,
    /// local center of mass.
    pub local_cm: Point3<f32>,
    pub inv_mass: f32,
    /// local moment of inertia.
    pub inv_local_inertia: Matrix3<f32>,
    pub ty: RigidType,
    pub vel: Vector3<f32>,
    pub ang_vel: Vector3<f32>,
}

impl Deref for PVoxels {
    type Target = CVoxels;
    fn deref(&self) -> &Self::Target {
        &self.cvoxels
    }
}

impl DerefMut for PVoxels {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.cvoxels
    }
}

impl PVoxels {
    pub fn from_cvoxels(cvoxels: CVoxels, density: f32, ty: RigidType) -> Self {
        let mut temp = Self {
            cvoxels,
            local_cm: Point3::new(0.0, 0.0, 0.0),
            inv_mass: 0.0,
            inv_local_inertia: Matrix3::zeros(),
            ty,
            vel: Vector3::zeros(),
            ang_vel: Vector3::zeros(),
        };
        temp.recal_local_cm();
        temp.recal_mass_and_inertia(density);
        temp
    }

    /// recalculate local_cm
    pub fn recal_local_cm(&mut self) {
        let mut result = Point3::new(0.0, 0.0, 0.0);
        let mut counts = 0.0;

        for k in 0..self.shape.z {
            let k_base = k * self.area;
            let z = (k as f32 + 0.5) * self.dx - self.half_size.z;
            for j in 0..self.shape.y {
                let jk_base = k_base + j * self.shape.x;
                let y = (j as f32 + 0.5) * self.dx - self.half_size.y;
                for i in 0..self.shape.x {
                    let index = jk_base + i;
                    if self.data[index] == CVoxelType::Air {
                        continue;
                    }
                    counts += 1.0;

                    let x = (i as f32 + 0.5) * self.dx - self.half_size.x;
                    result += Vector3::new(x, y, z);
                }
            }
        }
        result *= 1.0 / counts;
        self.local_cm = result;
    }

    /// recalculate inv_mass and inv_local_inertia
    pub fn recal_mass_and_inertia(&mut self, density: f32) -> Option<()> {
        let mut inertia = Matrix3::zeros();
        let mut mass = 0.0;

        for k in 0..self.shape.z {
            let k_base = k * self.area;
            let z_to_cm = (k as f32 + 0.5) * self.dx - self.half_size.z - self.local_cm.z;
            for j in 0..self.shape.y {
                let jk_base = k_base + j * self.shape.x;
                let y_to_cm = (j as f32 + 0.5) * self.dx - self.half_size.y - self.local_cm.y;
                for i in 0..self.shape.x {
                    let index = jk_base + i;
                    if self.data[index] == CVoxelType::Air {
                        continue;
                    }
                    let x_to_cm = (i as f32 + 0.5) * self.dx - self.half_size.x - self.local_cm.x;

                    let x2 = x_to_cm * x_to_cm;
                    let y2 = y_to_cm * y_to_cm;
                    let z2 = z_to_cm * z_to_cm;
                    let nxy = -x_to_cm * y_to_cm;
                    let nyz = -y_to_cm * z_to_cm;
                    let nzx = -z_to_cm * x_to_cm;

                    inertia[(0, 0)] += y2 + z2;
                    inertia[(1, 1)] += z2 + x2;
                    inertia[(2, 2)] += x2 + y2;
                    inertia[(0, 1)] += nxy;
                    inertia[(0, 2)] += nzx;
                    inertia[(1, 0)] += nxy;
                    inertia[(1, 2)] += nyz;
                    inertia[(2, 0)] += nzx;
                    inertia[(2, 1)] += nyz;

                    mass += 1.0;
                }
            }
        }

        let unit_m = density * self.dx * self.dx * self.dx;
        mass *= unit_m;
        inertia *= unit_m;

        self.inv_mass = 1.0 / mass;
        self.inv_local_inertia = inertia.try_inverse()?;
        Some(())
    }

    pub fn step_dt(&mut self, dt: f32) {
        let rot = UnitQuaternion::from_scaled_axis(self.ang_vel * dt);
        let trans = self.vel * dt;
        self.transform.rotation *= rot;
        self.transform.translation.vector += trans;
    }


    pub fn gen_contacts(&self, rhs: &PVoxels, contacts: &mut Vec<Contact>) {
        // 1. utilities
        #[inline]
        fn component_floor(p: &Point3<f32>) -> Point3<usize> {
            Point3::new(p.x as usize, p.y as usize, p.z as usize)
        }
        #[inline]
        fn component_ceil(p: &Point3<f32>) -> Point3<usize> {
            Point3::new(
                p.x.ceil() as usize,
                p.y.ceil() as usize,
                p.z.ceil() as usize,
            )
        }

        use std::ops::Range;
        #[inline]
        fn axis_range(dimensionless: f32, len: usize) -> Range<usize> {
            if dimensionless < -0.5 {
                Range { start: 0, end: 0 }
            } else if dimensionless < 0.5 {
                Range { start: 0, end: 1 }
            } else if dimensionless < len as f32 - 0.5 {
                let start = dimensionless.round() as usize - 1;
                Range {
                    start,
                    end: start + 2,
                }
            } else if dimensionless < len as f32 + 0.5 {
                Range {
                    start: len - 1,
                    end: len,
                }
            } else {
                Range {
                    start: len,
                    end: len,
                }
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

        // 2. process
        let intersection_aabb = if let Some(aabb) = self.intersection_aabb(rhs) {
            aabb
        } else {
            return;
        };
        let intersection_part = (intersection_aabb + self.half_size) / self.dx;
        let start = component_floor(&intersection_part.min);
        let end = component_ceil(&intersection_part.max).inf(&self.shape.into());

        let to_rhs_transform = rhs.transform.inverse() * self.transform;
        let rhs_inv_dx = 1.0 / rhs.dx;

        for k in start.z..end.z {
            let z_base = k * self.area;
            let z = (k as f32 + 0.5) * self.dx - self.half_size.z;
            for j in start.y..end.y {
                let yz_base = z_base + j * self.shape.x;
                let y = (j as f32 + 0.5) * self.dx - self.half_size.y;
                for i in start.x..end.x {
                    let index = yz_base + i;
                    if self.data[index] == CVoxelType::Air {
                        continue;
                    }
                    let x = (i as f32 + 0.5) * self.dx - self.half_size.x;
                    let coord = Point3::new(x, y, z);
                    let coord_in_rhs = to_rhs_transform * coord;
                    let unit_coord = (coord_in_rhs + rhs.half_size) * rhs_inv_dx;
                    let [range_x, range_y, range_z] = axis_ranges(&unit_coord, &rhs.shape);

                    for rk in range_z {
                        let rz_base = rk * rhs.area;
                        for rj in range_y.clone() {
                            let rzy_base = rz_base + rj * rhs.shape.x;
                            for ri in range_x.clone() {
                                let rindex = rzy_base + ri;
                                if rhs.data[rindex] != CVoxelType::Air {
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

pub struct Contact {
    pub point: Point3<f32>,
    pub normal: Vector3<f32>,
}