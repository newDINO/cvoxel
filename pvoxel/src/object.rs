use std::ops::{Deref, DerefMut};

use cvoxel::{ty_of_data, CVoxelType, CVoxels};
use nalgebra::{Matrix3, Point3, UnitQuaternion, Vector3};


#[derive(PartialEq, Eq, Clone, Copy)]
pub enum RigidType {
    Dynamic,
    Fixed,
}

/// Unlike CVoxels, the transform of a PVoxels object is around its center of mass.
pub struct PVoxels {
    pub cvoxels: CVoxels,
    pub ty: RigidType,

    /// local center of mass.
    pub local_cm: Point3<f32>,
    pub inv_mass: f32,
    /// local moment of inertia.
    pub inv_local_inertia: Matrix3<f32>,

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
                    if ty_of_data(self.data[index]) == CVoxelType::Air {
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
                    if ty_of_data(self.data[index]) == CVoxelType::Air {
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
}