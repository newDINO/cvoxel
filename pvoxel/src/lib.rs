use std::ops::{Deref, DerefMut};

use cvoxel::{ty_of_data, CVoxelType, CVoxels, FaceDir, VoxelData};
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

#[derive(Clone, Copy, Debug)]
pub struct Contact {
    pub point: Point3<f32>,
    pub normal: Vector3<f32>,
    pub i1: usize,
    pub i2: usize,
}

pub struct PhysicsWorld {
    pub objects: Vec<PVoxels>,
    pub contacts: Vec<Contact>,
}

impl PhysicsWorld {
    pub fn step_dt(&mut self, dt: f32) {
        for obj in &mut self.objects {
            if obj.ty == RigidType::Fixed {
                continue;
            }
            obj.step_dt(dt);
        }
    }

    pub fn resolve_contacts(&mut self, iters: usize) {
        for _ in 0..iters {
            for i in 0..self.contacts.len() {
                self.resolve_one_contact(i);
            }
        }
        self.contacts.clear();
    }

    /// This will apply the impulse after it is solved.
    pub fn resolve_one_contact(&mut self, index: usize) {
        let contact = &self.contacts[index];
        let (obj1, obj2) = if contact.i1 < contact.i2 {
            let (part1, part2) = self.objects.split_at_mut(contact.i2);
            (&mut part1[contact.i1], &mut part2[0])
        } else {
            let (part1, part2) = self.objects.split_at_mut(contact.i1);
            (&mut part2[0], &mut part1[contact.i2])
        };

        fn solve_dynamic_fixed(dynamic: &mut PVoxels, contact: &Contact) {
            let cm = dynamic.transform * dynamic.local_cm;
            let r = contact.point - cm;

            let vp = dynamic.vel + dynamic.ang_vel.cross(&r);
            let vpn = vp.dot(&contact.normal);
            if vpn > 0.0 {
                return;
            }
            let dv = -2.0 * 0.7 * vpn;

            let rot = dynamic.transform.rotation.to_rotation_matrix();
            let inv_inertia = rot * dynamic.inv_local_inertia * rot.transpose();

            let cross_n = r.cross(&contact.normal);

            let m_eff_n = 1.0 / (dynamic.inv_mass + cross_n.dot(&(inv_inertia * cross_n)));
            let lambda_n = dv * m_eff_n;

            // friction
            if 1.0 - contact.normal.dot(&vp).abs() < 1e-3 {
                let impluse = lambda_n * contact.normal;
                dynamic.ang_vel += inv_inertia * r.cross(&impluse);
                dynamic.vel += dynamic.inv_mass * impluse;
                return;
            }
            let t = contact
                .normal
                .cross(&(vp.cross(&contact.normal)))
                .normalize();
            let vpt = vp.dot(&t);
            let cross_t = r.cross(&t);
            let m_eff_t = 1.0 / (dynamic.inv_mass + cross_t.dot(&(inv_inertia * cross_t)));
            let lambda_t = (lambda_n * 0.7).min(-vpt * m_eff_t);

            // apply impulse
            let impluse = lambda_n * contact.normal + lambda_t * t;
            dynamic.ang_vel += inv_inertia * r.cross(&impluse);
            dynamic.vel += dynamic.inv_mass * impluse;
        }

        if obj1.ty == RigidType::Dynamic {
            if obj2.ty == RigidType::Fixed {
                solve_dynamic_fixed(obj1, contact);
            }
        } else if obj1.ty == RigidType::Fixed {
            if obj2.ty == RigidType::Dynamic {
                // normal must point from the fixed to the dynamic for solve_dynamic_fixed()
                let mut reverted_contact = contact.clone();
                reverted_contact.normal = -reverted_contact.normal;
                solve_dynamic_fixed(obj2, contact);
            }
        }
    }

    pub fn gen_contacts(&mut self) {
        for i1 in 0..self.objects.len() {
            for i2 in 0..self.objects.len() {
                if i1 == i2 {
                    continue;
                }
                self.gen_contact_for(i1, i2);
            }
        }
    }

    pub fn gen_contact_for(&mut self, i1: usize, i2: usize) {
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

        fn sphere_contact(
            voxels1: &PVoxels,
            voxels2: &PVoxels,
            p1_in1: &Point3<f32>,
            p2_in2: &Point3<f32>,
            contacts: &mut Vec<Contact>,
            i1: usize,
            i2: usize,
        ) {
            let p1 = voxels1.transform * p1_in1;
            let p2 = voxels2.transform * p2_in2;
            let n = p1 - p2;
            let dist = n.magnitude();
            let overlap = (voxels1.dx + voxels2.dx) * 0.5 - dist;
            if overlap < 0.0 {
                return;
            }
            let n = 1.0 / dist * n;
            let point = p1.lerp(&p2, 0.5);
            contacts.push(Contact {
                point,
                normal: n,
                i1,
                i2,
            })
        }

        fn face_contact(
            voxels1: &PVoxels,
            voxels2: &PVoxels,
            p1_in1: &Point3<f32>,
            contacts: &mut Vec<Contact>,
            i1: usize,
            i2: usize,
            face_dir: FaceDir,
        ) {
            const DIRS: [Vector3<f32>; 6] = [
                Vector3::new(-1.0, 0.0, 0.0),
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(0.0, -1.0, 0.0),
                Vector3::new(0.0, 1.0, 0.0),
                Vector3::new(0.0, 0.0, -1.0),
                Vector3::new(0.0, 0.0, 1.0),
            ];
            let n = voxels2.transform.rotation * DIRS[face_dir as u8 as usize];
            let p1 = voxels1.transform * p1_in1;
            let point = p1 - voxels1.dx * 0.5 * n;
            contacts.push(Contact {
                point,
                normal: n,
                i1,
                i2,
            })
        }

        // 2. process
        let obj1 = &self.objects[i1];
        let obj2 = &self.objects[i2];

        let intersection_aabb = if let Some(aabb) = obj1.intersection_aabb(obj2) {
            aabb
        } else {
            return;
        };
        let intersection_part = (intersection_aabb + obj1.half_size) / obj1.dx;
        let start = component_floor(&intersection_part.min);
        let end = component_ceil(&intersection_part.max).inf(&obj1.shape.into());

        let to_obj2_transform = obj2.transform.inverse() * obj1.transform;
        let obj2_inv_dx = 1.0 / obj2.dx;

        for k in start.z..end.z {
            let z_base = k * obj1.area;
            let z = (k as f32 + 0.5) * obj1.dx - obj1.half_size.z;
            for j in start.y..end.y {
                let yz_base = z_base + j * obj1.shape.x;
                let y = (j as f32 + 0.5) * obj1.dx - obj1.half_size.y;
                for i in start.x..end.x {
                    let index = yz_base + i;

                    let ty1 = ty_of_data(obj1.data[index]);

                    if ty1 == CVoxelType::Edge {
                        let x = (i as f32 + 0.5) * obj1.dx - obj1.half_size.x;
                        let coord = Point3::new(x, y, z);
                        let coord_in_obj2 = to_obj2_transform * coord;
                        let unit_coord = (coord_in_obj2 + obj2.half_size) * obj2_inv_dx;
                        let [range_x, range_y, range_z] = axis_ranges(&unit_coord, &obj2.shape);

                        for rk in range_z {
                            let rz_base = rk * obj2.area;
                            let rz = (rk as f32 + 0.5) * obj2.dx - obj2.half_size.z;
                            for rj in range_y.clone() {
                                let rzy_base = rz_base + rj * obj2.shape.x;
                                let ry = (rj as f32 + 0.5) * obj2.dx - obj2.half_size.y;
                                for ri in range_x.clone() {
                                    let rindex = rzy_base + ri;

                                    let ty2 = ty_of_data(obj2.data[rindex]);

                                    if ty2 == CVoxelType::Edge {
                                        let rx = (ri as f32 + 0.5) * obj2.dx - obj2.half_size.x;
                                        let p2_in2 = Point3::new(rx, ry, rz);
                                        sphere_contact(
                                            obj1,
                                            obj2,
                                            &coord,
                                            &p2_in2,
                                            &mut self.contacts,
                                            i1,
                                            i2,
                                        );
                                    }
                                }
                            }
                        }
                    } else if ty1 == CVoxelType::Corner {
                        let x = (i as f32 + 0.5) * obj1.dx - obj1.half_size.x;
                        let p1_in1 = Point3::new(x, y, z);
                        let coord_in_obj2 = to_obj2_transform * p1_in1;
                        let unit_coord = (coord_in_obj2 + obj2.half_size) * obj2_inv_dx;
                        let [range_x, range_y, range_z] = axis_ranges(&unit_coord, &obj2.shape);

                        for rk in range_z {
                            let rz_base = rk * obj2.area;
                            let rz = (rk as f32 + 0.5) * obj2.dx - obj2.half_size.z;
                            for rj in range_y.clone() {
                                let rzy_base = rz_base + rj * obj2.shape.x;
                                let ry = (rj as f32 + 0.5) * obj2.dx - obj2.half_size.y;
                                for ri in range_x.clone() {
                                    let rindex = rzy_base + ri;

                                    let voxel_data2 = VoxelData::from_u8(obj2.data[rindex]);

                                    if voxel_data2.ty == CVoxelType::Face {
                                        face_contact(
                                            obj1,
                                            obj2,
                                            &p1_in1,
                                            &mut self.contacts,
                                            i1,
                                            i2,
                                            voxel_data2.dir,
                                        );
                                    } else if voxel_data2.ty == CVoxelType::Edge
                                        || voxel_data2.ty == CVoxelType::Corner
                                    {
                                        let rx = (ri as f32 + 0.5) * obj2.dx - obj2.half_size.x;
                                        let p2_in2 = Point3::new(rx, ry, rz);
                                        sphere_contact(
                                            obj1,
                                            obj2,
                                            &p1_in1,
                                            &p2_in2,
                                            &mut self.contacts,
                                            i1,
                                            i2,
                                        );
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
