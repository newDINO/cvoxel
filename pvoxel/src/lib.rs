use cvoxel::{ty_of_data, CVoxelType, FaceDir, VoxelData};
use math::LMatrix;
use nalgebra::{Point3, SVector, Vector3};

pub mod math;
mod object;
pub use object::*;

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
        for _iter in 0..iters {
            for i in 0..self.contacts.len() {
                // let restitution = 1.0 - iter as f32 * 0.5 / (iters as f32 - 1.0);
                self.resolve_one_contact(i, 1.0);
            }
        }
        self.contacts.clear();
    }

    pub fn resolve_group_normal(&mut self, start: usize) {
        // 1. consts
        const DIM: usize = 4;
        const N: usize = DIM * (DIM + 1) / 2;
        const RESTITUTION: f32 = 0.7;

        // 2. utils
        fn dynamic_fixed_params(
            i: usize,
            end: usize,
            b_data: &mut [f32; DIM],
            lmatrix: &mut LMatrix<N>,
            objects: &[PVoxels],
            contacts: &[Contact],
            dynamic_i: usize,
            contact: &Contact
        ) {
            let dynamic = &objects[dynamic_i];
            let cm = dynamic.transform * dynamic.local_cm;
            let r = contact.point - cm;

            let vp = dynamic.vel + dynamic.ang_vel.cross(&r);
            let vpn = vp.dot(&contact.normal);
            
            if vpn > 0.0 {
                b_data[i] = 0.0;
            } else {
                b_data[i] = -2.0 * RESTITUTION * vpn;
            }

            let rot = dynamic.transform.rotation.to_rotation_matrix();
            let inv_inertia = rot * dynamic.inv_local_inertia * rot.transpose();

            let cross = r.cross(&contact.normal);

            for j in i..end {
                if dynamic_i == contacts[j].i1 {
                    let object_j = &objects[contacts[j].i1];
                    let nj = contacts[j].normal;
                    let rj = contact.point - object_j.transform * object_j.local_cm;
                    lmatrix[(j, i)] = dynamic.inv_mass * contact.normal.dot(&nj) + cross.dot(&(inv_inertia * rj.cross(&nj)));
                } else if dynamic_i == contacts[j].i2 {
                    let object_j = &objects[contacts[j].i1];
                    let nj = -contacts[j].normal;
                    let rj = contact.point - object_j.transform * object_j.local_cm;
                    lmatrix[(j, i)] = dynamic.inv_mass * contact.normal.dot(&nj) + cross.dot(&(inv_inertia * rj.cross(&nj)));
                }
            }
        }
        fn apply_impulse(object: &mut PVoxels, contact: &Contact, impulse: &Vector3<f32>) {
            let cm = object.transform * object.local_cm;
            let r = contact.point - cm;
            let rot = object.transform.rotation.to_rotation_matrix();
            let inv_inertia = rot * object.inv_local_inertia * rot.transpose();
            object.vel += object.inv_mass * impulse;
            object.ang_vel += inv_inertia * r.cross(&impulse);
        }

        // 3. process
        let mut b_data = [0.0; DIM];
        let mut lmatrix = LMatrix { data: [0.0; N] };

        let end = start + DIM;

        for i in start..end {
            let mut contact = self.contacts[i];
            let voxeli1 = &self.objects[contact.i1];
            let voxeli2 = &self.objects[contact.i2];
            if voxeli1.ty == RigidType::Dynamic {
                if voxeli2.ty == RigidType::Dynamic {

                } else if voxeli2.ty == RigidType::Fixed {
                    dynamic_fixed_params(i, end, &mut b_data, &mut lmatrix, &self.objects, &self.contacts, contact.i1, &contact);
                }
            } else if voxeli1.ty == RigidType::Fixed {
                if voxeli2.ty == RigidType::Dynamic {
                    contact.normal = -contact.normal;
                    dynamic_fixed_params(i, end, &mut b_data, &mut lmatrix, &self.objects, &self.contacts, contact.i2, &contact);
                }
            }
        }

        let mut lambdas = SVector::from(b_data);
        lmatrix.solve_mut(&mut lambdas);
        for i in 0..DIM {
            let contact = &self.contacts[i + start];
            let impulse = contact.normal * lambdas[i];
            
            let obj1 = &mut self.objects[contact.i1];
            if obj1.ty == RigidType::Dynamic {
                apply_impulse(obj1, contact, &impulse);
            }
            let obj2 = &mut self.objects[contact.i2];
            if obj2.ty == RigidType::Dynamic {
                apply_impulse(obj2, contact, &-impulse);
            }
        }
    }

    /// This will apply the impulse after it is solved.
    pub fn resolve_one_contact(&mut self, index: usize, restitution: f32) {
        // 1. utils
        fn solve_dynamic_dynamic(d1: &mut PVoxels, d2: &mut PVoxels, contact: &Contact, restitution: f32) {
            let cm1 = d1.transform * d1.local_cm;
            let r1 = contact.point - cm1;
            let cm2 = d2.transform * d2.local_cm;
            let r2 = contact.point - cm2;

            let n1 = contact.normal;
            let n2 = -n1;

            let vp1 = d1.vel + d1.ang_vel.cross(&r1);
            let vpn1 = vp1.dot(&n1);
            let vp2 = d2.vel + d2.ang_vel.cross(&r2);
            let vpn2 = vp2.dot(&n2);

            if vpn1 + vpn2 > 0.0 {
                return;
            }
            let dv = -2.0 * restitution * (vpn1 + vpn2);

            let rot1 = d1.transform.rotation.to_rotation_matrix();
            let inv_inertia1 = rot1 * d1.inv_local_inertia * rot1.transpose();
            let rot2 = d2.transform.rotation.to_rotation_matrix();
            let inv_inertia2 = rot2 * d2.inv_local_inertia * rot2.transpose();

            let cross_n1 = r1.cross(&n1);
            let cross_n2 = r2.cross(&n2);

            let kn1 = d1.inv_mass + cross_n1.dot(&(inv_inertia1 * cross_n1));
            let kn2 = d2.inv_mass + cross_n2.dot(&(inv_inertia2 * cross_n2));
            let m_eff_n = 1.0 / (kn1 + kn2);
            let lambda_n = dv * m_eff_n;

            let impulse_n1 = lambda_n * n1;
            let impulse_n2 = lambda_n * n2;
            d1.vel += d1.inv_mass * impulse_n1;
            d1.ang_vel += inv_inertia1 * r1.cross(&impulse_n1);
            d2.vel += d2.inv_mass * impulse_n2;
            d2.ang_vel += inv_inertia2 * r2.cross(&impulse_n2);
        }

        fn solve_dynamic_fixed(dynamic: &mut PVoxels, contact: &Contact, restitution: f32) {
            let cm = dynamic.transform * dynamic.local_cm;
            let r = contact.point - cm;

            let vp = dynamic.vel + dynamic.ang_vel.cross(&r);
            let vpn = vp.dot(&contact.normal);
            if vpn > 0.0 {
                return;
            }
            let dv = -2.0 * restitution * vpn;

            let rot = dynamic.transform.rotation.to_rotation_matrix();
            let inv_inertia = rot * dynamic.inv_local_inertia * rot.transpose();

            let cross_n = r.cross(&contact.normal);

            let m_eff_n = 1.0 / (dynamic.inv_mass + cross_n.dot(&(inv_inertia * cross_n)));
            let lambda_n = dv * m_eff_n;

            // friction
            let t = vp.cross(&contact.normal).cross(&contact.normal);

            let t_mag = t.magnitude();
            if t_mag < 0.05 {
                let impluse = lambda_n * contact.normal;
                println!("apply impulse: {:?}", impluse);
                dynamic.ang_vel += inv_inertia * r.cross(&impluse);
                dynamic.vel += dynamic.inv_mass * impluse;
                return;
            }
            let t = 1.0 / t_mag * t;

            let vpt = vp.dot(&t);
            let cross_t = r.cross(&t);
            let m_eff_t = 1.0 / (dynamic.inv_mass + cross_t.dot(&(inv_inertia * cross_t)));
            let lambda_t = (lambda_n * 0.7).min(-vpt * m_eff_t);

            // apply impulse
            let impluse = lambda_n * contact.normal + lambda_t * t;
            println!("apply impulse: {:?}", impluse);
            dynamic.ang_vel += inv_inertia * r.cross(&impluse);
            dynamic.vel += dynamic.inv_mass * impluse;
        }

        // 2. process
        let contact = &self.contacts[index];
        let (obj1, obj2) = if contact.i1 < contact.i2 {
            let (part1, part2) = self.objects.split_at_mut(contact.i2);
            (&mut part1[contact.i1], &mut part2[0])
        } else {
            let (part1, part2) = self.objects.split_at_mut(contact.i1);
            (&mut part2[0], &mut part1[contact.i2])
        };


        if obj1.ty == RigidType::Dynamic {
            if obj2.ty == RigidType::Fixed {
                solve_dynamic_fixed(obj1, contact, restitution);
            } else if obj2.ty == RigidType::Dynamic {
                solve_dynamic_dynamic(obj1, obj2, contact, restitution);
            }
        } else if obj1.ty == RigidType::Fixed {
            if obj2.ty == RigidType::Dynamic {
                // normal must point from the fixed to the dynamic for solve_dynamic_fixed()
                let mut reverted_contact = contact.clone();
                reverted_contact.normal = -reverted_contact.normal;
                solve_dynamic_fixed(obj2, contact, restitution);
            }
        }
    }

    pub fn gen_contacts(&mut self) {
        for i1 in 0..self.objects.len() {
            for i2 in 0..self.objects.len() {
                if i1 == i2
                    || (self.objects[i1].ty == RigidType::Dynamic
                        && self.objects[i2].ty == RigidType::Dynamic)
                {
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
