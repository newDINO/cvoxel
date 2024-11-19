use nalgebra::Vector3;

use crate::{CVoxelType, CVoxels};

const NZ_VERTICES: [Vector3<f32>; 6] = [
    Vector3::new(0.0, 0.0, 0.0),
    Vector3::new(0.0, 1.0, 0.0),
    Vector3::new(1.0, 1.0, 0.0),
    Vector3::new(1.0, 1.0, 0.0),
    Vector3::new(1.0, 0.0, 0.0),
    Vector3::new(0.0, 0.0, 0.0),
];
const PZ_VERTICES: [Vector3<f32>; 6] = [
    Vector3::new(0.0, 0.0, 1.0),
    Vector3::new(1.0, 0.0, 1.0),
    Vector3::new(1.0, 1.0, 1.0),
    Vector3::new(1.0, 1.0, 1.0),
    Vector3::new(0.0, 1.0, 1.0),
    Vector3::new(0.0, 0.0, 1.0),
];
const NY_VERTICES: [Vector3<f32>; 6] = [
    Vector3::new(0.0, 0.0, 0.0),
    Vector3::new(1.0, 0.0, 0.0),
    Vector3::new(1.0, 0.0, 1.0),
    Vector3::new(1.0, 0.0, 1.0),
    Vector3::new(0.0, 0.0, 1.0),
    Vector3::new(0.0, 0.0, 0.0),
];
const PY_VERTICES: [Vector3<f32>; 6] = [
    Vector3::new(0.0, 1.0, 0.0),
    Vector3::new(0.0, 1.0, 1.0),
    Vector3::new(1.0, 1.0, 1.0),
    Vector3::new(1.0, 1.0, 1.0),
    Vector3::new(1.0, 1.0, 0.0),
    Vector3::new(0.0, 1.0, 0.0),
];
const NX_VERTICES: [Vector3<f32>; 6] = [
    Vector3::new(0.0, 0.0, 0.0),
    Vector3::new(0.0, 0.0, 1.0),
    Vector3::new(0.0, 1.0, 1.0),
    Vector3::new(0.0, 1.0, 1.0),
    Vector3::new(0.0, 1.0, 0.0),
    Vector3::new(0.0, 0.0, 0.0),
];
const PX_VERTICES: [Vector3<f32>; 6] = [
    Vector3::new(1.0, 0.0, 0.0),
    Vector3::new(1.0, 1.0, 0.0),
    Vector3::new(1.0, 1.0, 1.0),
    Vector3::new(1.0, 1.0, 1.0),
    Vector3::new(1.0, 0.0, 1.0),
    Vector3::new(1.0, 0.0, 0.0),
];

impl CVoxels {
    // Create a trimesh of the surface. You need to apply the transform of the CVoxels to show the mesh at the right position.
    pub fn surface_mesh(&self) -> Vec<[f32; 3]> {
        let mut vertices = Vec::new();
        let half_size = self.size() * 0.5;

        let area = self.area();
        for k in 0..self.shape.z {
            let k_part = k * area;
            let z_base = k as f32 * self.dx - half_size.z;
            for j in 0..self.shape.y {
                let j_part = j * self.shape.x;
                let y_base = j as f32 * self.dx - half_size.y;
                for i in 0..self.shape.x {
                    let index = k_part + j_part + i;
                    let x_base = i as f32 * self.dx - half_size.x;
                    let base = Vector3::new(x_base, y_base, z_base);

                    if self.data[index] == CVoxelType::Air {
                        continue;
                    }

                    if k == 0 || self.data[index - area] == CVoxelType::Air {
                        for mut vertex in NZ_VERTICES {
                            vertex = self.dx * vertex + base;
                            vertices.push(vertex.data.0[0]);
                        }
                    }
                    if k == self.shape.z - 1 || self.data[index + area] == CVoxelType::Air {
                        for mut vertex in PZ_VERTICES {
                            vertex = self.dx * vertex + base;
                            vertices.push(vertex.data.0[0]);
                        }
                    }
                    if j == 0 || self.data[index - self.shape.x] == CVoxelType::Air {
                        for mut vertex in NY_VERTICES {
                            vertex = self.dx * vertex + base;
                            vertices.push(vertex.data.0[0]);
                        }
                    }
                    if j == self.shape.y - 1 || self.data[index + self.shape.x] == CVoxelType::Air {
                        for mut vertex in PY_VERTICES {
                            vertex = self.dx * vertex + base;
                            vertices.push(vertex.data.0[0]);
                        }
                    }
                    if i == 0 || self.data[index - 1] == CVoxelType::Air {
                        for mut vertex in NX_VERTICES {
                            vertex = self.dx * vertex + base;
                            vertices.push(vertex.data.0[0]);
                        }
                    }
                    if i == self.shape.x - 1 || self.data[index + 1] == CVoxelType::Air {
                        for mut vertex in PX_VERTICES {
                            vertex = self.dx * vertex + base;
                            vertices.push(vertex.data.0[0]);
                        }
                    }
                }
            }
        };

        vertices
    }
}
