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
    // Create a trimesh of the surface. You need to apply the transform of the CVoxel to show the mesh at the right position.
    pub fn surface_mesh(&self) -> Vec<[f32; 3]> {
        let mut vertices = Vec::new();

        let area = self.area();
        for k in 0..self.shape.z {
            let k_part = k * area;
            let kf32 = k as f32;
            for j in 0..self.shape.y {
                let j_part = j * self.shape.x;
                let jf32 = j as f32;
                for i in 0..self.shape.x {
                    let index = k_part + j_part + i;
                    let if32 = i as f32;
                    let coord = Vector3::new(if32, jf32, kf32);

                    if self.data[index] == CVoxelType::Air {
                        continue;
                    }

                    if k == 0 || self.data[index - area] == CVoxelType::Air {
                        for mut vertex in NZ_VERTICES {
                            vertex = self.dx * (vertex + coord);
                            vertices.push(vertex.data.0[0]);
                        }
                    }
                    if k == self.shape.z - 1 || self.data[index + area] == CVoxelType::Air {
                        for mut vertex in PZ_VERTICES {
                            vertex = self.dx * (vertex + coord);
                            vertices.push(vertex.data.0[0]);
                        }
                    }
                    if j == 0 || self.data[index - self.shape.x] == CVoxelType::Air {
                        for mut vertex in NY_VERTICES {
                            vertex = self.dx * (vertex + coord);
                            vertices.push(vertex.data.0[0]);
                        }
                    }
                    if j == self.shape.y - 1 || self.data[index + self.shape.x] == CVoxelType::Air {
                        for mut vertex in PY_VERTICES {
                            vertex = self.dx * (vertex + coord);
                            vertices.push(vertex.data.0[0]);
                        }
                    }
                    if i == 0 || self.data[index - 1] == CVoxelType::Air {
                        for mut vertex in NX_VERTICES {
                            vertex = self.dx * (vertex + coord);
                            vertices.push(vertex.data.0[0]);
                        }
                    }
                    if i == self.shape.x - 1 || self.data[index + 1] == CVoxelType::Air {
                        for mut vertex in PX_VERTICES {
                            vertex = self.dx * (vertex + coord);
                            vertices.push(vertex.data.0[0]);
                        }
                    }
                }
            }
        };

        vertices
    }
}
