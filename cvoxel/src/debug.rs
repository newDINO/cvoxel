use nalgebra::Vector3;

use crate::{ty_of_data, CVoxelType, CVoxels};

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

pub struct SufMesh {
    pub position: Vec<[f32; 3]>,
    pub color: Vec<[f32; 4]>,
}

impl CVoxels {
    // Create a trimesh of the surface. You need to apply the transform of the CVoxels to show the mesh at the right position.
    pub fn surface_mesh(&self) -> SufMesh {
        let mut vertices = Vec::new();
        let mut colors = Vec::new();

        for k in 0..self.shape.z {
            let k_part = k * self.area;
            let z_base = k as f32 * self.dx - self.half_size.z;
            for j in 0..self.shape.y {
                let j_part = j * self.shape.x;
                let y_base = j as f32 * self.dx - self.half_size.y;
                for i in 0..self.shape.x {
                    let index = k_part + j_part + i;
                    let x_base = i as f32 * self.dx - self.half_size.x;
                    let base = Vector3::new(x_base, y_base, z_base);

                    let color = match ty_of_data(self.data[index]) {
                        CVoxelType::Corner => [0.0, 1.0, 1.0, 1.0],
                        CVoxelType::Edge => [1.0, 1.0, 0.0, 1.0],
                        CVoxelType::Face => [1.0, 1.0, 1.0, 1.0],
                        _ => continue,
                    };

                    let push_vertex = &mut |ty: &[Vector3<f32>; 6]| {
                        for vertex in ty {
                            let vertex = self.dx * vertex + base;
                            vertices.push(vertex.data.0[0]);
                            colors.push(color);
                        }
                    };

                    if k == 0 || ty_of_data(self.data[index - self.area]) == CVoxelType::Air {
                        push_vertex(&NZ_VERTICES);
                    }
                    if k == self.shape.z - 1
                        || ty_of_data(self.data[index + self.area]) == CVoxelType::Air
                    {
                        push_vertex(&PZ_VERTICES);
                    }
                    if j == 0 || ty_of_data(self.data[index - self.shape.x]) == CVoxelType::Air {
                        push_vertex(&NY_VERTICES);
                    }
                    if j == self.shape.y - 1
                        || ty_of_data(self.data[index + self.shape.x]) == CVoxelType::Air
                    {
                        push_vertex(&PY_VERTICES);
                    }
                    if i == 0 || ty_of_data(self.data[index - 1]) == CVoxelType::Air {
                        push_vertex(&NX_VERTICES);
                    }
                    if i == self.shape.x - 1 || ty_of_data(self.data[index + 1]) == CVoxelType::Air
                    {
                        push_vertex(&PX_VERTICES);
                    }
                }
            }
        }
        SufMesh {
            position: vertices,
            color: colors,
        }
    }
}
