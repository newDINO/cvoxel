// use crate::{CVoxelType, CVoxels};

// pub struct TriMesh {
//     pub vertices: Vec<[f32; 3]>,
// }

// impl CVoxels {
//     // Create a trimesh of the surface. You need to apply the transform of the CVoxel to show the mesh at the right position.
//     pub fn surface_mesh(&self) -> TriMesh {
//         let area = self.area();
//         for k in 0..self.shape.z {
//             let k_part = k * area;
//             for j in 0..self.shape.y {
//                 let j_part = j * self.shape.x;
//                 for i in 0..self.shape.x {
//                     let coord = k_part + j_part + i;

//                     if self.data[coord] == CVoxelType::Air {
//                         continue;
//                     }

//                     if k == 0 || self.data[coord - area] == CVoxelType::Air {

//                     }

//                 }
//             }
//         };
//     }
// }
