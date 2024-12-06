use nalgebra::SVector;
use std::ops::{Index, IndexMut};

#[derive(Debug, Clone, Copy)]
pub struct LMatrix<const N: usize> {
    pub data: [f32; N],
}
impl<const N: usize> Index<(usize, usize)> for LMatrix<N> {
    type Output = f32;
    #[inline]
    fn index(&self, index: (usize, usize)) -> &Self::Output {
        &self.data[(index.0 * (index.0 + 1) >> 1) + index.1]
    }
}

impl<const N: usize> IndexMut<(usize, usize)> for LMatrix<N> {
    #[inline]
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        &mut self.data[(index.0 * (index.0 + 1) >> 1) + index.1]
    }
}

impl<const N: usize> LMatrix<N> {
    pub fn solve_mut<const DIM: usize>(&mut self, b: &mut SVector<f32, DIM>) {
        let mut ts = [0.0; DIM];
        let mut inv_ds = [0.0; DIM];
        inv_ds[0] = 1.0 / self.data[0];
        for i in 1..DIM {
            // aij -> tij
            for j in 0..i {
                let mut tij_part = 0.0;
                for k in 0..j {
                    tij_part += ts[k] * self[(j, k)];
                }
                ts[j] = self[(i, j)] - tij_part;
            }
            // tij -> lij
            for j in 0..i {
                self[(i, j)] = ts[j] * inv_ds[j];
            }
            let mut di_part = 0.0;
            // di
            for k in 0..i {
                di_part += ts[k] * self[(i, k)];
            }
            let di = self[(i, i)] - di_part;
            inv_ds[i] = 1.0 / di;
        }

        for i in 1..DIM {
            let mut part_yi = 0.0;
            for k in 0..i {
                part_yi += self[(i, k)] * b[k];
            }
            b[i] = b[i] - part_yi;
        }
        let mut i = DIM - 1;
        loop {
            let mut part_xi = 0.0;
            for k in i + 1..DIM {
                part_xi += self[(k, i)] * b[k];
            }
            b[i] = b[i] * inv_ds[i] - part_xi;
            if i == 0 {
                break;
            }
            i -= 1;
        }
    }
}
