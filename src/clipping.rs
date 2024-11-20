use nalgebra::{Point2, UnitVector2, UnitVector3, Vector2};

pub struct Plane {
    pub normal: UnitVector3<f32>,
    /// Distance to the origin.
    pub distance: f32,
}

/// The origin point in the plane is the coords of the plane normal in 3d.
pub struct PlaneWithBasis {
    pub plane: Plane,
    pub x_basis: UnitVector3<f32>,
    pub y_basis: UnitVector3<f32>,
}

impl PlaneWithBasis {
    /// Returns None `if self.plane.normal.cross(&rhs.plane.normal).magnitude_squared() < epsilon`, 
    /// else return the intersection line in each representation.
    pub fn intersection_line(&self, rhs: &PlaneWithBasis, epsilon: f32) -> Option<(Line2d, Line2d)> {
        let v = self.plane.normal.cross(&rhs.plane.normal);
        if v.magnitude_squared() < epsilon {
            return None;
        }
        let n1 = v.cross(&self.plane.normal);
        let d1 = rhs.plane.distance * n1.dot(&rhs.plane.normal);
        let n2 = rhs.plane.normal.cross(&v);
        let d2 = self.plane.distance * n2.dot(&self.plane.normal);

        let n1_2d = Vector2::new(self.x_basis.dot(&n1), self.y_basis.dot(&n1));
        let n2_2d = Vector2::new(rhs.x_basis.dot(&n2), rhs.y_basis.dot(&n2));
        Some((
            Line2d {
                normal: UnitVector2::new_unchecked(n1_2d),
                distance: d1,
            },
            Line2d {
                normal: UnitVector2::new_unchecked(n2_2d),
                distance: d2,
            }
        ))
    }
}


pub struct Line2d {
    pub normal: UnitVector2<f32>,
    /// Distance to the origin.
    pub distance: f32,
}

impl Line2d {
    pub fn intersect_line_segment(&self, p1: Point2<f32>, p2: Point2<f32>) -> Option<Point2<f32>> {
        let d1 = p1.coords.dot(&self.normal);
        let d2 = p2.coords.dot(&self.normal);
        if (d1 - self.distance) * (d2 - self.distance) < 0.0 {
            Some(p1.lerp(&p2, (self.distance - d1) / (d2 - d1)))
        } else {
            None
        }
    }
}
