use nalgebra::{Point2, UnitVector2, UnitVector3, Vector2, Vector3};

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
    /// Create one configuration.
    pub fn from_plane(plane: Plane) -> Self {
        let up1 = Vector3::y_axis();
        let up2 = Vector3::z_axis();
        
        let up = if plane.normal.dot(&up1) > plane.normal.dot(&up2) {
            up1
        } else {
            up2
        };
        
        let x_basis = UnitVector3::new_unchecked(plane.normal.cross(&up));
        let y_basis = UnitVector3::new_unchecked(plane.normal.cross(&x_basis));
        Self {
            plane, x_basis, y_basis
        }
    }

    /// Returns None `if self.plane.normal.cross(&rhs.plane.normal).magnitude_squared() < epsilon`, 
    /// else return the intersection line in each representation.
    pub fn intersection_lines(&self, rhs: &PlaneWithBasis, epsilon: f32) -> Option<(Line2d, Line2d)> {
        let v = self.plane.normal.cross(&rhs.plane.normal);
        if v.magnitude_squared() < epsilon {
            return None;
        }
        let n1 = v.cross(&self.plane.normal);
        let n2 = rhs.plane.normal.cross(&v);

        let cos = self.plane.normal.dot(&rhs.plane.normal);
        let sin = (1.0 - cos * cos).sqrt();
        let sin_inv = 1.0 / sin;
        let cot = cos * sin_inv;
        let d1 = rhs.plane.distance * sin_inv - self.plane.distance * cot;
        let d2 = self.plane.distance * sin_inv - rhs.plane.distance * cot;

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

    /// Project a 3d vector to the plane.
    pub fn vertor_in_plane(&self, vector: &Vector3<f32>) -> Vector2<f32> {
        Vector2::new(self.x_basis.dot(&vector), self.y_basis.dot(&vector))
    }
}


pub struct Line2d {
    pub normal: UnitVector2<f32>,
    /// Distance to the origin.
    pub distance: f32,
}

impl Line2d {
    /// Currently unused by other functions.
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

pub struct PolygonPlane {
    pub plane: PlaneWithBasis,
    pub vertices: Vec<Point2<f32>>,
}

impl PolygonPlane {
    pub fn clip_each_other(&mut self, rhs: &mut PolygonPlane, out1: bool, out2: bool, epsilon: f32) {
        let lines = self.plane.intersection_lines(&rhs.plane, epsilon);
        let (line1, line2) = if let Some(lines) = lines {
            lines
        } else {
            return;
        };

        let np1 = if out1 {
            self.plane.plane.normal
        } else {
            -self.plane.plane.normal
        };
        let np2 = if out2 {
            rhs.plane.plane.normal
        } else {
            -rhs.plane.plane.normal
        }; 

        let np2_in_1 = self.plane.vertor_in_plane(&np2);
        let np1_in_2 = rhs.plane.vertor_in_plane(&np1);

        let line1_out = np2_in_1.dot(&line1.normal) > 0.0;
        let line2_out = np1_in_2.dot(&line2.normal) > 0.0;

        self.vertices = line_clip_polygon(&line1, &self.vertices, line1_out);
        rhs.vertices = line_clip_polygon(&line2, &rhs.vertices, line2_out);
    }
}

/// If `out` is true, then `point.coord.dot(&line.normal) > line.distance` will be clipped,
/// else `point.coord.dot(&line.normal) < line.distance` will be clipped.
/// 
/// Panic if the polygon contains no vertices.
pub fn line_clip_polygon(line: &Line2d, polygon: &Vec<Point2<f32>>, out: bool) -> Vec<Point2<f32>> {
    let mut result = Vec::new();

    let mut vp = polygon[0];
    let mut dp = vp.coords.dot(&line.normal) - line.distance;
    let mut vp_reserve = out && dp < 0.0 || !out && dp > 0.0;
    for i in 1..polygon.len() {
        let vn = polygon[i];
        let dn = vn.coords.dot(&line.normal) - line.distance;
        let vn_reserve = out && dn < 0.0 || !out && dn > 0.0;

        if vp_reserve {
            result.push(vp);
            if !vn_reserve {
                result.push(vp.lerp(&vn, -dp / (dn - dp)));
            }
        } else {
            if vn_reserve {
                result.push(vp.lerp(&vn, -dp / (dn - dp)));
            }
        }

        vp = vn;
        dp = dn;
        vp_reserve = vn_reserve;
    }
    result
}