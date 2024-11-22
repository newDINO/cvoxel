use nalgebra::{Point2, UnitVector2, UnitVector3, Vector2, Vector3};

pub struct Line2dWithSide {
    pub normal: UnitVector2<f32>,
    /// Distance to the origin.
    pub distance: f32,
    pub out: bool,
}

impl Line2dWithSide {
    /// Currently unused by other functions.
    pub fn intersect_line_segment(&self, p1: Point2<f32>, p2: Point2<f32>) -> Option<Point2<f32>> {
        let d1 = p1.coords.dot(&self.normal) - self.distance;
        let d2 = p2.coords.dot(&self.normal) - self.distance;
        if d1 * d2 < 0.0 {
            Some(p1.lerp(&p2, -d1 / (d2 - d1)))
        } else {
            None
        }
    }

    pub fn intersect_line(&self, p1: Point2<f32>, p2: Point2<f32>) -> Point2<f32> {
        let d1 = p1.coords.dot(&self.normal) - self.distance;
        let d2 = p2.coords.dot(&self.normal) - self.distance;
        p1.lerp(&p2, -d1 / (d2 - d1))
    }

    pub fn intersect_ray(&self, origin: &Point2<f32>, point_at: &Point2<f32>) -> Option<Point2<f32>> {
        let d_start = origin.coords.dot(&self.normal) - self.distance;
        let d_dir = point_at.coords.dot(&self.normal) - self.distance;
        let weight = -d_start / (d_dir - d_start);
        if weight > 0.0 {
            Some(origin.lerp(&point_at, weight))
        } else {
            None
        }
    }

    /// If `out` is true, then `point.coord.dot(&line.normal) > line.distance` will be clipped,
    /// else `point.coord.dot(&line.normal) < line.distance` will be clipped.
    ///
    /// Panic if the polygon contains no vertices.
    pub fn clip_polygon(&self, polygon: &[Point2<f32>]) -> Vec<Point2<f32>> {
        let mut result = Vec::new();

        let mut vp = polygon[0];
        let mut dp = vp.coords.dot(&self.normal) - self.distance;
        let mut vp_reserve = self.out && dp < 0.0 || !self.out && dp > 0.0;
        for i in 1..polygon.len() {
            let vn = polygon[i];
            let dn = vn.coords.dot(&self.normal) - self.distance;
            let vn_reserve = self.out && dn < 0.0 || !self.out && dn > 0.0;

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
}

#[derive(Clone, Copy)]
pub struct Plane {
    pub normal: UnitVector3<f32>,
    /// Distance to the origin.
    pub distance: f32,
}

/// The origin point in the plane is the coords of the plane normal in 3d.
#[derive(Clone, Copy)]
pub struct PlaneWithBasis {
    pub plane: Plane,
    pub x_basis: UnitVector3<f32>,
    pub y_basis: UnitVector3<f32>,
    pub out: bool,
}

impl PlaneWithBasis {
    /// Create one configuration.
    pub fn from_plane(plane: Plane, out: bool) -> Self {
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
            plane, x_basis, y_basis, out
        }
    }

    /// Returns None `if self.plane.normal.cross(&rhs.plane.normal).magnitude_squared() < epsilon`,
    /// else return the intersection line in each representation.
    pub fn intersection_lines(&self, rhs: &PlaneWithBasis, epsilon: f32) -> Option<(Line2dWithSide, Line2dWithSide)> {
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
            Line2dWithSide {
                normal: UnitVector2::new_unchecked(n1_2d),
                distance: d1,
                out: rhs.out,
            },
            Line2dWithSide {
                normal: UnitVector2::new_unchecked(n2_2d),
                distance: d2,
                out: self.out,
            }
        ))
    }

    /// Project a 3d vector to the plane.
    pub fn vertor_in_plane(&self, vector: &Vector3<f32>) -> Vector2<f32> {
        Vector2::new(self.x_basis.dot(&vector), self.y_basis.dot(&vector))
    }
}

pub struct PolygonPlane {
    pub plane: PlaneWithBasis,
    pub vertices: Vec<Point2<f32>>,
}

impl PolygonPlane {
    pub fn clip_with_polygon(&mut self, rhs: &mut PolygonPlane, epsilon: f32) {
        let lines = self.plane.intersection_lines(&rhs.plane, epsilon);
        let (line1, line2) = if let Some(lines) = lines {
            lines
        } else {
            return;
        };
        self.vertices = line1.clip_polygon(&self.vertices);
        rhs.vertices = line2.clip_polygon(&rhs.vertices);
    }
    pub fn clip_with_clippable(&mut self, rhs: &mut Clippable, epsilon: f32) {
        match rhs {
            Clippable::Empty(plane) => {
                let lines = self.plane.intersection_lines(plane, epsilon);
                let (line1, line2) = if let Some(lines) = lines {
                    lines
                } else {
                    return;
                };
                self.vertices = line1.clip_polygon(&self.vertices);
                *rhs = Clippable::Line(LineInPlane {
                    plane: *plane, line: line2,
                })
            },
            Clippable::Line(line_in_plane) => {
                let lines = self.plane.intersection_lines(&line_in_plane.plane, epsilon);
                let (line1, line2) = if let Some(lines) = lines {
                    lines
                } else {
                    return;
                };
                self.vertices = line1.clip_polygon(&self.vertices);

                let dir2 = Vector2::new(-line2.normal.y, line2.normal.x);
                if dir2.dot(&line_in_plane.line.normal) < epsilon {
                    return;
                }
                let p_start = Point2::from(line_in_plane.line.normal.into_inner() * line_in_plane.line.distance);
                let p_end = Point2::from(line2.normal.into_inner() * line2.distance);
                let p_temp = p_end + dir2;
                let cross = line_in_plane.line.intersect_line(p_end, p_temp);
                *rhs = Clippable::Open(OpenLinesInPlane {
                    plane: line_in_plane.plane,
                    points: vec![p_start, cross, p_end],
                })
            },
            Clippable::Open(open_lines) => {
                let lines = self.plane .intersection_lines(&open_lines.plane, epsilon);
                let (line1, line2) = if let Some(lines) = lines {
                    lines
                } else {
                    return;
                };
                self.vertices = line1.clip_polygon(&self.vertices);

                let start_origin = open_lines.points[1];
                let d_start_origin = start_origin.coords.dot(&line2.normal) - line2.distance;
                let start_origin_in = d_start_origin > 0.0 && line2.out || d_start_origin < 0.0 && !line2.out;

                let end_origin = open_lines.points[open_lines.points.len() - 2];
                let d_end_origin = end_origin.coords.dot(&line2.normal) - line2.distance;
                let end_origin_in = d_end_origin > 0.0 && line2.out || d_end_origin < 0.0 && !line2.out;

                let start_cross = line2.intersect_ray(&start_origin, &open_lines.points[1]);
                let end_cross = line2.intersect_ray(&end_origin, open_lines.points.last().unwrap());

                if start_origin_in {
                    if let Some(start_cross) = start_cross {
                        if end_origin_in {
                            if let Some(end_cross) = end_cross {

                            }
                        }
                    }
                }
            }
            _ => {}
        }
    }
}

pub struct LineInPlane {
    pub plane: PlaneWithBasis,
    pub line: Line2dWithSide,
}

pub struct OpenLinesInPlane {
    pub plane: PlaneWithBasis,
    /// The first point and the last point are not vertices, they indicate the direction of the ray.
    pub points: Vec<Point2<f32>>,
}

pub enum Clippable {
    Empty(PlaneWithBasis),
    Line(LineInPlane),
    Open(OpenLinesInPlane),
    Polygon(PolygonPlane),
}
