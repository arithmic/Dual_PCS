use curve_traits::CurveParams;
use traits::traits::Field;

// Definition of the Line Struct
// line is in form c0 + c1*x + c2*y
#[derive(Clone, Copy, Debug)]
pub struct Line<C: CurveParams> {
    pub c0: C::FieldElement,
    pub c1: C::FieldElement,
    pub c2: C::FieldElement,
}

// Implemetation on the Line struct
impl<C: CurveParams> Line<C> {
    pub const fn new(c0: C::FieldElement, c1: C::FieldElement, c2: C::FieldElement) -> Self {
        return Self { c0, c1, c2 };
    }

    // Function to compute the tangent line
    pub fn tangent(p: C::ProjectivePoint) -> Self {
        let mut c0 = C::FieldElement::ZERO;
        let mut c1 = C::FieldElement::ZERO;
        let mut c2 = C::FieldElement::ZERO;
        if p.x.is_zero() && p.z.is_zero() {
            c2 = C::FieldElement::ONE;
        }
        //tangent line at point (x,y,z) has coefficients: c0=-3x^2, c1=2yz, c2=-(-y^2+12z^2)
        else {
            c0 = -C::FieldElement::from(3 as u128) * p.x.square();
            c1 = C::FieldElement::from(2 as u128) * (p.y * p.z);
            c2 = -(-p.y.square()
                + ((C::FieldElement::from(3 as u128) * <C as CurveParams>::EQUATION_B)
                    * p.z.square())); //
        }
        return Line::new(c0, c1, c2);
    }

    //function to compute line ax + by + c =0
    pub fn line_fn(p: C::ProjectivePoint, q: C::ProjectivePoint) -> Self {
        let mut c0 = C::FieldElement::ZERO;
        let mut c1 = C::FieldElement::ZERO;
        let mut c2 = C::FieldElement::ZERO;
        //if either of the point is point at infinity, line: ax+c=0
        if p.x.is_zero() && p.z.is_zero() {
            c0 = q.z;
            c2 = -q.x;
        } else if q.x.is_zero() && q.z.is_zero() {
            c0 = p.z;
            c2 = -p.x;
        } else if p == q {
            return Self::tangent(p);
        } else {
            if p.x == q.x {
                c0 = p.z;
                c2 = -p.x;
            }
            //line joining points (x1,y1,z1) and (x2,y2,z2) has coefficients: c0=(y1*z2-y2*z1), c1=-(x1*z2-x2*z1), c2=(x1*y2-x2*y1)
            else {
                c0 = (p.y * q.z) - (q.y * p.z);
                c1 = -((p.x * q.z) - (q.x * p.z));
                c2 = (p.x * q.y) - (p.y * q.x);
            }
        }
        return Line::new(c0, c1, c2);
    }
}


