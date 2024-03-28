// Computes error bound for the in_sphere predicate
// Inspired by CGAL source code
// Adapted to single-precision numbers.
// 
// Compile with g++ -I. main.cpp
// 
// Test computed bound for doubles and compare with CGAL:
//  CGAL/internal/Static_filters/Side_of_oriented_sphere_3.h:
//    1.2466136531027298e-13    OK
//    
//  Computed value for float:
//    6.6876506e-05

#include <CGALmini/Static_filter_error.h>

template <class T> inline T det2x2(
    T a11, T a12,
    T a21, T a22
) {
    return a11*a22 - a12*a21;
}

template <class T> inline T det3x3(
    T a11, T a12, T a13,
    T a21, T a22, T a23,
    T a31, T a32, T a33
) {
    return a11*det2x2(a22, a23, a32, a33) - a21*det2x2(a12, a13, a32, a33) + a31*det2x2(a12, a13, a22, a23);
}

template <class T> inline T det4x4(
    T a11, T a12, T a13, T a14,
    T a21, T a22, T a23, T a24,               
    T a31, T a32, T a33, T a34,  
    T a41, T a42, T a43, T a44  
) {
    T m12 = a21*a12 - a11*a22;
    T m13 = a31*a12 - a11*a32;
    T m14 = a41*a12 - a11*a42;
    T m23 = a31*a22 - a21*a32;
    T m24 = a41*a22 - a21*a42;
    T m34 = a41*a32 - a31*a42;
    
    T m123 = m23*a13 - m13*a23 + m12*a33;
    T m124 = m24*a13 - m14*a23 + m12*a43;
    T m134 = m34*a13 - m14*a33 + m13*a43;
    T m234 = m34*a23 - m24*a33 + m23*a43;
    
    return (m234*a14 - m134*a24 + m124*a34 - m123*a44);
}   

void compute_in_sphere_3d_bound_double() {
    typedef CGAL::Static_filter_error F;
    F t1 = F(1,F::ulp()/2);    // First translation
    F sq = t1*t1+t1*t1+t1*t1;  // squares
    F det = det4x4(t1, t1, t1, sq,
		   t1, t1, t1, sq,
		   t1, t1, t1, sq,
		   t1, t1, t1, sq); // Full det
    double err = det.error();
    err += err * 3 * F::ulp(); // Correction due to "eps * maxx * ...".
    std::cout.precision(std::numeric_limits<double>::max_digits10);
    std::cout << "Bound for double: " << err << std::endl;
}

void compute_in_sphere_3d_bound_float() {
    typedef CGAL::Static_filter_error_float F;    
    F t1 = F(1,F::ulp()/2);    // First translation
    F sq = t1*t1+t1*t1+t1*t1;  // squares
    F det = det4x4(t1, t1, t1, sq,
		   t1, t1, t1, sq,
		   t1, t1, t1, sq,
		   t1, t1, t1, sq); // Full det
    float err = det.error();
    err += err * 3 * F::ulp(); // Correction due to "eps * maxx * ...".
    std::cout.precision(std::numeric_limits<float>::max_digits10);    
    std::cout << "Bound for float: " << err << std::endl;
}

int main() {
    compute_in_sphere_3d_bound_double();
    compute_in_sphere_3d_bound_float();    
    return 0;
}
