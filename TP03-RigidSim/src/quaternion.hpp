

#ifndef _QUATERNION_H_
#define _QUATERNION_H_

#include "typedefs.hpp"
#include "Vector3.hpp"
#include "Matrix3x3.hpp"
#include "math.h"

class Quaternion {
public:
  tReal qs, qx, qy, qz;

//public:
  //Quaternion() : qs(1.0), qx(0.0), qy(0.0), qz(0.0) {}

  Quaternion(tReal s, tReal x, tReal y, tReal z) : qs(s), qx(x), qy(y), qz(z) {}

  Quaternion(const Vec3f vec, const tReal s) {
    qx = vec.x;
    qy = vec.y;
    qz = vec.z;
    qs = s;
  }

  Quaternion operator/= (Quaternion &div) const {
    return Quaternion(*this) /= div;
  }

  Quaternion operator* (const Quaternion &prod) {
    qs *= prod.qs;
    qx *= prod.qx;
    qy *= prod.qy;
    qz *= prod.qz;

    return *this;
  }

  Quaternion operator* (const tReal &prod) {
    qs *= prod;
    qx *= prod;
    qy *= prod;
    qz *= prod;

    return *this;
  }

  Quaternion operator/ (const tReal &div) {
    qs /= div;
    qx /= div;
    qy /= div;
    qz /= div;

    return *this;
  }

  Quaternion operator+ (const tReal &sum) {
    qs += sum;
    qx += sum;
    qy += sum;
    qz += sum;

    return *this;
  }

  Quaternion operator+ (const Quaternion &sum) {
    qs += sum.qs;
    qx += sum.qx;
    qy += sum.qy;
    qz += sum.qz;

    return *this;
  }

  Mat3f rotationMatrix() const {
    return Mat3f(
        1 - 2*qy*qy - 2*qz*qz,  2*qx*qy - 2*qs*qz,      2*qx*qz + 2*qs*qy,
        2*qx*qy + 2*qs*qz,      1 - 2*qx*qx - 2*qz*qz,  2*qy*qz - 2*qs*qx,
        2*qx*qz - 2*qs*qy,      2*qy*qz + 2*qs*qx,      1 - 2*qx*qx - 2*qy*qy
    );
  }

  Quaternion product(Quaternion &qua) const {
    return Quaternion(
      qs*qua.qs - qx*qua.qx - qy*qua.qy - qz*qua.qz,
      qs*qua.qx + qx*qua.qs + qy*qua.qz - qz*qua.qy,
      qs*qua.qy - qx*qua.qz + qy*qua.qs + qz*qua.qx,
      qs*qua.qz + qx*qua.qy - qy*qua.qx + qz*qua.qs
    );
  }

  tReal normSqrt() {
    return sqrt(qx*qx + qy*qy + qz*qz + qs*qs);
  }

  Quaternion normalize() {
    return Quaternion(*this) / normSqrt();
  }

};

#endif  /* _QUATERNION_H_ */
