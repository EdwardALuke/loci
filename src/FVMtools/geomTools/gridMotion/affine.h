//#############################################################################
//#
//# Copyright 2014-2025, Mississippi State University
//#
//# The GridMover is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The GridMover software is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the GridMover software. If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#ifndef AFFINE_H
#define AFFINE_H


#include "gridTypes.h"

namespace gridMotion {

  class AffineMatrix {
    /** Helper class - Matrix Row */
    typedef Array<double,4> Row;

  public:
    Array<Row,4> M;

  public:
    /** Default constructor - Identity transformation */
    AffineMatrix() {
      for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
          M[i][j] = ( i==j ? 1.0 : 0.0 );
        }
      }
    }
    /** Conversion constructor - Translation transformation */
    AffineMatrix(const vect3d & v) {
      for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
          M[i][j] = ( i==j ? 1.0 : 0.0 );
        }
      }
      M[0][3] = v.x;
      M[1][3] = v.y;
      M[2][3] = v.z;
    }
    /** Conversion Constructor - Rotation transformation */
    AffineMatrix(const tens3d & R) {
      for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
          M[i][j] = ( i==j ? 1.0 : 0.0 );
        }
      }
      M[0][0] = R.x.x;
      M[0][1] = R.x.y;
      M[0][2] = R.x.z;
      M[1][0] = R.y.x;
      M[1][1] = R.y.y;
      M[1][2] = R.y.z;
      M[2][0] = R.z.x;
      M[2][1] = R.z.y;
      M[2][2] = R.z.z;
    }

  public:
    /** Array index accessor - non-const */
    Row & operator[](const int i) {
      return M[i];
    }
    /** Array index accessor - const */
    const Row & operator[](const int i) const {
      return M[i];
    }

  public:
    /** Addition with assign operator */
    AffineMatrix & operator+=( const AffineMatrix & N ) {
      for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
          M[i][j] += N[i][j];
        }
      }
      return (*this);
    }
    /** Multiplication with assign operator */
    AffineMatrix & operator*=( const AffineMatrix & N ) {
      double A[4][4];
      for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
          A[i][j] = 0.0;
          for (int k=0; k<4; ++k) {
            A[i][j] += M[i][k]*N[k][j];
          }
        }
      }
      for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
          M[i][j] = A[i][j];
        }
      }
      return (*this);
    }
    /** Multiplication by a scalar with assign */
    AffineMatrix & operator*=( const double s ) {
      for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
          M[i][j] *= s;
        }
      }
      return (*this);
    }
    /** Division by a scalar with assign */
    AffineMatrix & operator/=( const double s ) {
      for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
          M[i][j] /= s;
        }
      }
      return (*this);
    }

  public:
    /** Translation with assign */
    AffineMatrix & translate(const vect3d v) {
      M[0][3] += v.x;
      M[1][3] += v.y;
      M[2][3] += v.z;
      return (*this);
    }
    /** Rotation with assign
     *
     * NOTE: To properly apply rotations, the rotation matrix must be
     * left multiplied into the current matrix.
     */
    AffineMatrix & rotate(const tens3d R) {
      AffineMatrix N(R);
      double A[4][4];
      for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
          A[i][j] = 0.0;
          for (int k=0; k<4; ++k) {
            A[i][j] += N[i][k]*M[k][j];
          }
        }
      }
      for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
          M[i][j] = A[i][j];
        }
      }
      return (*this);
    }

  public:
    const tens3d rotation() const {
      return tens3d(vect3d(M[0][0], M[0][1], M[0][2]),
                    vect3d(M[1][0], M[1][1], M[1][2]),
                    vect3d(M[2][0], M[2][1], M[2][2]));
    }
    const vect3d translation() const {
      return vect3d(M[0][3], M[1][3], M[2][3]);
    }
  };

  /** Matrix-Matrix Sum */
  inline AffineMatrix operator+( const AffineMatrix & M,
                          const AffineMatrix & N ) {
    AffineMatrix A(M);
    A += N;
    return A;
  }
  /** Matrix-Matrix Product */
  inline AffineMatrix operator*( const AffineMatrix & M,
                          const AffineMatrix & N ) {
    AffineMatrix A(M);
    A *= N;
    return A;
  }
  /** Scalar-Matrix Product */
  inline AffineMatrix operator*( const double s,
                          const AffineMatrix & M ) {
    AffineMatrix A(M);
    A *= s;
    return A;
  }
  /** Matrix-Scalar Product */
  inline AffineMatrix operator*( const AffineMatrix & M,
                          const double s ) {
    AffineMatrix A(M);
    A *= s;
    return A;
  }
  /** Matrix-Scalar Division */
  inline AffineMatrix operator/( const AffineMatrix & M,
                          const double s ) {
    AffineMatrix A(M);
    A /= s;
    return A;
  }
  /** Matrix-Vector Product
   *
   *  NOTE: We are implicitly promoting the vect3d to a homogeneous
   *  vector.
   */
  inline vect3d operator*( const AffineMatrix & A, const vect3d & v ) {
    return vect3d(A[0][0]*v.x + A[0][1]*v.y + A[0][2]*v.z + A[0][3],
                  A[1][0]*v.x + A[1][1]*v.y + A[1][2]*v.z + A[1][3],
                  A[2][0]*v.x + A[2][1]*v.y + A[2][2]*v.z + A[2][3]);
  }

}


#endif // AFFINE_H
