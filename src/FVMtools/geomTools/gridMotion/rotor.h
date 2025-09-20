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
#ifndef ROTOR_H
#define ROTOR_H
#include <iostream>
#include <utility>
#include <cstdlib>
#include <vector>
#include "gridTypes.h"

namespace gridMotion {
  typedef vector3d<double> vect3d;
  typedef tensor3d<double> tens3d;
  typedef double real;

  using std::cerr;
  using std::endl;
  using std::rand;
  using std::vector;
  using std::pair;
  using std::make_pair ;
  inline
  double randf() {
    return (double)rand()/(double)RAND_MAX;
  }

  struct Rotor {
    float alpha;
    vector3d<float> beta;

  public:
    Rotor()
      : alpha(1), beta(0,0,0) {
    }
    Rotor(const double a, const vect3d & b)
      : alpha(a), beta(b) {
    }
    /**
     * Input: 3x3 rotation matrix
     *
     * Output: A unit rotor which approximates the rotation component
     * of the given matrix.
     *
     * NOTE: Least squares approximation is utilized to help minimize
     * the error associated with this change in representation.  There
     * are nine equations and four unknowns.
     */
    Rotor(const tens3d & M);
    /**
     * Input: Two vectors.
     *
     * Output: A rotor which will rotate the first vector into the
     * orientation of the second.  More generally, any vector with a
     * component oriented parallel to the plane spanned by these two
     * vectors will have that component rotated in the same direction
     * and by the same angle as between the first and second vectors.
     */
    Rotor(const vect3d & u, const vect3d & v);
    /**
     * Input: Two edge vectors and two displacement vectors.  Here we
     * assume that the two edge vectors have a common origin and when the
     * displacements are added to the edge vectors, we get the two
     * displaced edges, also sharing a common origin.
     */
    Rotor( const pair<vect3d,vect3d> & u,
           const pair<vect3d,vect3d> & v );
    /**
     * Input: three nodes and their displacements.  Here b is the pivot.
     *
     * We transform to a reference frame where b is stationary.  To
     * accomplish this, we subtract the delta-b from all displacements.
     * Next, we determine the axis and angle required to rotate the facet
     * from its initial orientation to its transformed orientation.
     * Finally, we compute the quaternion components which generate this
     * rotation.
     *
     * NOTE: If the angle between the edges should change, then there is
     * no guarantee that the obtained rotor will be able to reproduce both
     * of the desired orientation changes exactly.  This is to be expected
     * since the nodes are obviously not undergoing rigid body motion.
     */
    Rotor(const vect3d & a1, const vect3d & b1, const vect3d & c1,
          const vect3d & da, const vect3d & db, const vect3d & dc) {
      Rotor R(make_pair(a1-b1,da-db), make_pair(c1-b1,dc-db));
      alpha = R.alpha;
      beta = R.beta;
    }
    /**
     * Input: position and displacment of the current node under
     * consideration, and a list of positions and displacments of the
     * adjacent nodes relative to the current node and its displacement.
     * This list consists of all edges adjacent to a given node.  The
     * positions are the locations of the adjacent nodes relative to the
     * current node, and the displacements are the nodal displacements
     * minus the displacement of the current node.
     */
    Rotor(const vector< pair<vect3d, vect3d> > & edgeData);
    /**
     * Input: position and displacment of the current node under
     * consideration, and a list of positions and displacments of the
     * adjacent nodes relative to the current node and its
     * displacement.  This list consists of all edges adjacent to a
     * given node.  The positions are the locations of the adjacent
     * nodes relative to the current node, and the displacements are
     * the nodal displacements minus the displacement of the current
     * node.  This method also includes a reference axis about which
     * the rotation is to be constrained.
     */
    Rotor(const vector< pair<vect3d, vect3d> > & edgeData,
          const vect3d & axis);

  public:
    /**
     * Complement operator
     */
    Rotor operator~() const {
      return Rotor(alpha, -1.0*beta);
    }

  public:
    /**
     * Rotate vector using chem's quaternion implementation
     */
    vect3d operator()(const vect3d & rhs) const {
      gridMotion::Quaternion Q(beta.x, beta.y, beta.z, alpha);
      return Q*rhs;
    }

  public:
    /**
     * Compose two Rotors into one, with assign.
     */
    Rotor & operator*=(const Rotor & rhs) {
      alpha = alpha*rhs.alpha - dot(beta, rhs.beta);
      beta  = alpha*rhs.beta + beta*rhs.alpha + cross(beta, rhs.beta);
      return (*this);
    }
    /**
     * Addition with assign
     */
    Rotor & operator+=(const Rotor & rhs) {
      double dotProd = dot(beta,rhs.beta);
      if (dotProd < 0) {
        alpha -= rhs.alpha;
        beta   = beta - rhs.beta;
      }
      else {
        alpha += rhs.alpha;
        beta   = beta + rhs.beta;
      }
      if (alpha < 0) {
        alpha = -alpha;
        beta  = -1.0*beta;
      }
      return (*this);
    }

  public:
    /**
     * Scale the rotation angle by a constant factor.
     */
    Rotor & scaleAngle(const double s) {
      // Note: With the quaternion representation,
      //        alpha  = cos(theta/2)
      //        |beta| = sin(theta/2).
      const double cosHalfThetaSq = alpha*alpha;
      const double sinHalfThetaSq = dot(beta,beta);
      const double theta(sinHalfThetaSq < 1e-30 ? 0.0 :
                         (cosHalfThetaSq < 1e-30 ? M_PI :
                          2*atan2(sqrt(sinHalfThetaSq),sqrt(cosHalfThetaSq))));
      alpha = cos(0.5*s*theta);
      beta = sin(0.5*s*theta)*(beta/norm(beta));
      return (*this);
    }
    /**
     * Return the rotation angle.
     */
    double angle() const {
      // Note: With the quaternion representation,
      //        alpha  = cos(theta/2)
      //        |beta| = sin(theta/2).
      const double cosHalfThetaSq = alpha*alpha;
      const double sinHalfThetaSq = dot(beta,beta);
      const double theta(sinHalfThetaSq < 1e-30 ? 0.0 :
                         (cosHalfThetaSq < 1e-30 ? M_PI :
                          2*atan2(sqrt(sinHalfThetaSq),sqrt(cosHalfThetaSq))));
      return theta;
    }
    /**
     * Return the rotation axis.
     */
    vect3d axis() const {
      const double sinHalfThetaSq = dot(beta,beta);
      if (sinHalfThetaSq < 1e-30) {
        // Rotation axis is undefined, assume no rotation.
        return vect3d(0,0,0);
      }
      else {
        return beta/norm(beta);
      }
    }
    /**
     * Return the rotation matrix
     */
    tens3d matrix() const {
      const double & a  = alpha;
      const double & bx = beta.x;
      const double & by = beta.y;
      const double & bz = beta.z;
      const tens3d R(vect3d((a*a + bx*bx - by*by - bz*bz),
                            (2*bx*by - 2*a*bz),
                            (2*bx*bz + 2*a*by)),
                     vect3d((2*bx*by + 2*a*bz),
                            (a*a - bz*bz + by*by - bx*bx),
                            (2*by*bz - 2*a*bx)),
                     vect3d((2*bx*bz - 2*a*by),
                            (2*by*bz + 2*a*bx),
                            (a*a - bx*bx - by*by + bz*bz)));
      return R;
    }

  protected:
    /**
     * Input: Axis of rotation, two vectors.
     *
     * Output: A rotor with the prescribed axis of rotation which will
     * rotate the perpendicular component of v1 into the same
     * direction as the perpendicular component of v2.
     */
    Rotor findOrthogonalRotor(const vect3d & axis, const vect3d & v1,
                              const vect3d & v2);
  };

  /**
   * Rotor norm
   */
  inline
  double norm(const Rotor & R) {
    const double lenSq = pow(R.alpha,2) + dot(R.beta,R.beta);
    return ( lenSq > 0 ? sqrt(lenSq) : 0 );
  }

  /**
   * Dot product
   */
  inline
  double dot(const Rotor & lhs, const Rotor & rhs) {
    return (lhs.alpha*rhs.alpha - lhs.beta.x*rhs.beta.x -
            lhs.beta.y*rhs.beta.y - lhs.beta.z*rhs.beta.z);
  }

  /**
   * Compare two rotors - Equality
   */
  inline
  bool operator==(const Rotor & lhs, const Rotor & rhs) {
    // Note that when quaternions are used as rotors, both q and -q
    // will give the same rotations.
    return ( ( ( lhs.alpha  - rhs.alpha  < 1e-12 ) &&
               ( lhs.beta.x - rhs.beta.x < 1e-12 ) &&
               ( lhs.beta.y - rhs.beta.y < 1e-12 ) &&
               ( lhs.beta.z - rhs.beta.z < 1e-12 ) ) ||
             ( ( lhs.alpha  + rhs.alpha  < 1e-12 ) &&
               ( lhs.beta.x + rhs.beta.x < 1e-12 ) &&
               ( lhs.beta.y + rhs.beta.y < 1e-12 ) &&
               ( lhs.beta.z + rhs.beta.z < 1e-12 ) ) );
  }

  /**
   * Compare two rotors - Inequality
   */
  inline
  bool operator!=(const Rotor & lhs, const Rotor & rhs) {
    // Note that when quaternions are used as rotors, both q and -q
    // will give the same rotations.
    return ( ( ( lhs.alpha  - rhs.alpha  >= 1e-12 ) &&
               ( lhs.alpha  + rhs.alpha  >= 1e-12 ) ) ||
             ( ( lhs.beta.x - rhs.beta.x >= 1e-12 ) &&
               ( lhs.beta.x + rhs.beta.x >= 1e-12 ) ) ||
             ( ( lhs.beta.y - rhs.beta.y >= 1e-12 ) &&
               ( lhs.beta.y + rhs.beta.y >= 1e-12 ) ) ||
             ( ( lhs.beta.z - rhs.beta.z >= 1e-12 ) &&
               ( lhs.beta.z + rhs.beta.z >= 1e-12 ) ) );
  }

  /**
   * Scale the Rotor by a constant factor.
   */
  inline
  Rotor operator*(const double lhs, const Rotor & rhs) {
    return Rotor(lhs*rhs.alpha, lhs*rhs.beta);
  }

  /**
   * Scale the Rotor by a constant factor.
   */
  inline
  Rotor operator*(const Rotor & lhs, const double rhs) {
    return Rotor(rhs*lhs.alpha, rhs*lhs.beta);
  }

  /**
   * Divide by a constant factor.
   */
  inline
  Rotor operator/(const Rotor & lhs, const double rhs) {
    return Rotor(lhs.alpha/rhs, lhs.beta/rhs);
  }

  /**
   * Compose two rotations.
   */
  inline
  Rotor operator*(const Rotor & lhs, const Rotor & rhs) {
    return Rotor( lhs.alpha*rhs.alpha - dot(lhs.beta, rhs.beta),
                  lhs.alpha*rhs.beta + lhs.beta*rhs.alpha +
                  cross(lhs.beta, rhs.beta) );
  }

  /**
   * Add two rotors.
   */
  inline
  Rotor operator+(const Rotor & lhs, const Rotor & rhs) {
    double dotProd = dot(lhs.beta, rhs.beta);
    Rotor R;
    if (dotProd < 0) {
      R.alpha = lhs.alpha - rhs.alpha;
      R.beta  = lhs.beta - rhs.beta;
    }
    else {
      R.alpha = lhs.alpha + rhs.alpha;
      R.beta  = lhs.beta + rhs.beta;
    }
    return R;
  }

  /**
   * Difference between two rotations.
   */
  inline
  Rotor operator-(const Rotor & lhs, const Rotor & rhs) {
    double dotProd = dot(lhs.beta, rhs.beta);
    Rotor R;
    if (dotProd < 0) {
      R.alpha = lhs.alpha + rhs.alpha;
      R.beta  = lhs.beta + rhs.beta;
    }
    else {
      R.alpha = lhs.alpha - rhs.alpha;
      R.beta  = lhs.beta - rhs.beta;
    }
    return R;
  }

  /**
   * Apply the rotation represented by the Rotor to the vector.  When
   * the vector is the rhs, the rotation is applied in a positive
   * sense.
   */
  inline
  vect3d operator*(const Rotor & lhs, const vect3d & rhs) {
    const tens3d R = lhs.matrix();
    return dot(R,rhs);
  }

}


#include <iostream>
#include <iomanip>

inline
std::ostream & operator<<(std::ostream & lhs, gridMotion::Rotor & rhs) {
  lhs << rhs.alpha << " " << rhs.beta.x << " " << rhs.beta.y
      << " " << rhs.beta.z;
  return lhs;
}

inline
std::istream & operator>>(std::istream & lhs, gridMotion::Rotor & rhs) {
  lhs >> rhs.alpha >> rhs.beta.x >> rhs.beta.y >> rhs.beta.z;
  return lhs;
}

#endif // ROTOR_H
