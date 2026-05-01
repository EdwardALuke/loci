/** ****************************************************************************
 * @file      mathUtils.h
 * @author    Mark A. Hunt (CFDRC)
 * @brief
 * @date      2022-11-20
 * @copyright CFDRC Copyright (c) 2022
 * @version   0.1
 * @details
 ******************************************************************************/
#ifndef MATH_UTILS_H
#define MATH_UTILS_H


#include <iostream>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort

using namespace std;

/** ****************************************************************************
 * @brief
 * @tparam T
 * @param v
 * @return vector<size_t>
 ******************************************************************************/
template <typename T> vector<size_t> sort_indices(const vector<T> &v)
{
  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  stable_sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2)
                                      {
                                        return v[i1] < v[i2];
                                      });

  return idx;
} // End of sort_indices()


/** ****************************************************************************
 * @brief
 * @param a
 * @param b
 * @param epsilon
 * @return true
 * @return false
 ******************************************************************************/
inline bool bApproximatelyEqual(real_t a, real_t b, real_t epsilon)
{
  return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

/** ****************************************************************************
 * @brief
 * @param a
 * @param b
 * @param epsilon
 * @return true
 * @return false
 ******************************************************************************/
inline bool bEssentiallyEqual(real_t a, real_t b, real_t epsilon)
{
  return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

/** ****************************************************************************
 * @brief
 * @param a
 * @param b
 * @param epsilon
 * @return true
 * @return false
 ******************************************************************************/
inline bool bDefinitelyGreaterThan(real_t a, real_t b, real_t epsilon)
{
  return (a - b) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

/** ****************************************************************************
 * @brief
 * @param a
 * @param b
 * @param epsilon
 * @return true
 * @return false
 ******************************************************************************/
inline bool bDefinitelyLessThan(real_t a, real_t b, real_t epsilon)
{
  return (b - a) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

/** ****************************************************************************
 * @brief
 * @param u
 * @param v
 * @param p
 * @return true
 * @return false
 ******************************************************************************/
inline bool bPointOnEdge(vector3d<real_t> u, vector3d<real_t> v, vector3d<real_t> p)
{
  vector3d<real_t> normal    = cross(v-u, p-u);
  real_t   r1dot2    = dot(  v-u, p-u);
  real_t   rDist1to2 = pow(norm(v - u),2);

  if(norm(normal) < 1.0e-9 && r1dot2 < rDist1to2)
    return true;
  return false;
}

/** ****************************************************************************
 * @brief
 * @param u
 * @param v
 * @return vector3d<real_t>
 ******************************************************************************/
inline vector3d<real_t> v3dVectorProduct(vector3d<real_t> u, vector3d<real_t> v)
{
  return vector3d<real_t>(u.x*v.x, u.y*v.y, u.z*v.z);
}


/** ****************************************************************************
 * @brief
 * @param u
 * @param v
 * @return vector3d<real_t>
 ******************************************************************************/
inline vector3d<real_t> v3dGetNormalToPlane(vector3d<real_t> u, vector3d<real_t> v, vector3d<real_t> w)
{
  return cross(v-u,w-u);
}


/** ****************************************************************************
 * @brief Check if an array of 3D vectors has a particular vector.
 * @param vecU
 * @param v
 * @param epsilon
 * @return true
 * @return false
 ******************************************************************************/
inline bool bHasVector(vector<vector3d<real_t>> vecU, vector3d<real_t> v, real_t epsilon)
{
  int vs = vecU.size();
  for(int iK=0; iK<vs; iK++)
  {
    if(bApproximatelyEqual(vecU[iK].x, v.x, epsilon) &&
       bApproximatelyEqual(vecU[iK].y, v.y, epsilon) &&
       bApproximatelyEqual(vecU[iK].z, v.z, epsilon))
    {
      return true;
    }
  }
  return false;
} // End of bHasVector()

#endif // MATH_UTILS_H
