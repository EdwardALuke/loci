//#############################################################################
//#
//# Copyright 2008-2025, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#include <Loci.h>
#include "vogtools.h"
#include <iostream>
using std::cout ;
using std::cerr ;
using std::endl ;
using std::ios ;
#include <fstream>
using std::ifstream ;
using std::ofstream ;

#include <utility>
#include <algorithm>
using std::swap ;
using std::max ;
using std::min ;
using std::max_element ;
using std::min_element ;
using std::sort ;
using std::reverse ;
using std::unique ;
using std::rotate ;
using std::find ;
using std::transform ;
#include <string>
using std::string ;
#include <vector>
using std::vector ;
#include <list>
using std::list ;
#include <map>
using std::map ;
#include <functional>
using std::plus ;
//using std::bind2nd ;
#include <limits>

#ifdef NO_CMATH
#include <math.h>
#else
#include <cmath>
#endif

#ifdef NO_CSTDLIB
#include <stdlib.h>
#include <ctype.h>
#else
#include <cstdlib>
#include <cctype>
#endif

#include <sys/time.h>

#include <sstream>
using std::stringstream ;


/*************************************************************
 * include for use of J. R. Shewchuk's robust geometric      *
 * predicates and one additional predicate by myself         *
 *************************************************************/
#define REAL double
#ifdef __cplusplus
extern "C" {
#endif
REAL orient2dfast(REAL*, REAL*, REAL*) ;
REAL orient2dexact(REAL*, REAL*, REAL*) ;
REAL orient2d(REAL*, REAL*, REAL*) ;
REAL cross2dfast(REAL*, REAL*, REAL*, REAL*) ;
REAL cross2dexact(REAL*, REAL*, REAL*, REAL*) ;
void exactinit() ;
#ifdef __cplusplus
}
#endif


/**************************************************************
 *                      CODE BEGINS HERE                      *
 **************************************************************/
//#define DEBUG
#define NORMAL
//#define FACE_SPLIT_SELF_CHECK
//#define DEBUG_SPLIT_SHOWME
//#define DEBUG_POLY_AND_POLY
//#define SHOW_QT_PROPERTY
//#define FACE_SPLIT_PROPERTY_CHECK
//#define UNSTABLE_SEG_SEG_INT_CODE
//#define POINT_SHIFT_CONVEXITY_CHECK
//#define POLY_AND_POLY_CONVEXITY_CHECK

typedef vector3d<double> vec3d ;
typedef vector2d<double> vec2d ;
typedef vector3d<vec3d> ten3d ;

/**
***************************************************************
*                    Some Global Variables                    *
**/
inline double operator-(const timeval& t1, const timeval& t2) {
  double dt1 = t1.tv_sec + t1.tv_usec*1e-6 ;
  double dt2 = t2.tv_sec + t2.tv_usec*1e-6 ;
  return dt1-dt2 ;
}
// variables for timing report
timeval time_prog_start ;
timeval time_prog_end ;

timeval time_essential_start ;
timeval time_essential_end ;
double prog_essential_time = 0 ;

timeval time_grid_read_start ;
timeval time_grid_read_end ;

timeval time_data_structure_start ;
timeval time_data_structure_end ;
double data_structure_time = 0 ;

timeval time_shift_point_total_start ;
timeval time_shift_point_total_end ;
timeval time_shift_point_start ;
timeval time_shift_point_end ;
timeval time_shift_point_qt_start ;
timeval time_shift_point_qt_end ;

timeval time_face_split_total_start ;
timeval time_face_split_total_end ;
timeval time_face_split_start ;
timeval time_face_split_end ;
timeval time_face_split_qt_start ;
timeval time_face_split_qt_end ;
timeval time_face_split_frm_start ;
timeval time_face_split_frm_end ;

timeval time_BC1_reconstruct_start ;
timeval time_BC1_reconstruct_end ;

timeval time_BC2_reconstruct_start ;
timeval time_BC2_reconstruct_end ;

timeval time_boundary_nodes_matching_start ;
timeval time_boundary_nodes_matching_end ;

timeval time_3D_to_2D_projection_start ;
timeval time_3D_to_2D_projection_end ;

timeval time_2D_to_3D_projection_start ;
timeval time_2D_to_3D_projection_end ;

timeval time_new_grid_write_start ;
timeval time_new_grid_write_end ;

// usually set in range (0.0,1.0)
double SHIFT_RANGE_THRESHOLD = 0.05 ;
const double PI =
3.141592653589793238462643383279502884197169399375105820974944592308 ;
inline double deg2rad(double d) {
  return d * PI / 180.0 ;
}
// flag used to supress terminal output in functions other than main()
bool global_verbose = true ;

// the xor defined for bool
inline bool xor_bool(bool a, bool b) {
  return ( (a^b) == 1) ;
}
/**
***************************************************************
*                    Some vector utilities                    *
**/
// test if a vec2d is (0,0)
inline bool zero_vec2d(const vec2d& v) {
  return ( (v.x==0.0) && (v.y==0.0)) ;
}

const double NORMALIZE_ZERO_THRESHOLD = 0.0 ;
// normalize the given vector
inline void normalize(vec2d& v) {
  // zero vectors testing 
  if( (fabs(v.x) <= NORMALIZE_ZERO_THRESHOLD) &&
      (fabs(v.y) <= NORMALIZE_ZERO_THRESHOLD)) {
    cerr << "WARNING: normalizing zero vector, nothing done!" << endl ;
    return ;
  }
  double t = sqrt(v.x*v.x + v.y*v.y) ;
  v.x /= t ;
  v.y /= t ;
}

inline void normalize(vec3d& v) {
  if( (fabs(v.x) <= NORMALIZE_ZERO_THRESHOLD) &&
      (fabs(v.y) <= NORMALIZE_ZERO_THRESHOLD) &&
      (fabs(v.z) <= NORMALIZE_ZERO_THRESHOLD)) {
    cerr << "WARNING: normalizing zero vector, nothing done!" << endl ;
    return ;
  }
  double t = sqrt(v.x*v.x + v.y*v.y + v.z*v.z) ;
  v.x /= t ;
  v.y /= t ;
  v.z /= t ;
}

// distance between two points
inline double dist(const vec2d& v1, const vec2d& v2) {
  return norm(v1-v2) ;
}
inline double dist(const vec3d& v1, const vec3d& v2) {
  return norm(v1-v2) ;
}
// define the == operator for vec2d
inline bool operator==(const vec2d& v1, const vec2d& v2) {
  return ( (v1.x==v2.x) && (v1.y==v2.y)) ;
}
// define the != operator for vec2d
inline bool operator!=(const vec2d& v1, const vec2d& v2) {
  return !(v1==v2) ;
}
/**
***************************************************************
*                     Some global Threshold                   *
**/
// this is the threshold used in comparing points in the grid.
// any two point below this threshold will be considered the
// same. NOTE: this threshold will be reset when we get the
// the minimum edge length in both grid. It will be set as
// 1e-12 * the minimum edge lentgh
// NOTE: this threshold is NOT in active use
double GRID_POINT_COMPARE_THRESHOLD = 1e-12 ;
// this threshold controls the amount of tolerance in checking
// the area sum after split. Basically it means that the
// tolerated area difference is to be within
// boundary area / AREA_SUM_CHECK_THRESHOLD
double AREA_SUM_CHECK_THRESHOLD = 1e8 ;
// this is the threshold for on line segment test in
// function point_on_edge (the dot product)
// 1e-12 is about +/-(1e-4) degree difference
//const double ON_LINE_SEG_DOT_THRESHOLD = 1e-12 ;
// this value is approximately +/- 1e-3 degree
//static double ON_LINE_SEG_DOT_THRESHOLD = 2e-10 ;
static double ON_LINE_SEG_DOT_THRESHOLD = 2e-5 ;

/************************************************
 * Some floating point functions based on       *
 * threshold values. They are not used directly *
 * inside the face splitting routine, hence     *
 * they will not affect the robustness of the   *
 * face splitting. Instead, they are used in    *
 * the pre-processing step such as point        *
 * shifting, or for the approximate geometric   *
 * predicates for boundary outer contour        *
 * geometry processing.                         *
 ***********************************************/
inline bool pointsEqual(const vec2d& v1, const vec2d& v2) {
  return (v1 == v2) ;
}
// this function determines the sign of a-b based on the
// passed in threshold. It returns 0 if a b are equal,
// +1, if a > b, and -1 if a < b
// NOTE: threshold should always be positive
inline int twoDoubleNumSign(const double& a, const double& b,
                            const double& threshold) {
  double d = a-b ;
  if(abs(d) <= threshold)
    return 0 ;
  else if(d > 0)
    return 1 ;
  else // d < 0
    return -1 ;
}
// the followings are double number comparison functions
// that based on the passed in threshold
inline bool doubleNumEqual(const double& d1, const double& d2,
                           const double& threshold) {
  return (twoDoubleNumSign(d1,d2,threshold) == 0) ;
}
inline bool doubleNumGreat(const double& d1, const double& d2,
                           const double& threshold) {
  return (twoDoubleNumSign(d1,d2,threshold) == 1) ;
}
inline bool doubleNumLess(const double& d1, const double& d2,
                          const double& threshold) {
  return (twoDoubleNumSign(d1,d2,threshold) == -1) ;
}

/**
*******************
* utility functions for floating point number summation 
**/
inline bool abs_double_descend(const double& d1, const double& d2) {
  return (abs(d1) > abs(d2)) ;
}
inline double addWithError(double& err, const double& a, const double& b) {
  double sum = a+b;
  err = b-(sum-a);
  return sum;
}
// summation using "doubly compensated summation algorithm"
// NOTE: this function changes the passed in vector!!!
// We could use without reference here, but for efficiency
// concern, we would use reference
double double_summation(vector<double>& v)
{
  if(v.empty())
    return 0. ;
  sort(v.begin(),v.end(),abs_double_descend) ;
  double sum = v[0];
  double sum_err = 0.;
  vector<double>::size_type size = v.size() ;
  for(vector<double>::size_type i=1;i<size;++i) {
    double addend_err;
    double addend = addWithError(addend_err, sum_err, v[i]);
    double sum1_err;
    double sum1 = addWithError(sum1_err, sum, addend);
    sum = addWithError(sum_err, sum1, addend_err + sum1_err);
  }
  return sum ;
}

/***********************************************************
 * Here are some robust, non-robust, and approximate       *
 * geometric predicates. Most of the robust one uses       *
 * underline predicates from the predicates.c file.        *
 ***********************************************************/
// this function returns if point p is to the left of
// of the segment (head,tail). NOTE: the segment is from
// tail to head. As the vector head-tail
// this is a non-robust version
inline bool left_fast(const vec2d& head, const vec2d& tail,
                      const vec2d& p) {
  vec2d seg_v = head - tail ;
  vec2d point_v = p - tail ;
  return (cross(seg_v,point_v) > 0.0) ;
}
// the robust version
inline bool left_exact(const vec2d& head, const vec2d& tail,
                       const vec2d& p) {
  double hc[2] = {head.x, head.y} ;
  double tc[2] = {tail.x, tail.y} ;
  double pc[2] = {p.x, p.y} ;
  return (orient2d(tc,hc,pc) > 0.0) ;
}
// the non-robust, but fast version of left on.
// all the definitions are the same as above.
inline bool leftOn_fast(const vec2d& head, const vec2d& tail,
                        const vec2d& p) {
  vec2d seg_v = head - tail ;
  vec2d point_v = p - tail ;
  return (cross(seg_v,point_v) >= 0.0) ;
}
// the robust version
inline bool leftOn_exact(const vec2d& head, const vec2d& tail,
                         const vec2d& p) {
  double hc[2] = {head.x, head.y} ;
  double tc[2] = {tail.x, tail.y} ;
  double pc[2] = {p.x, p.y} ;
  return (orient2d(tc,hc,pc) >= 0.0) ;
}
// returns the sign of the cross product of two vectors
// defined by four points. v1b, v1e define vector v1=v1b-v1e
// v2b, v2e define vector v2=v2b-v2e. This function returns
// +1 if cross(v1,v2) is positive; -1 if cross(v1,v2) is
// negative; and 0 if cross(v1,v2) is zero.
// this is a non-robust version
inline int cross2d_fast(const vec2d& v1b, const vec2d& v1e,
                        const vec2d& v2b, const vec2d& v2e) {
  vec2d v1=v1b-v1e ;
  vec2d v2=v2b-v2e ;
  double c = cross(v1,v2) ;
  if(c > 0.0)
    return 1 ;
  else if(c < 0.0)
    return -1 ;
  else // c == 0.0
    return 0 ;
}
// this is the robust (but not yet adaptive) version
inline int cross2d_exact(const vec2d& v1b, const vec2d& v1e,
                         const vec2d& v2b, const vec2d& v2e) {
  double a[2] = {v1b.x, v1b.y} ;
  double b[2] = {v1e.x, v1e.y} ;
  double c[2] = {v2b.x, v2b.y} ;
  double d[2] = {v2e.x, v2e.y} ;
  double sign = cross2dexact(a,b,c,d) ;
  if(sign > 0.0)
    return 1 ;
  else if(sign < 0.0)
    return -1 ;
  else // sign == 0.0
    return 0 ;
}
// predicate for collinearity of point a, b, and c
// this is the fast, but non-robust version
inline bool collinear3_fast(const vec2d& a, const vec2d& b, const vec2d& c) {
  vec2d v1=b-a ;
  vec2d v2=c-a ;
  return (cross(v1,v2) == 0.0) ;
}
// the robust version
inline bool collinear3_exact(const vec2d& a, const vec2d& b,
                             const vec2d& c) {
  double pa[2] = {a.x, a.y} ;
  double pb[2] = {b.x, b.y} ;
  double pc[2] = {c.x, c.y} ;
  return (orient2d(pa,pb,pc) == 0.0) ;
}
// this is the approximate version. In that, we allow an error
// threshold of the collinearity. We compute the largest angle
// of the triangle formed by three points passed in. If the
// largest angle is within 180degree +/- ON_LINE_SEG_DOT_THRESHOLD
// they are still considered to be collinear
inline bool collinear3_approx(const vec2d& a, const vec2d& b,
                              const vec2d& c) {
  // we'll need to first test if any two points coincide, if
  // true, then they are surely collinear
  if( (a==b) || (a==c) || (b==c))
    return true ;
  // we determine the collinearity by find the largest angle in
  // the triangle defined by a,b,c. The largest angle should
  // be close to 180 degree. Thus this equivalents to compute
  // the minimum cosine value of the three internal angles.
  // We do this by computing the dot product.
  vec2d v1, v2 ;
  // first angle A
  v1=b-a ; v2=c-a ; normalize(v1) ; normalize(v2) ;
  double da = dot(v1,v2) ;
  // then angle B
  v1=a-b ; v2=c-b ; normalize(v1) ; normalize(v2) ;
  double db = dot(v1,v2) ;
  // finally angle C
  v1=a-c ; v2=b-c ; normalize(v1) ; normalize(v2) ;
  double dc = dot(v1,v2) ;
  // find the minimum
  double d = min(min(da,db),dc) ;
  // if d is -1, then a,b,c are collinear
  return (doubleNumEqual(d,-1.0,ON_LINE_SEG_DOT_THRESHOLD)) ;
}
// this function returns whether the point c lies on
// the segment defined by a and b. it is assumed that
// a, b, c are collinear already
inline bool between(const vec2d& a, const vec2d& b, const vec2d& c) {
  // if a b not vertical, the check x value, otherwise check y value
  if(a.x!=b.x)
    return ( ((a.x <= c.x) && (b.x >= c.x)) ||
             ((a.x >= c.x) && (b.x <= c.x))) ;
  else
    return ( ((a.y <= c.y) && (b.y >= c.y)) ||
             ((a.y >= c.y) && (b.y <= c.y))) ;
}
// determin if the given point lies on the line segments defined
// by two points. This is the exact version. Formally, it computes
// whether the given point is in [line_start,line_end]
inline bool pointOnEdge_exact(const vec2d& p,
                              const vec2d& line_start,
                              const vec2d& line_end) {
  if(!collinear3_exact(line_end,line_start,p))
    return false ;
  return between(line_start,line_end,p) ;
}
// this is the approximate version. It computers the angle
// by segments (p->line_start) and (p->line_end). If it is
// nearly 180 degree, then we consider p on the edge
inline bool pointOnEdge_approx(const vec2d& p,
                               const vec2d& line_start,
                               const vec2d& line_end) {
  // first check if p coincides with either line_start or line_end
  if( (p == line_start) || (p == line_end))
    return true ;
  vec2d v1 = line_start - p ;
  vec2d v2 = line_end - p ;
  normalize(v1) ; normalize(v2) ;
  double a = dot(v1,v2) ;
  // a should be -1 if p is on the line segments
  if(doubleNumEqual(a,-1.0,ON_LINE_SEG_DOT_THRESHOLD))
    return true ;
  else
    return false ;
}

typedef enum {POS_NOTPAR,POS_PAR,POS_COLLINEAR} seg_seg_pos_type ;
// computes the relative position of two line segments.
// possible results are Not parallel, parallel but not
// collinear, or collinear. The line segments are defined
// as (l1s,l1e), (l2s,l2e), l1s and l2s are heads,
// and l1e, l2e are tails. the direction is from
// tail to head.

// this is the approximate version
inline seg_seg_pos_type
seg_seg_pos_approx(const vec2d& l1s, const vec2d& l1e,
                   const vec2d& l2s, const vec2d& l2e) {
  bool l1s_collinear = collinear3_approx(l1s,l2s,l2e) ;
  bool l1e_collinear = collinear3_approx(l1e,l2s,l2e) ;
  
  if(l1s_collinear && l1e_collinear) {
    return POS_COLLINEAR ;
  } else if(!l1s_collinear && l1e_collinear) {
    return POS_NOTPAR ;
  } else if(l1s_collinear && !l1e_collinear) {
    return POS_NOTPAR ;
  } else {
    // neither l1s nor l1e collinear with segment (l2s,l2e)
    // therefore there are two possibilities: they are parallel
    // or they are not parallel, thus we compute the cross2d
    // if it is zero, then they are parallel, otherwise, they
    // are not parallel
    if(cross2d_fast(l1s,l1e,l2s,l2e) == 0)
      return POS_PAR ;
    else
      return POS_NOTPAR ;
  }
}
// this is the exact version
inline seg_seg_pos_type
seg_seg_pos_exact(const vec2d& l1s, const vec2d& l1e,
                  const vec2d& l2s, const vec2d& l2e) {
  if(cross2d_exact(l1s,l1e,l2s,l2e) == 0) {
    // they could be parallel or collinear
    if(collinear3_exact(l1s,l2s,l2e))
      return POS_COLLINEAR ;
    else
      return POS_PAR ;
  }else
    return POS_NOTPAR ;
}
/**************the end of predicates************************/

// a global table for storing betweenness of points.
// i.e., if point p lies on edge defined by a and b (this also means
// that p,a,b are collinear), then record (p,a,b) will be
// in the table. when storing the record (p,a,b), (a,b) will
// be sorted in increasing order because that saves recording
// space. Because if p lies between a,b, so does b,a.
// The semantic of "on edge" is the same as that in the
// pointOnEdge_exact predicate
// This table is mainly created for handling the inexactness
// caused by the point-->edge shifting process.
class point_on_edge_lookup_table {
  // levels of lookup table
  struct final {
    entitySet domain ;
    map<Entity,entitySet> record ;
    bool inSet(Entity id) const {return domain.inSet(id) ;}
    void clear() {domain=EMPTY ; record.clear() ;}
    void add_record(Entity lb, Entity le) {
      if(!inSet(lb))
        domain += lb ;
      record[lb] += le ;
    }
    bool has_record(Entity lb, Entity le) const {
      if(!inSet(lb)) // no record yet
        return false ;
      else {
        // we'll need to see whether le is in or not
        return (record.find(lb)->second).inSet(le) ;
      }
    }
  } ;
  
  struct first {
    entitySet domain ;
    map<Entity,final> record ;
    bool inSet(Entity id) const {return domain.inSet(id) ;}
    void clear() {domain=EMPTY ; record.clear() ;}
    void add_record(Entity p, Entity lb, Entity le) {
      if(!inSet(p)) {
        domain += p ;
      }
      record[p].add_record(lb,le) ;
    }
    bool has_record(Entity p, Entity lb, Entity le) const {
      if(inSet(p))
        return (record.find(p)->second).has_record(lb,le) ;
      else
        return false ;
    }
  } ;

  // the table data
  first table ;
  // The following are public interfaces
public:
  void clear() {
    table.clear() ;
  }
  void add_record(Entity p, Entity lb, Entity le) {
    // first sort segment's id
    Entity new_lb,new_le ;
    new_lb = min(lb,le) ; new_le = max(lb,le) ;
    table.add_record(p,new_lb,new_le) ;
  }
  bool has_record(Entity p, Entity lb, Entity le) const {
    // first sort segment's id
    Entity new_lb,new_le ;
    new_lb = min(lb,le) ; new_le = max(lb,le) ;
    return table.has_record(p,new_lb,new_le) ;
  }
} ;
// now we declare a global table for point betweenness query
point_on_edge_lookup_table global_pbetween_ltable ;
//////////end of point betweenness lookup table////////////

/**
 * Following are some geometric predicates based on the
 * point betweenness lookup table. In that, they would
 * first query the lookup table, if the results are already
 * in, then just return what is found. Otherwise, call the
 * exact predicates to predict the result and possibly
 * store the results back in. NOTE: all such routines end
 * with the suffix _tl
 **/
// The first is the collinear3 test
inline bool collinear3_tl(const vec2d& a, const vec2d& b, const vec2d& c,
                          Entity aid, Entity bid, Entity cid) {
  // first lookup the table
  if(global_pbetween_ltable.has_record(cid,aid,bid))
    return true ;
  // otherwise call collinear3_exact
  return collinear3_exact(a,b,c) ;
}
inline bool collinear3_approx_tl
(const vec2d& a, const vec2d& b, const vec2d&c,
 Entity aid, Entity bid, Entity cid) {
  if(global_pbetween_ltable.has_record(cid,aid,bid))
    return true ;
  return collinear3_approx(a,b,c) ;
}
// then the table lookup version of between
inline bool between_tl(const vec2d& a, const vec2d& b, const vec2d& c,
                       Entity aid, Entity bid, Entity cid) {
  // first lookup the table
  if(global_pbetween_ltable.has_record(cid,aid,bid))
    return true ;
  // otherwise call between
  return between(a,b,c) ;
}
// This one is table lookup version of the pointOnEdge
inline bool pointOnEdge_tl(const vec2d& p, const vec2d& lb, const vec2d& le,
                           Entity pid, Entity lb_id, Entity le_id) {
  if(global_pbetween_ltable.has_record(pid,lb_id,le_id))
    return true ;
  if(pointOnEdge_exact(p,lb,le)) {
    global_pbetween_ltable.add_record(pid,lb_id,le_id) ;
    return true ;
  }else
    return false ;
}
// this is the table lookup version of the segment
// position predicate
inline seg_seg_pos_type
seg_seg_pos_tl(const vec2d& l1s, const vec2d& l1e,
               const vec2d& l2s, const vec2d& l2e,
               Entity l1s_id, Entity l1e_id,
               Entity l2s_id, Entity l2e_id) {
  bool l1s_collinear = collinear3_tl(l1s,l2s,l2e,l1s_id,l2s_id,l2e_id) ;
  bool l1e_collinear = collinear3_tl(l1e,l2s,l2e,l1e_id,l2s_id,l2e_id) ;

  if(l1s_collinear && l1e_collinear) {
    return POS_COLLINEAR ;
  } else if(!l1s_collinear && l1e_collinear) {
    return POS_NOTPAR ;
  } else if(l1s_collinear && !l1e_collinear) {
    return POS_NOTPAR ;
  } else {
    // neither l1s nor l1e collinear with segment (l2s,l2e)
    // therefore there are two possibilities: they are parallel
    // or they are not parallel, thus we compute the cross2d
    // if it is zero, then they are parallel, otherwise, they
    // are not parallel
    if(cross2d_exact(l1s,l1e,l2s,l2e) == 0)
      return POS_PAR ;
    else
      return POS_NOTPAR ;
  }
}

/**
 ***********************************************************
 *               Some utility functions                    *
 ***********************************************************
 */

// followings are utilities for number, string conversion
// and integer and floating number parsing
template<typename T>
inline string num2str(T n) {
  stringstream ss ;
  ss << n ;
  return ss.str() ;
}

inline int chars2int(char* s) {
  stringstream ss ;
  ss << s ;
  int ret ;
  ss >> ret ;

  return ret ;
}

inline double chars2double(char* s) {
  stringstream ss ;
  ss << s ;
  double ret ;
  ss >> ret ;

  return ret ;  
}

inline int str2int(string s) {
  stringstream ss ;
  ss << s ;
  int ret ;
  ss >> ret ;

  return ret ;
}

inline double str2double(string s) {
  stringstream ss ;
  ss << s ;
  double ret ;
  ss >> ret ;

  return ret ;  
}

// determine whether range [b,e) contains only digits
// [b,e) must be a valid range and contains characters
template <class ForwardIter>
inline bool alldigits(ForwardIter b, ForwardIter e) {
  for(;b!=e;++b)
    if(!isdigit(*b))
      return false ;

  return true ;
}
// determine whether s contains a valid integer (including the sign)
inline bool valid_int(const string& s) {
  if(s.empty())
    return false ;
  if( (s[0] == '-') || (s[0] == '+')) {
    if(s.size() > 1)
      return alldigits(s.begin()+1,s.end()) ;
    else
      return false ;
  }else
    return alldigits(s.begin(),s.end()) ;
}
// determine whether s contains a valid double (including the sign)
inline bool valid_double(const string& s) {
  if(s.empty())
    return false ;
  // find the "." place
  string::const_iterator si = find(s.begin(),s.end(),'.') ;
  // if no "." then we test for things before "e" or "E" if has
  if(si == s.end())
    si = find(s.begin(),s.end(),'e') ;
  if(si == s.end())
    si = find(s.begin(),s.end(),'E') ;
  if(!valid_int(string(s.begin(),si)))
    return false ;
  
  // find the "e" or "E" if any
  string::const_iterator si2 = find(s.begin(),s.end(),'e') ;
  if(si2 == s.end())
    si2 = find(s.begin(),s.end(),'E') ;
  if(si<si2)
    if(!alldigits(si+1,si2))
      return false ;
  if(si2 != s.end())
    if(!valid_int(string(si2+1,s.end())))
      return false ;
  
  return true ;
}

// distance between a point and a line through two points
inline double point2lineDist(const vec2d& p,
                             const vec2d& line_start,
                             const vec2d& line_end) {
  // first we check if two points that define the line coincide
  // this degenerates to the distance of two points
  if(line_start==line_end) {
    return dist(p,line_start) ;
  }
  double line_mag = dist(line_end,line_start) ;
  double u = ( (p.x - line_start.x)*(line_end.x - line_start.x) +
               (p.y - line_start.y)*(line_end.y - line_start.y))
    / (line_mag*line_mag) ;
  vec2d intersect = line_start + u*(line_end-line_start) ;
  return dist(p,intersect) ;
}
// computes the intersection of the given line segment (ls)
// and a line passes the given point p and perpendicular to ls.
// return true if the intersection is within ls, false
// if intersection is out of ls range. When return true,
// it also set the intersection point and the distance
inline bool point2lineInt(const vec2d& p,
                          const vec2d& line_start,
                          const vec2d& line_end,
                          vec2d& intersect,
                          double& distance) {
  // first we check if two points that define the line coincide
  // if so, we return false ;
  if(line_start == line_end) {
    return false ;
  }
  double line_mag = dist(line_end,line_start) ;
  double u = ( (p.x - line_start.x)*(line_end.x - line_start.x) +
               (p.y - line_start.y)*(line_end.y - line_start.y))
    / (line_mag*line_mag) ;
  if(doubleNumLess(u,0.0,1e-12) || // if((u<0.0)||(u>1.0))
     doubleNumGreat(u,1.0,1e-12))
    return false ;
  intersect = line_start + u*(line_end-line_start) ;
  distance = dist(p,intersect) ;
  return true ;
}

// given a vector in 3d space, computes two orthogonal of its vectors
// adapted from the metrics.cc (orthogonal_coords class) in Chem/src
void orthogonal_coords(const vec3d& n, vec3d& u, vec3d& v) {
  // Attempt to minimize cancellation error when finding orthogonal vector.
  // Find the coordinate direction which is largest and base orthognality
  // on that direction. (The three following cases are for x largest,
  // then y largest, and finally z largest.
  if(abs(n.x)>abs(n.y) && abs(n.x)>abs(n.z)) {
    if(abs(n.y)>abs(n.z)) {
      u.y = n.y ;
      u.z = -n.x ;
    } else {
      u.y = -n.x ;
      u.z = n.z ;
    }
    u.x = -(u.y * n.y + u.z * n.z) / n.x ;
  } else if(abs(n.y)>abs(n.x) && abs(n.y)>abs(n.z)) {
    if(abs(n.x)>abs(n.z)) {
      u.x = n.x ;
      u.z = -n.y ;
    } else {
      u.x = -n.y ;
      u.z = n.z ;
    }
    u.y = -(u.x * n.x + u.z * n.z) / n.y ;
  } else {
    if(abs(n.x)>abs(n.y)) {
      u.x = n.x ;
      u.y = -n.z ;
    } else {
      u.x = -n.z ;
      u.y = n.y ;
    }
    u.z = -(u.x * n.x + u.y * n.y) / n.z ;
  }
  
  const double usr = 1./sqrt(dot(u,u)) ;
  u *= usr ;  //normalize the vector
  v = cross(u,n) ;  // compute v at last
}

// do an axis-radius projection. axis is the "axis",
// u v is the radius plane, orig is the origin of rotation
// npos is the set of nodes' 3d position that needs to be projected
// NOTE: u,v,axis should all be normalized.
void axis_radius_projection(const vec3d& axis,
                            const vec3d& u,
                            const vec3d& v,
                            const vec3d& orig,
                            const store<vec3d>& npos,
                            store<vec2d>& proj_pos) {
  proj_pos.allocate(EMPTY) ;
  entitySet nodes = npos.domain() ;
  proj_pos.allocate(nodes) ;

  for(entitySet::const_iterator ei=nodes.begin();
      ei!=nodes.end();++ei) {
    vec3d pv = npos[*ei] - orig ;
    double x = dot(pv,axis) ;
    double u_proj = dot(pv,u) ;
    double v_proj = dot(pv,v) ;
    double r = sqrt(u_proj*u_proj + v_proj*v_proj) ;
    proj_pos[*ei] = vec2d(x,r) ;
  }
}
// do an orthogonal projection
// u, v is the projection plane, npos is the set
// of points to be projected.
// NOTE: u,v should normalized.
void orthogonal_projection(const vec3d& u,
                           const vec3d& v,
                           const store<vec3d>& npos,
                           store<vec2d>& proj_pos) {
  proj_pos.allocate(EMPTY) ;
  entitySet nodes = npos.domain() ;
  proj_pos.allocate(nodes) ;

  for(entitySet::const_iterator ei=nodes.begin();
      ei!=nodes.end();++ei) {
    const vec3d& pv = npos[*ei] ;
    double u_proj = dot(pv,u) ;
    double v_proj = dot(pv,v) ;
    proj_pos[*ei] = vec2d(v_proj,u_proj) ;
  }
}

// get all the nodes on a boundary
entitySet get_boundary_nodes(const multiMap& face2node,
                                    const entitySet& bc_faces) {
  entitySet nodes ;
  for(entitySet::const_iterator fi=bc_faces.begin();
      fi!=bc_faces.end();++fi) {
    int face_node_num = face2node.end(*fi) - face2node.begin(*fi) ;
    for(int i=0;i<face_node_num;++i)
      nodes += face2node[*fi][i] ;
  }
  return nodes ;
}

// get the position of the passed in nodes
Loci::storeRepP get_nodes_pos(const store<vec3d>& pos,
                                  const entitySet& nodes) {
  store<vec3d> subpos ;
  subpos.allocate(nodes) ;
  for(entitySet::const_iterator ei=nodes.begin();
      ei!=nodes.end();++ei)
    subpos[*ei] = pos[*ei] ;

  return subpos.Rep() ;
}

// get all the elements from the passed in index in a multiMap
entitySet get_multiMap_elems(const multiMap& m, int index) {
  entitySet ret ;
  int num_elems = m.num_elems(index) ;
  for(int i=0;i<num_elems;++i)
    ret += m[index][i] ;

  return ret ;
}
// get all teh elements in order from the passed in index in a multiMap
vector<int> get_multiMap_elems_inorder(const multiMap& m,
                                       int index) {
  vector<int> elems ;
  int num_elems = m.num_elems(index) ;
  for(int i=0;i<num_elems;++i)
    elems.push_back(m[index][i]) ;
  return elems ;
}
// given a face and face2nodes map, get all the nodes in order
vector<int> get_face_nodes(const multiMap& face2node,
                           int face_id) {
  return get_multiMap_elems_inorder(face2node,face_id) ;
}

// get all the edges on the passed in boundary mesh
void get_edges_info(const multiMap& face2node,
                    const multiMap& node2face,
                    const entitySet& bc_faces,
                    const Map& cr, dMap& N1, dMap& N2,
                    dMap& El, dMap& Er,
                    int& edge_num) {
  edge_num = 0 ;
  store<bool> visited ;
  visited.allocate(face2node.domain()) ;
  for(entitySet::const_iterator ei=visited.domain().begin();
      ei!=visited.domain().end();++ei)
    visited[*ei] = false ;

  for(entitySet::const_iterator ei=bc_faces.begin();
      ei!=bc_faces.end();++ei) {
    visited[*ei] = true ;
    // get the number of nodes on this face
    int nf = face2node.end(*ei) - face2node.begin(*ei) ;
    // loop over all the edges on this face
    for(int i=0;i<nf;++i) {
      // get the nodes sequentially
      int n1 = face2node[*ei][(i==0)?(nf-1):(i-1)] ;
      int n2 = face2node[*ei][i] ;
      // get the two faces that share the edge (n1 n2)
      entitySet n1_faces = get_multiMap_elems(node2face,n1) ;
      entitySet n2_faces = get_multiMap_elems(node2face,n2) ;
      // intersection will give us the faces
      entitySet neighbor_faces = n1_faces & n2_faces ;
      neighbor_faces &= bc_faces ;
      int nfsize = neighbor_faces.size() ;
      // there should only have two faces
      if( (nfsize == 0) || (nfsize > 2)) {
        cerr << "Data inconsistent in getting edges..." << endl ;
        cerr << "neighbor_faces.size() = " << nfsize << endl ;
        Loci::Abort() ;
      }
      if(nfsize == 1) {
        // then this edge must sit on the boundary
        Entity f = *(neighbor_faces.begin()) ;
        // and we check this
        if(f != *ei) {
          cerr << "Data inconsistent in getting edges..." << endl ;
          Loci::Abort() ;
        }
        // then collect the data
        El[edge_num] = *ei ;
        Er[edge_num] = cr[*ei] ;
        N1[edge_num] = n1 ;
        N2[edge_num] = n2 ;
        ++edge_num ;
      }else {
        // this is an internal edge
        // we record it if it is not seen before
        neighbor_faces -= *ei ;
        if(neighbor_faces.size() != 1) {
          cerr << "Data inconsistent in getting edges..." << endl ;
          Loci::Abort() ;
        }
        Entity f = *(neighbor_faces.begin()) ;
        if(!visited[f]) {
          El[edge_num] = *ei ;
          Er[edge_num] = f ;
          N1[edge_num] = n1 ;
          N2[edge_num] = n2 ;
          ++edge_num ;
        }
      }
    } // end of for(edges)
  } // end of for(bc_faces)
}

// output projections in 2dgv format
void TwodgvOutput(const char* fname, const store<vec2d>& npos,
                  const dMap& N1, const dMap& N2,
                  const dMap& El, const dMap& Er,
                  int edge_num, const entitySet& bc_faces) {
  
  std::ofstream out(fname,ios::out) ;
  out.precision(16) ;
  out << "general" << endl ;

  entitySet npos_dom = npos.domain() ;
  int node_num = npos_dom.size() ;
  out << node_num << ' ' << '1' << endl ;

  for(entitySet::const_iterator ei=npos_dom.begin();
      ei!=npos_dom.end();++ei)
    out << npos[*ei].x << ' ' << npos[*ei].y << endl ;

  // then build a node and face index store
  store<int> node_indx ;
  node_indx.allocate(npos_dom) ;
  int cnt = 1 ;
  for(entitySet::const_iterator ei=npos_dom.begin();
      ei!=npos_dom.end();++ei,++cnt)
    node_indx[*ei] = cnt ;
  
  store<int> face_indx ;
  face_indx.allocate(bc_faces) ;
  cnt = 1 ;
  for(entitySet::const_iterator ei=bc_faces.begin();
      ei!=bc_faces.end();++ei,++cnt)
    face_indx[*ei] = cnt ;

  // then output edges
  out << edge_num << " 1 "  << cnt << " 1" << endl ;
  for(int i=0;i<edge_num;++i) {
    out << node_indx[N1[i]] << ' '
        << node_indx[N2[i]] << ' '
        << face_indx[El[i]] << ' ' ;
    if(Er[i] < 0)
      out << Er[i] << endl ;
    else
      out << face_indx[Er[i]] << endl ;
  }

  // finally output the values associated with nodes
  // for the purpose of projection visualization, we
  // just put random values here
  // we first init the random seed using the time
  timeval t ;
  gettimeofday(&t,NULL) ;
  srand48(t.tv_usec) ;
  for(int i=0;i<node_num;++i)
    out << drand48() << endl ;

  // close the output file
  out.close() ;
}

// output both projections in 2dgv format
void TwodgvOutputBOTH(const char* fname,
                      const store<vec2d>& npos1,
                      const store<vec2d>& npos2,
                      const dMap& bc1_N1, const dMap& bc1_N2,
                      const dMap& bc1_El, const dMap& bc1_Er,
                      const dMap& bc2_N1, const dMap& bc2_N2,
                      const dMap& bc2_El, const dMap& bc2_Er,
                      int bc1_edge_num,
                      int bc2_edge_num,
                      const entitySet& bc1_faces,
                      const entitySet& bc2_faces) {
  std::ofstream out(fname,ios::out) ;
  out.precision(16) ;
  out << "general" << endl ;

  entitySet npos1_dom = npos1.domain() ;
  entitySet npos2_dom = npos2.domain() ;
  int node_num = npos1_dom.size() + npos2_dom.size() ;
  out << node_num << ' ' << '1' << endl ;

  for(entitySet::const_iterator ei=npos1_dom.begin();
      ei!=npos1_dom.end();++ei)
    out << npos1[*ei].x << ' ' << npos1[*ei].y << endl ;
  for(entitySet::const_iterator ei=npos2_dom.begin();
      ei!=npos2_dom.end();++ei)
    out << npos2[*ei].x << ' ' << npos2[*ei].y << endl ;

  // then build a node and face index store
  store<int> node1_indx ;
  node1_indx.allocate(npos1_dom) ;
  int cnt = 1 ;
  for(entitySet::const_iterator ei=npos1_dom.begin();
      ei!=npos1_dom.end();++ei,++cnt)
    node1_indx[*ei] = cnt ;

  store<int> node2_indx ;
  node2_indx.allocate(npos2_dom) ;
  for(entitySet::const_iterator ei=npos2_dom.begin();
      ei!=npos2_dom.end();++ei,++cnt)
    node2_indx[*ei] = cnt ;
  
  store<int> face1_indx ;
  face1_indx.allocate(bc1_faces) ;
  cnt = 1 ;
  for(entitySet::const_iterator ei=bc1_faces.begin();
      ei!=bc1_faces.end();++ei,++cnt)
    face1_indx[*ei] = cnt ;

  store<int> face2_indx ;
  face2_indx.allocate(bc2_faces) ;
  for(entitySet::const_iterator ei=bc2_faces.begin();
      ei!=bc2_faces.end();++ei,++cnt)
    face2_indx[*ei] = cnt ;

  // then output edges
  int edge_num = bc1_edge_num + bc2_edge_num ;
  out << edge_num << " 1 "
      << bc1_faces.size()+bc2_faces.size()
      << " 1" << endl ;
  for(int i=0;i<bc1_edge_num;++i) {
    out << node1_indx[bc1_N1[i]] << ' '
        << node1_indx[bc1_N2[i]] << ' '
        << face1_indx[bc1_El[i]] << ' ' ;
    if(bc1_Er[i] < 0)
      out << bc1_Er[i] << endl ;
    else
      out << face1_indx[bc1_Er[i]] << endl ;
  }
  for(int i=0;i<bc2_edge_num;++i) {
    out << node2_indx[bc2_N1[i]] << ' '
        << node2_indx[bc2_N2[i]] << ' '
        << face2_indx[bc2_El[i]] << ' ' ;
    if(bc2_Er[i] < 0)
      out << bc2_Er[i] << endl ;
    else
      out << face2_indx[bc2_Er[i]] << endl ;
  }

  // finally output the values associated with nodes
  // for the purpose of projection visualization, we
  // just put random values here
  // we first init the random seed using the time
  timeval t ;
  gettimeofday(&t,NULL) ;
  srand48(t.tv_usec) ;
  for(int i=0;i<node_num;++i)
    out << drand48() << endl ;

  // close the output file
  out.close() ;
}

/////////////////////////////////////////////////////////////
/////            point & polygon algorithms             /////
/////////////////////////////////////////////////////////////
// first this is an implementation of a rectangle(with double
// precision) and various utilities for rectangle and points
bool point_in_poly(const vec2d&,const vector<int>&,
                   const store<vec2d>&) ;
class rectangle {
public:
  // l is the upper left corner,
  // r is the lower right corner.
  vec2d l, r ; // it is defined by its left & right corner
  rectangle(const vec2d& lc=vec2d(0.0,0.0),
            const vec2d& rc=vec2d(0.0,0.0)):l(lc),r(rc) {}
  // determin if the rectangle is valid or not
  // i.e., the l corner is to the left and top
  // of the right corner. the degenerated case
  // is a point only
  bool is_valid() const
  {return ( (l.x<=r.x) && (l.y>=r.y)) ;}
  // check valid and stop if not valid
  void check_valid() const {
    if(!is_valid()) {
      cerr << "ERROR: rectangle definition invalid!" << endl ;
      cerr << "left corner: " << l << endl ;
      cerr << "right corner: " << r << endl ;
      Loci::Abort() ;  
    }
  }
  // determine if the given point is inside.
  // only points on upper and left edges are
  // counted in the rectangle. points on the
  // bottom and right edges are not inside
  bool is_point_in(const vec2d& p) const ;
  // see comment below above the implementation
  int point_location(const vec2d& p) const ;
  // determine if the given segment intersect with or lie in
  // the rectangle, this method only detects boolean
  // intersection, it does not really computes the
  // intersection points
  bool is_seg_int(const vec2d& b, const vec2d& e) const ;
  // determine if the given rectangle intersects
  // and if so, computes the intersected rectangle
  // this function does not count the degenerated
  // cases as intersection (i.e., the intersection
  // is an edge or a point)
  bool rectangle_int(const rectangle& rc,
                     rectangle& new_rc) const ;
  // the same as the above one, except this one
  // does not compute the intersection rectangle
  bool rectangle_int(const rectangle& rc) const ;
  // determin if a circle intersects the rectangle.
  // if the circle barely touches the bottom and the right
  // edges of the rectangle, it is considered NOT intersected!
  bool is_circle_int(const vec2d& center, double radius) const ;
  // determin if a convex polygon intersects the rectangle
  // the polygon is defined by a vector of its nodes x and y value
  bool is_convex_polygon_int(const vector<int>& poly_nodes,
                             const store<vec2d>& node_pos) const ;
} ;
inline bool rectangle::is_point_in(const vec2d& p) const {
  // first we check if the rectangle is valid
  check_valid() ;
  // only points on upper and left edges are
  // counted in the rectangle. points on the
  // bottom and right edges are not inside
  return ( (p.x >= l.x) && (p.x < r.x) &&
           (p.y <= l.y) && (p.y > r.y)) ;
}
inline bool rectangle::rectangle_int(const rectangle& rc,
                                     rectangle& new_rc) const {
  // first we check if the rectangle is valid
  check_valid() ;
  rc.check_valid() ;
  // first we determin if they intersect
  if( (r.x <= rc.l.x) || (l.x >= rc.r.x) ||
      (r.y >= rc.l.y) || (l.y <= rc.r.y))
    return false ;
  // then we compute the new intersected rectangle
  // since now the two rectangles are surely intersected
  // the point positions should satisfy the following conditions
  if(l.x >= rc.l.x)
    new_rc.l.x = l.x ;
  else
    new_rc.l.x = rc.l.x ;

  if(l.y <= rc.l.y)
    new_rc.l.y = l.y ;
  else
    new_rc.l.y = rc.l.y ;

  if(r.x <= rc.r.x)
    new_rc.r.x = r.x ;
  else
    new_rc.r.x = rc.r.x ;

  if(r.y >= rc.r.y)
    new_rc.r.y = r.y ;
  else
    new_rc.r.y = rc.r.y ;

  return true ;
}
inline bool rectangle::rectangle_int(const rectangle& rc) const {
  // first we check if the rectangle is valid
  check_valid() ;
  rc.check_valid() ;
  // first we determin if they intersect
  if( (r.x <= rc.l.x) || (l.x >= rc.r.x) ||
      (r.y >= rc.l.y) || (l.y <= rc.r.y))
    return false ;

  return true ;
}
inline bool rectangle::is_circle_int(const vec2d& center,
                                     double radius) const {
  check_valid() ;
  // if the center of the circle is in the rectangle
  // then the circle of course intersects with the rectangle
  if(is_point_in(center))
    return true ;
  // test if any of the four rectangle edges intersects the circle
  if(fabs(center.x - l.x) <= radius)
    return true ;
  if(fabs(center.x - r.x) <= radius)
    return true ;
  if(fabs(center.y - l.y) <= radius)
    return true ;
  if(fabs(center.y - r.y) <= radius)
    return true ;
  return false ;
}
// this utility function returns the location of a point in
// the plane divided by the rectangle.
//                0 | 1 | 2
//                ---------
//                3 | 4 | 5
//                ---------
//                6 | 7 | 8
// region 4 is the rectangle. NOTE the inequality relations
// in the function are set in careful.
inline int rectangle::point_location(const vec2d& p) const {
  if(is_point_in(p))
    return 4 ;
  else if( (p.x < l.x) && (p.y > l.y))
    return 0 ;
  else if( ((p.x >= l.x) && (p.x < r.x)) && (p.y > l.y))
    return 1 ;
  else if( (p.x >= r.x) && (p.y > l.y))
    return 2 ;
  else if( (p.x < l.x) && ((p.y <= l.y)&&(p.y > r.y)))
    return 3 ;
  else if( (p.x >= r.x) && ((p.y <= l.y)&&(p.y > r.y)))
    return 5 ;
  else if( (p.x < l.x) && (p.y <= r.y))
    return 6 ;
  else if( ((p.x >= l.x) && (p.x < r.x)) && (p.y <= r.y))
    return 7 ;
  else
    return 8 ;
}
bool rectangle::is_seg_int(const vec2d& b, const vec2d& e) const {
  // quick test first
  if(is_point_in(b) || is_point_in(e))
    return true ;
  double min_x, max_x, min_y, max_y ;
  min_x = min(b.x,e.x) ; max_x = max(b.x,e.x) ;
  min_y = min(b.y,e.y) ; max_y = max(b.y,e.y) ;
  if( (max_x <= l.x) || (min_x >= r.x) ||
      (min_y >= l.y) || (max_y <= r.y))
    return false ;
  bool left = false, right = false ;
  left_exact(b,e,l)?left=true:right=true ;
  left_exact(b,e,r)?left=true:right=true ;
  left_exact(b,e,vec2d(l.x,r.y))?left=true:right=true ;
  left_exact(b,e,vec2d(r.x,l.y))?left=true:right=true ;
  if( (!left && right) || (left && !right))
    return false ;
  else
    return true ;
}
bool
rectangle::is_convex_polygon_int(const vector<int>& poly_nodes,
                                 const store<vec2d>& node_pos) const {
  check_valid() ;
  // we test if the four vertex of the rectangle falls
  // in the polygon, if any, then they are intersected
  if(point_in_poly(l,poly_nodes,node_pos))
    return true ;
  if(point_in_poly(r,poly_nodes,node_pos))
    return true ;
  if(point_in_poly(vec2d(l.x,r.y),poly_nodes,node_pos))
    return true ;
  if(point_in_poly(vec2d(r.x,l.y),poly_nodes,node_pos))
    return true ;

  // then we test if any edge intersects with the
  // rectangle
  vector<int>::size_type nd = poly_nodes.size() ;
  for(vector<int>::size_type i=0;i!=nd;++i) {
    // get the nodes sequentially
    const vec2d& n1 = node_pos[poly_nodes[(i==0)?(nd-1):(i-1)]] ;
    const vec2d& n2 = node_pos[poly_nodes[i]] ;
    if(is_seg_int(n1,n2))
      return true ;
  }
  return false ;
}
// define the output
std::ostream & operator<<(std::ostream& s, const rectangle& r)
{
  s << "(<" << r.l.x << ", " << r.l.y << ">" << " | "
    << "<" << r.r.x << ", " << r.r.y << ">)" ;
  return s ;
}
//////////////end of rectangle implementation///////////////////

/*********************************************************
 * Functions that compute Point and Polygon Relations    *
 ********************************************************/
// given a point and all the nodes that define the polygon
// return if the point is on any edge of the polygon
bool point_on_poly_edge(const vec2d& p,
                        const vector<double>& x,
                        const vector<double>& y) {
  typedef vector<double>::size_type size ;
  size xs = x.size() ;
  size ys = y.size() ;
  if(xs != ys) {
    cerr << "ERROR in point_in_poly, polygon ill defined!" << endl ;
    Loci::Abort() ;
  }
  size i, j ;
  for(i=0,j=xs-1;i<xs;j=i++) {
    vec2d p1 = vec2d(x[j],y[j]) ;
    vec2d p2 = vec2d(x[i],y[i]) ;
    // here we are using the exact version of the
    // point on edge test.
    // this function is solely used in the building
    // of the quadtree, so we wish to compute the
    // relations exactly.
    if(pointOnEdge_exact(p,p1,p2))
      return true ;
  }
  return false ;
}

// given a point and a polygon, return if the point is inside
// the polygon (from a well known algorithm that derives from
// the Jordan curve theorem) The code is actually from
// Randolph Franklin at RPI
bool point_in_poly(const vec2d& p,
                   const vector<double>& x,
                   const vector<double>& y) {
  typedef vector<double>::size_type size ;
  size xs = x.size() ;
  size ys = y.size() ;
  if(xs != ys) {
    cerr << "ERROR in point_in_poly, polygon ill defined!" << endl ;
    Loci::Abort() ;
  }
  
  size i, j ;
  bool in = false ;
  for(i=0,j=xs-1;i<xs;j=i++) {
    if ((((y[i] <= p.y) && (p.y < y[j])) ||
         ((y[j] <= p.y) && (p.y < y[i]))) &&
        (p.x < (x[j] - x[i]) * (p.y - y[i]) / (y[j] - y[i]) + x[i]))
      in = !in ;
  }
  return in ;
}

// a wrapper of point_in_poly, the polygon is defined
// by a nodeSet and a store of position
// but at here, we first check to see if the point locates
// on any vertex of the polygon, if it is, we don't consider
// them inside the polygon in this program and we would
// return false immediately

// NOTE: this is also the first version of point_in_poly,
// in this function, if a point is on a vertex of the polygon,
// we do NOT regart it in. But if the point is on an edge
// of the polygon, we consider it inside.
bool point_in_poly(const vec2d& p,
                   const vector<int>& poly_nodes,
                   const store<vec2d>& node_pos) {
  vector<double> poly_x, poly_y ;
  for(vector<int>::const_iterator ei=poly_nodes.begin();
      ei!=poly_nodes.end();++ei) {
    const vec2d& ppoly = node_pos[*ei] ;
    if(p==ppoly) {
      return false ;
    }
    const double& x=ppoly.x ;
    const double& y=ppoly.y ;
    poly_x.push_back(x) ;
    poly_y.push_back(y) ;
  }
  // the we see if the point is on any edge
  // if so, we regarded it in the polygon
  if(point_on_poly_edge(p,poly_x,poly_y)) {
    return true ;
  }
  return point_in_poly(p,poly_x,poly_y) ;
}
// a different version of point_in_poly, points locate
// on vertex and edges are also considered inside
// and polygon is defined by a vector of vec2d
bool point_in_poly2(const vec2d& p,
                    const vector<vec2d>& poly_nodes) {
  vector<double> poly_x, poly_y ;
  for(vector<vec2d>::const_iterator ei=poly_nodes.begin();
      ei!=poly_nodes.end();++ei) {
    if(p == *ei) {
      return true ;
    }
    const double& x=ei->x ;
    const double& y=ei->y ;
    poly_x.push_back(x) ;
    poly_y.push_back(y) ;
  }
  // the we see if the point is on any edge
  // if so, we regarded it in the polygon
  if(point_on_poly_edge(p,poly_x,poly_y)) {
    return true ;
  }
  return point_in_poly(p,poly_x,poly_y) ;
}
// here is another version for point_in_poly2.
// points located on vertices and edges are also
// inside the polygon. But the polygon is defined
// by a node index and a store of vec2d
bool point_in_poly2b(const vec2d& p,
                     const vector<int>& poly_nodes,
                     const store<vec2d>& node_pos) {
  vector<double> poly_x, poly_y ;
  for(vector<int>::const_iterator ei=poly_nodes.begin();
      ei!=poly_nodes.end();++ei) {
    const vec2d& ppoly = node_pos[*ei] ;
    if(p==ppoly) {
      return true ;
    }
    const double& x=ppoly.x ;
    const double& y=ppoly.y ;
    poly_x.push_back(x) ;
    poly_y.push_back(y) ;
  }
  // the we see if the point is on any edge
  // if so, we regarded it in the polygon
  if(point_on_poly_edge(p,poly_x,poly_y)) {
    return true ;
  }
  return point_in_poly(p,poly_x,poly_y) ;
}
// Yet another version of point_in_poly. This time,
// point on vertex and edges are NOT considered inside
bool point_in_poly3(const vec2d& p,
                    const vector<int>& poly_nodes,
                    const store<vec2d>& node_pos) {
  vector<double> poly_x, poly_y ;
  for(vector<int>::const_iterator ei=poly_nodes.begin();
      ei!=poly_nodes.end();++ei) {
    const vec2d& ppoly = node_pos[*ei] ;
    if(p == ppoly) {
      return false ;
    }
    const double& x=ppoly.x ;
    const double& y=ppoly.y ;
    poly_x.push_back(x) ;
    poly_y.push_back(y) ;
  }
  // the we see if the point is on any edge
  // if so, we regarded it in the polygon
  if(point_on_poly_edge(p,poly_x,poly_y)) {
    return false ;
  }
  return point_in_poly(p,poly_x,poly_y) ;
}

/***************************************************
 *   2D polygon area computation routines          *
 **************************************************/
// this function gives the area of a general polygon
// uses a signed area algorithm from Joseph O'Rourke's
// "Computational Geometry in C (2nd ed.)"
double polygon_area(const vector<vec2d>& P) {
  double area = 0.0 ;
  vector<vec2d>::size_type p_size = P.size() ;
  for(vector<vec2d>::size_type i=0;i!=p_size;++i) {
    vector<vec2d>::size_type di = i ;
    vector<vec2d>::size_type dip1 = (i==p_size-1)?(0):(i+1) ;
    area += (P[di].x + P[dip1].x) * (P[dip1].y - P[di].y) ;
  }
  area /= 2.0 ;
  if(area < 0.0) {
    area = -area ;
  }
  
  return area ;
}
// an overloaded version for convenient purpose
double polygon_area(const vector<Entity>& nodes,
                    const store<vec2d>& pos) {
  vector<vec2d> poly ;
  for(vector<int>::const_iterator vi=nodes.begin();
      vi!=nodes.end();++vi)
    poly.push_back(pos[*vi]) ;

  return (polygon_area(poly)) ;
}

/**************************************************
 *    Bounding box computations for a point set   *
 *************************************************/
// given a set of nodes position, return the bounding box
void nodes_bbox(const store<vec2d>& pos,
                rectangle& b) {
  vector<double> xcoord, ycoord ;
  entitySet dom = pos.domain() ;
  if(dom == EMPTY) {
    b.l = vec2d(0.0,0.0) ;
    b.r = vec2d(0.0,0.0) ;
    return ;
  }
  entitySet::const_iterator ei = dom.begin() ;
  
  double min_x = pos[*ei].x ;
  double max_x = min_x ;
  double min_y = pos[*ei].y ;
  double max_y = min_y ;

  for(++ei;ei!=dom.end();++ei) {
    double x = pos[*ei].x ;
    double y = pos[*ei].y ;
    if(x<min_x) min_x = x ;
    if(x>max_x) max_x = x ;
    if(y<min_y) min_y = y ;
    if(y>max_y) max_y = y ;
  }
  b.l = vec2d(min_x,max_y) ;
  b.r = vec2d(max_x,min_y) ;
}
// another overloaded form
void nodes_bbox(const vector<vec2d>& pos,
                rectangle& b) {
  vector<double> xcoord, ycoord ;
  if(pos.empty()) {
    b.l = vec2d(0.0,0.0) ;
    b.r = vec2d(0.0,0.0) ;
    return ;
  }
  vector<vec2d>::const_iterator ei=pos.begin() ;
  
  double min_x = ei->x ;
  double max_x = min_x ;
  double min_y = ei->y ;
  double max_y = min_y ;
  
  for(++ei;ei!=pos.end();++ei) {
    double x = ei->x ;
    double y = ei->y ;
    if(x<min_x) min_x = x ;
    if(x>max_x) max_x = x ;
    if(y<min_y) min_y = y ;
    if(y>max_y) max_y = y ;
  }

  b.l = vec2d(min_x,max_y) ;
  b.r = vec2d(max_x,min_y) ;
}
// another overloaded form
void nodes_bbox(const vector<int>& nodes,
                const store<vec2d>& pos,
                rectangle& b) {
  if(nodes.empty())
    return ;

  double min_x = pos[nodes[0]].x ;
  double max_x = min_x ;
  double min_y = pos[nodes[0]].y ;
  double max_y = min_y ;

  for(size_t i=1;i!=nodes.size();++i) {
    double x = pos[nodes[i]].x ;
    double y = pos[nodes[i]].y ;
    if(x<min_x) min_x = x ;
    if(x>max_x) max_x = x ;
    if(y<min_y) min_y = y ;
    if(y>max_y) max_y = y ;
  }

  b.l = vec2d(min_x,max_y) ;
  b.r = vec2d(max_x,min_y) ;
}

/**
*************************************************************
*   General convex polygon - polygon intersection routine   *
*   Based on the linear algorithm by Joseph O'Rourke et al  *
*   In Journal of Computer Graphics and Image Processing    *
*   19, 384-391 (1982) (the poly_and_poly function)         *
*************************************************************
**/
// returns if the vector contains duplicated entities
inline bool detect_dup(vector<Entity> v) {
  sort(v.begin(),v.end()) ;
  vector<Entity>::iterator new_end = unique(v.begin(),v.end()) ;
  return (new_end != v.end()) ;
}
// this function computes if the passed in vector
// of points are collinear
// it calls the table lookup version of the collinear3
// function, so we'd have to pass in the ids for
// the points also.
// the length of p and pid should be equal.
// pid[i] should be the id for point in p[i]
bool collinear(const vector<vec2d>& p,
               const vector<Entity>& pid) {
  // if p just contains zero, one, or two points,
  // then they are collinear
  vector<vec2d>::size_type s = p.size() ;
  if(s<=2)
    return true ;
  // if all points in p are collinear, p is collinear
  // if any three points are NOT collinear, then p
  // is NOT collinear
  vector<vec2d>::size_type i,ic,ic1,ic2 ;
  for(i=0;i<s;++i) {
    ic = i ; ic1 = (i+1)%s ; ic2 = (i+2)%s ;
    if(!collinear3_tl(p[ic],p[ic1],p[ic2],pid[ic],pid[ic1],pid[ic2]))
      return false ;
  }
  return true ;
}

/***************************************************************
 * Type definition for two segments position                   *
 * PARALLEL   :  two segments are parallel but not collinear   *
 * NOINT      :  two segments are NOT intersected, nor do      *
 *                   they parallel or collinear to each other  *
 * NORMALINT  :  two segments have a normal intesection point  *
 * COLLINEAR  :  two segments collinear                        *
 **************************************************************/
typedef enum {PARALLEL,NOINT,NORMALINT,COLLINEAR} seg_int_type ;
/**
*********************************************************************
*   Implementation of the global segment segment intersection
*   lookup table
**/
// a small object used in recording segment intersection result
struct ssi_result {
  seg_int_type intersect_type ; // intersection result
  vec2d intersect_pos ; // if intersected, the point
  Entity intersect_id ; // and the point's id
  ssi_result() {intersect_id = -1 ;}
  ssi_result(seg_int_type t, const vec2d& p, Entity id):
    intersect_type(t),intersect_pos(p),intersect_id(id) {}
} ;
class seg_int_lookup_table {
  // levels of lookup table
  struct final {
    entitySet domain ;
    map<Entity,ssi_result> record ;
    bool inSet(Entity id) const {return domain.inSet(id) ;}
    void clear() {domain=EMPTY ; record.clear() ;}
    void add_record(Entity l2e, const ssi_result& r) {
      domain += l2e ;
      record[l2e] = r ;
    }
    bool get_record(Entity l2e, ssi_result& r) const {
      if(!inSet(l2e)) // no record yet
        return false ;
      else {
        r = record.find(l2e)->second ;
      }
      return true ;
    }
    ssi_result& operator()(Entity l2e) {
      if(!inSet(l2e))
        domain += l2e ;
      return record[l2e] ;
    }
  } ;
  
  struct third {
    entitySet domain ;
    map<Entity,final> record ;
    bool inSet(Entity id) const {return domain.inSet(id) ;}
    void clear() {domain=EMPTY ; record.clear() ;}
    void add_record(Entity l2s, Entity l2e, const ssi_result& r) {
      if(inSet(l2s))
        record[l2s].add_record(l2e,r) ;
      else {
        final fin ; fin.add_record(l2e,r) ;
        domain += l2s ;
        record[l2s] = fin ;
      }
    }
    bool get_record(Entity l2s, Entity l2e, ssi_result& r) const {
      if(inSet(l2s)) 
        return (record.find(l2s)->second).get_record(l2e,r) ;
      else
        return false ;
    }
    ssi_result& operator()(Entity l2s, Entity l2e) {
      if(!inSet(l2s))
        domain += l2s ;
      return record[l2s](l2e) ;
    }
  } ;
  
  struct second {
    entitySet domain ;
    map<Entity,third> record ;
    bool inSet(Entity id) const {return domain.inSet(id) ;}
    void clear() {domain=EMPTY ; record.clear() ;}
    void add_record(Entity l1e, Entity l2s, Entity l2e,
                    const ssi_result& r) {
      if(inSet(l1e))
        record[l1e].add_record(l2s,l2e,r) ;
      else { // not in, create a new one
        third thi ; thi.add_record(l2s,l2e,r) ;
        domain += l1e ;
        record[l1e] = thi ;
      }
    }
    bool get_record(Entity l1e, Entity l2s, Entity l2e,
                    ssi_result& r) const {
      if(inSet(l1e))
        return (record.find(l1e)->second).get_record(l2s,l2e,r) ;
      else
        return false ;
    }
    bool find_record(Entity le) const {
      return inSet(le) ;
    }
    ssi_result& operator()(Entity l1e, Entity l2s, Entity l2e) {
      if(!inSet(l1e))
        domain += l1e ;
      return record[l1e](l2s,l2e) ;
    }
  } ;

  struct first {
    entitySet domain ;
    map<Entity,second> record ;
    bool inSet(Entity id) const {return domain.inSet(id) ;}
    void clear() {domain=EMPTY ; record.clear() ;}
    void add_record(Entity l1s, Entity l1e, Entity l2s, Entity l2e,
                    const ssi_result& r) {
      if(inSet(l1s)) {
        record[l1s].add_record(l1e,l2s,l2e,r) ;
      }else {
        // completely not in the record, create a whole new record
        second sec ; sec.add_record(l1e,l2s,l2e,r) ;
        domain += l1s ;
        record[l1s] = sec ;
      }
    }
    bool get_record(Entity l1s, Entity l1e, Entity l2s, Entity l2e,
                    ssi_result& r) const {
      if(inSet(l1s))
        return (record.find(l1s)->second).get_record(l1e,l2s,l2e,r) ;
      else
        return false ;
    }
    bool find_record(Entity ls, Entity le) const {
      if(inSet(ls))
        return (record.find(ls)->second).find_record(le) ;
      else
        return false ;
    }
    ssi_result& operator()(Entity l1s, Entity l1e, Entity l2s, Entity l2e) {
      if(!inSet(l1s))
        domain += l1s ;
      return record[l1s](l1e,l2s,l2e) ;
    }
  } ;

  // the table data
  first table ;
  // the table stores segment segment intersection information
  // in a particular order so that later testing can reuse the information.
  // basically, the two segments (four node ids) are sorted so that
  // each segment's head node is the node with smaller id number,
  // and also the segment with a smaller head node comes first.

  // this function takes two segments definition and produces the
  // reordered segments used for accessing the table
  void sort_segs(Entity l1s, Entity l1e, Entity l2s, Entity l2e,
                 Entity& new_l1s, Entity& new_l1e,
                 Entity& new_l2s, Entity& new_l2e) const {
    new_l1s = min(l1s,l1e) ; new_l1e = max(l1s,l1e) ;
    new_l2s = min(l2s,l2e) ; new_l2e = max(l2s,l2e) ;
    // NOTE: we cannot do this step if two boundaries are
    // split independently.
//     if(new_l1s > new_l2s) {
//       std::swap(new_l1s, new_l2s) ;
//       std::swap(new_l1e, new_l2e) ;
//     }
  }

  // The following are public interfaces
public:
  void clear() {
    table.clear() ;
  }
  void add_record(Entity l1s, Entity l1e, Entity l2s, Entity l2e,
                  const ssi_result& record) {
    // first sort each segment's id
    Entity new_l1s,new_l1e,new_l2s,new_l2e ;
    sort_segs(l1s,l1e,l2s,l2e,new_l1s,new_l1e,new_l2s,new_l2e) ;
    
    table.add_record(new_l1s,new_l1e,new_l2s,new_l2e,record) ;
  }
  bool find_record(Entity l1s, Entity l1e, Entity l2s, Entity l2e,
                   ssi_result& record) const {
    // first sort each segment's id
    Entity new_l1s,new_l1e,new_l2s,new_l2e ;
    sort_segs(l1s,l1e,l2s,l2e,new_l1s,new_l1e,new_l2s,new_l2e) ;
    
    return table.get_record(new_l1s,new_l1e,new_l2s,new_l2e,record) ;
  }
  // this one merely query if there are any record
  // associated with line segments (ls,le)
  bool find_record(Entity ls, Entity le) const {
    // first sort each segment's id
    Entity new_ls,new_le ;
    new_ls = min(ls,le) ; new_le = max(ls,le) ;
    
    return table.find_record(new_ls,new_le) ;
  }
  ssi_result& operator()(Entity l1s, Entity l1e, Entity l2s, Entity l2e) {
    // first sort each segment's id
    Entity new_l1s,new_l1e,new_l2s,new_l2e ;
    sort_segs(l1s,l1e,l2s,l2e,new_l1s,new_l1e,new_l2s,new_l2e) ;

    return table(new_l1s,new_l1e,new_l2s,new_l2e) ;
  }
} ;
// a global lookup table created for recording
// segment segment intersection results
seg_int_lookup_table global_ssint_ltable ;
/////////////////end of global seg int lookup table///////////////////

// This is the current new node index, it will be initialized
// in the face_split function and incremented in the seg_seg_int
// function each time a new intersection point is found.
// IMPORTANT:
// A NOTE to the way we generate mesh_new_node_index:
// When we split faces on a give boundary(master boundary), any node
// that is not
// from the boundary will be treated new node and hence they require
// a new index number. These new nodes come from two ways: the first,
// any newly generated intersections (I mean proper intersections,
// i.e., intersections that are in the interior of both segments)
// are clearly new nodes; second, any node from the other
// boundary (slave boundary)
// which forms a new split face is also new node to the master boundary.
// These nodes also require a new index number.
// The first place we generate new index number is in the segment
// intersection routine (seg_seg_int). First, the passed in arguments are
// segment 1 (refereed to as s1) and segment 2 (refereed to as s2).
// The order we pass in s1 and s2 guarantees that s1 is the
// segment on the master boundary and s2 is the segment on the slave
// boundary. Possible new index generation occurs for a proper
// intersection or an improper intersection that involves s2 end points.
// We now examine each cases in seg_seg_int. First, if s1 s2 parallel,
// then there is no intersection and there is no need to generate
// a new node index number. Second, if s1 s2 collinear, then they might
// intersect with each other, or may not. But collinear is a special
// case in our polygon-polygon intersection, we don't actually need
// the intersection points. Thus in this case, we don't need to generate
// a new node index also. And the actual index setting is unimportant.
// Third, s1 s2 intersect at any end points.
// then there is NO new node since the intersection is a node in
// the master boundary (because it is one of s1's end point). Hence,
// no need to generate new node index too in this case. We just
// need to reuse the existing index from s1 end point. Fourth case,
// one end point of s1 is on s2. Then the intersection point is
// NOT new, we don't need to have new node index for it. And we just
// need to reuse the existing index from s1. Fifth
// case, one end point of s2 is on s1. In this case, since the
// intersection point is from s2 which is from the slave boundary,
// then we DO NEED to generate a new node index for the intersection.
// The last case, a proper intersection. As discussed, in this case,
// we will need to generate a new index number for the intersection.
// Now the second place that we could generate a new node index
// is at the polygon intersection routine. In the polygon intersection
// routine, the passed in P polygon is on the master boundary and
// Q polygon in on the slave boundary. When we advance_Q, we may
// output the current point on Q if QIN flag is on. This also counts
// a new point for the master boundary. Therefore, we also need to
// generate a new node for that output point.
static Entity mesh_new_node_index ;
// mesh new face index is relative easy to maintain. It is only
// incremented after each successful face spliting.
static Entity mesh_new_face_index ;

// The following map records the new index for
// every node that is NOT originally on the master boundary
static dMap global_Q_id_remap ;
// These two global variables are used to indicate
// the current nodes in the master and slave boundaries
// that are on the boundary edges.
static entitySet masterbc_edge_nodes ;
static entitySet slavebc_edge_nodes ;

// this function computes the intersection of two lines,
// it returns NOINT if they do not intersect, NORMALINT if they
// intersect and set the intersection point. NOTE: if two
// lines are coincide, then we do not regard them intersect
// and will return a special code COLLINEAR.
// If one line has an end point that touches the other line,
// we consider this a valid intersection.
// NOTE: we used the global segment intersection
// lookup table in computing the intersection
// NOTE: lx_start is the head of segment lx
//       lx_end is the tail of segment lx

seg_int_type
seg_seg_int(const vec2d& l1_start, const vec2d& l1_end,
            const vec2d& l2_start, const vec2d& l2_end,
            Entity l1s, Entity l1e, Entity l2s, Entity l2e,
            vec2d& intersection, Entity& int_id) {
  // first we need to lookup the global table
  ssi_result si ;
  if(global_ssint_ltable.find_record(l1s,l1e,l2s,l2e,si)) {
    intersection = si.intersect_pos ;
    int_id = si.intersect_id ;
    return si.intersect_type ;
  }
  
  // we first use seg_seg_pos to determine the position
  // of the line segments. NOTE: here we are using
  // the table lookup version of seg_seg_pos predicate
  seg_seg_pos_type pos_type =
    seg_seg_pos_tl(l1_start,l1_end,l2_start,l2_end,
                   l1s,l1e,l2s,l2e) ;

  // if two segments are parallel, then they don't have
  // intersection
  if(pos_type == POS_PAR) {
    // add record to global table
    si.intersect_type = PARALLEL ;
    global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;

    return PARALLEL ;
  }
  // if they are collinear, then they may have valid
  // intersection point
  // NOTE: we are using the table lookup version
  // of the between predicate
  if(pos_type == POS_COLLINEAR) {
#ifdef UNSTABLE_SEG_SEG_INT_CODE
    // NOTE: the code in this ifdef is not stable
    // and is only experimental, and probably will be
    // removed in the future.
    
    // first we want to eliminate some special cases,
    // if two end points on two segments are the same,
    // and the other two end points are at the opposite
    // site, then we regard such a case normal intersection.
    // e.g., l1s<--------(l1e,l2s)<----------l2e
    // then we say the intersection point is l1e.

    // first case: l1e---------->(l1s,l2s)<----------l2e
    if(l1_start == l2_start) {
      if(between(l1_end, l2_end, l1_start)) {
        intersection = l1_start ;
        int_id = l1s ;
        // add record
        si.intersect_type = NORMALINT ;
        si.intersect_pos = intersection ;
        si.intersect_id = int_id ;
        
        global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;
        
        return NORMALINT ;
      }
    }

    // second case: l1e--------->(l1s,l2e)---------->l2s
    if(l1_start == l2_end) {
      if(between(l1_end, l2_start, l1_start)) {
        intersection = l1_start ;
        int_id = l1s ;
        // add record
        si.intersect_type = NORMALINT ;
        si.intersect_pos = intersection ;
        si.intersect_id = int_id ;
        
        global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;
        
        return NORMALINT ;
      }
    }

    // third case: l1s<----------(l1e,l2s)<----------l2e
    if(l1_end == l2_start) {
      if(between(l1_start, l2_end, l1_end)) {
        intersection = l1_end ;
        int_id = l1e ;
        // add record
        si.intersect_type = NORMALINT ;
        si.intersect_pos = intersection ;
        si.intersect_id = int_id ;
        
        global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;
        
        return NORMALINT ;
      }
    }

    // last case: l1s<-----------(l1e,l2e)---------->l2s
    if(l1_end == l2_start) {
      if(between(l1_start, l2_start, l1_end)) {
        intersection = l1_end ;
        int_id = l1e ;
        // add record
        si.intersect_type = NORMALINT ;
        si.intersect_pos = intersection ;
        si.intersect_id = int_id ;
        
        global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;
        
        return NORMALINT ;
      }
    }
#endif
    
    // intersections are l2_start & l2_end
    if(between_tl(l1_start,l1_end,l2_start,l1s,l1e,l2s) &&
       between_tl(l1_start,l1_end,l2_end,l1s,l1e,l2e)) {
      intersection = l2_start ;
      int_id = l2s ;
      // add record
      si.intersect_type = COLLINEAR ;
      si.intersect_pos = intersection ;
      si.intersect_id = int_id ;
      
      global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;
    
      return COLLINEAR ;
    }
    // intersections are l1_start & l1_end
    else if(between_tl(l2_start,l2_end,l1_start,l2s,l2e,l1s) &&
            between_tl(l2_start,l2_end,l1_end,l2s,l2e,l1e)) {
      intersection = l1_start ;
      int_id = l1s ;
      // add record
      si.intersect_type = COLLINEAR ;
      si.intersect_pos = intersection ;
      si.intersect_id = int_id ;
      
      global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;

      return COLLINEAR ;
    }
    // intersections are l2_start & l1_end
    else if(between_tl(l1_start,l1_end,l2_start,l1s,l1e,l2s) &&
            between_tl(l2_start,l2_end,l1_end,l2s,l2e,l1e)) {
      intersection = l2_start ;
      int_id = l2s ;
      // add record
      si.intersect_type = COLLINEAR ;
      si.intersect_pos = intersection ;
      si.intersect_id = int_id ;
      
      global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;

      return COLLINEAR ;
    }
    // intersections are l2_start & l1_start
    else if(between_tl(l1_start,l1_end,l2_start,l1s,l1e,l2s) &&
            between_tl(l2_start,l2_end,l1_start,l2s,l2e,l1s)) {
      intersection = l1_start ;
      int_id = l1s ;
      // add record
      si.intersect_type = COLLINEAR ;
      si.intersect_pos = intersection ;
      si.intersect_id = int_id ;
      
      global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;

      return COLLINEAR ;
    }
    // intersections are l2_end & l1_end
    else if(between_tl(l1_start,l1_end,l2_end,l1s,l1e,l2e) &&
            between_tl(l2_start,l2_end,l1_end,l2s,l2e,l1e)) {
      intersection = l1_end ;
      int_id = l1e ;
      // add record
      si.intersect_type = COLLINEAR ;
      si.intersect_pos = intersection ;
      si.intersect_id = int_id ;
      
      global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;

      return COLLINEAR ;
    }
    // intersections are l2_end & l1_start
    else if(between_tl(l1_start,l1_end,l2_end,l1s,l1e,l2e) &&
            between_tl(l2_start,l2_end,l1_start,l2s,l2e,l1s)) {
      intersection = l2_end ;
      int_id = l2e ;
      // add record
      si.intersect_type = COLLINEAR ;
      si.intersect_pos = intersection ;
      si.intersect_id = int_id ;
      
      global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;

      return COLLINEAR ;
    }
    // no intersections
    else {
      // add record
      si.intersect_type = NOINT ;
      global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;

      return NOINT ;
    }
  }

  // end points check
  if( (l1_start == l2_start) ||
      (l1_start == l2_end)) {
    intersection = l1_start ;
    int_id = l1s ;
    // add record
    si.intersect_type = NORMALINT ;
    si.intersect_pos = intersection ;
    si.intersect_id = int_id ;

    global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;
    
    return NORMALINT ;
  }
  if( (l1_end == l2_start) ||
      (l1_end == l2_end)) {
    intersection = l1_end ;
    int_id = l1e ;
    // add record
    si.intersect_type = NORMALINT ;
    si.intersect_pos = intersection ;
    si.intersect_id = int_id ;

    global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;
    
    return NORMALINT ;
  }

  // then we will need to check whether one end point
  // of one segment is on another segment
  // NOTE: we are using the table lookup version
  // of the point on edge predicate
  if(pointOnEdge_tl(l1_start,l2_start,l2_end,l1s,l2s,l2e)) {
    intersection = l1_start ;
    int_id = l1s ;
    // add record
    si.intersect_type = NORMALINT ;
    si.intersect_pos = intersection ;
    si.intersect_id = int_id ;

    global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;
    global_pbetween_ltable.add_record(int_id,l2s,l2e) ;

    return NORMALINT ;    
  }
  
  if(pointOnEdge_tl(l1_end,l2_start,l2_end,l1e,l2s,l2e)) {
    intersection = l1_end ;
    int_id = l1e ;
    // add record
    si.intersect_type = NORMALINT ;
    si.intersect_pos = intersection ;
    si.intersect_id = int_id ;

    global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;
    global_pbetween_ltable.add_record(int_id,l2s,l2e) ;

    return NORMALINT ;    
  }
  
  if(pointOnEdge_tl(l2_start,l1_start,l1_end,l2s,l1s,l1e)) {
    intersection = l2_start ;
    int_id = l2s ;
    // add record
    si.intersect_type = NORMALINT ;
    si.intersect_pos = intersection ;
    si.intersect_id = int_id ;

    global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;
    global_pbetween_ltable.add_record(int_id,l1s,l1e) ;

    return NORMALINT ;    
  }
  
  if(pointOnEdge_tl(l2_end,l1_start,l1_end,l2e,l1s,l1e)) {
    intersection = l2_end ;
    int_id = l2e ;
    // add record
    si.intersect_type = NORMALINT ;
    si.intersect_pos = intersection ;
    si.intersect_id = int_id ;

    global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;
    global_pbetween_ltable.add_record(int_id,l1s,l1e) ;

    return NORMALINT ;    
  }

  // OK by now, we have tested all the improper intersections
  // We need to test is there a proper intersection?
  // We first perform a boolean testing on the possibility
  // of intersection. This is based on the method described
  // in Section 1.5.4. of Joseph O'Rourke's book:
  // Computational Geometry in C (second edition)
  // Here is a quote from the book: 
  // if two segments ab and cd intersect in
  // their interiors, then c and d are split by the line
  // containing ab: c is to one side and d to the other.
  // And likewise, a and b are split by the line containing
  // cd. Both conditions must be true.
  // Here we use the exact predicate to do testing
  if(!(xor_bool(left_exact(l1_start,l1_end,l2_start),
                left_exact(l1_start,l1_end,l2_end)) &&
       xor_bool(left_exact(l2_start,l2_end,l1_start),
                left_exact(l2_start,l2_end,l1_end)))) {
    // no intersection
    si.intersect_type = NOINT ;
    global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;
    //global_ssint_ltable(l2s,l2e,l1s,l1e) = si ;

    return NOINT ;
  }
  // then we compute the intersection
  // and this is done in INexact floating-point computation
  double denominator =
    (l2_end.y - l2_start.y) * (l1_end.x - l1_start.x) -
    (l2_end.x - l2_start.x) * (l1_end.y - l1_start.y) ;

  double numerator1 =
    (l2_end.x - l2_start.x) * (l1_start.y - l2_start.y) -
    (l2_end.y - l2_start.y) * (l1_start.x - l2_start.x) ;
  /*
  double numerator2 =
    (l1_end.x - l1_start.x) * (l1_start.y - l2_start.y) -
    (l1_end.y - l1_start.y) * (l1_start.x - l2_start.x) ;
  */
  double u1 = numerator1 / denominator ;
  //double u2 = numerator2 / denominator ;

  intersection = l1_start + u1 * (l1_end - l1_start) ;
  // Now, we have a true valid NEW intersection point
  
  // Since this is a new node, we need to create a
  // new index for this node. And we need to increment
  // the mesh_new_node_index too.
  int_id = mesh_new_node_index++ ;
  
  // add record
  si.intersect_type = NORMALINT ;
  si.intersect_pos = intersection ;
  si.intersect_id = int_id ;
  
  global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;
  // Since this is a proper intersection, the intersection
  // point p lies on both segments. Therefore, we also
  // need to put the record into the global_pbetween_ltable
  global_pbetween_ltable.add_record(int_id,l1s,l1e) ;
  global_pbetween_ltable.add_record(int_id,l2s,l2e) ;

  return NORMALINT ;
}

// this implements a robust (exact) point in convex polygon
// test. In this version, the polygon vertices are assumed
// to be ordered in counter clockwise order. If the point
// locates on any edge, vertex of the polygon, it is also
// counted as inside.
bool point_in_convex_poly_exact1(const vec2d& p,
                                 const vector<vec2d>& poly) {
  vector<vec2d>::size_type i,j,size=poly.size() ;
  for(i=0,j=size-1;i<size;j=i++) {
    if(!leftOn_exact(poly[i],poly[j],p))
      return false ;
  }
  return true ;
}
// this function checks if every vertex of polygon P lies in
// polygon Q
bool poly_vertex_in(const vector<vec2d>& P,
                    const vector<vec2d>& Q) {
  for(vector<vec2d>::const_iterator vi=P.begin();
      vi!=P.end();++vi)
    if(!point_in_convex_poly_exact1(*vi,Q))
      return false ;
  return true ;
}
// this function determins whether a point "p" lies in the
// half plane (to the left) defined by the vector (p1 - p2)
inline int in_half_plane(const vec2d& p1, const vec2d& p2, const vec2d& p,
                         Entity p1_id, Entity p2_id, Entity p_id) {
  // we first lookup in the global_pbetween_ltable,
  // if the record exists, then p1,p2,p are collinear
  // and we need to return 0
  if(global_pbetween_ltable.has_record(p_id,p1_id,p2_id))
    return 0 ;
  // if no record, then we compute it by using the
  // orient2d predicate directly
  double hc[2] = {p1.x, p1.y} ;
  double tc[2] = {p2.x, p2.y} ;
  double pc[2] = {p.x, p.y} ;
  double d = orient2d(tc,hc,pc) ;
  if(d > 0.0)
    return 1 ;
  else if(d < 0.0)
    return -1 ;
  else // d==0.0
    return 0 ;
}

// typedefs useful in polygon intersection routine
typedef enum {UNKNOWN, PIN, QIN} inside_type ;
// advance function
inline void advance_P(inside_type inflag,
                      vector<vec2d>::size_type P_size,
                      const vec2d& P_cur,
                      Entity P_cur_id,
                      Entity first_int_id,
                      Entity& last_output_id,
                      vector<vec2d>::size_type& bug_P,
                      vector<vec2d>::size_type& P_adv,
                      vector<vec2d>& Ins,
                      vector<Entity>& Ins_index) {
#ifdef DEBUG_POLY_AND_POLY
  /////////////
  cerr << "Advance P" << endl ;
  /////////////
#endif
  if(inflag == PIN) {
    if( (P_cur_id != first_int_id) &&
        (P_cur_id != last_output_id)) {
      Ins.push_back(P_cur) ;
      Ins_index.push_back(P_cur_id) ;
      last_output_id = P_cur_id ;
#ifdef DEBUG_POLY_AND_POLY
      ////////////////////////////
      cerr << "Output P" << endl ;
      cerr << "output index: " << P_cur_id << endl ;
      ////////////////////////////
#endif
    }
  }
  bug_P = (bug_P + 1) % P_size ;
  ++ P_adv ;
#ifdef DEBUG_POLY_AND_POLY
  ///////////////////
  cerr << endl ;
  ///////////////////
#endif
}
inline void advance_Q(inside_type inflag,
                      vector<vec2d>::size_type Q_size,
                      const vec2d& Q_cur,
                      Entity Q_cur_id,
                      Entity first_int_id,
                      Entity& last_output_id,
                      vector<vec2d>::size_type& bug_Q,
                      vector<vec2d>::size_type& Q_adv,
                      vector<vec2d>& Ins,
                      vector<Entity>& Ins_index) {
#ifdef DEBUG_POLY_AND_POLY
  ///////////////////
  cerr << "Advance Q" << endl ;
  ///////////////////
#endif
  if(inflag == QIN) {
    if( (Q_cur_id != first_int_id) &&
        (Q_cur_id != last_output_id)) {
      Ins.push_back(Q_cur) ;
      Ins_index.push_back(Q_cur_id) ;
      last_output_id = Q_cur_id ;
      
#ifdef DEBUG_POLY_AND_POLY
      /////////////////////////////
      cerr << "Output Q" << endl ;
      cerr << "output index: " << Q_cur_id << endl ;
      /////////////////////////////
#endif
    }
  }
  bug_Q = (bug_Q + 1) % Q_size ;
  ++ Q_adv ;
#ifdef DEBUG_POLY_AND_POLY
  ///////////////////
  cerr << endl ;
  ///////////////////
#endif
}

// this function is used to handle the two segments collinear
// (but pointing to the same direction) situation. because an
// edge may be defined by several collinear segments, and we
// want to output all those intermediate points, we will need
// to specially handle this case
void
advance_collinear_adv(vector<vec2d>::size_type poly_size,
                      bool& continue_collinear,
                      bool firstIntersect,
                      Entity first_int_id,
                      Entity& last_output_id,
                      Entity& last_p_id,
                      vec2d& last_p_pos,
                      vector<vec2d>::size_type& bug,
                      vector<vec2d>::size_type& adv,
                      Entity bug_id,
                      const vec2d& bug_pos,
                      vector<vec2d>& Ins,
                      vector<Entity>& Ins_index) {
  // first decide whether to ouput last recorded point
  if(continue_collinear
     // note, this conditition is used to ensure
     // that a point is not output until a normal
     // segment intersection occurs
     && !firstIntersect
     && (last_p_id != last_output_id)
     && (last_p_id != first_int_id)) {
#ifdef DEBUG_POLY_AND_POLY
    ///////////////////
    cerr << "Outputing collinear point: " << last_p_id << endl ;
    ///////////////////
#endif
    // then output
    Ins.push_back(last_p_pos) ;
    Ins_index.push_back(last_p_id) ;
    last_output_id = last_p_id ;
  }
  // record node for next round
  
  // NOTE: this is to prevent the collapse of the
  // two segments at the same node and therefore would
  // output the node twice.
  // Thus in case that the advanced node is already
  // a recorded node, then we reset the collinear continue flag
  // to let it start over again and therefore it would not
  // produce the node in the next round.
  if(continue_collinear) {
    if(last_p_pos == bug_pos) {
      continue_collinear = false ;
#ifdef DEBUG_POLY_AND_POLY
      ///////////////////
      cerr << "Resetting collinear flag" << endl ;
      ///////////////////
#endif
      bug = (bug+1) % poly_size ;
      ++adv ;
      return ;
    }
  }
  
  continue_collinear = true ;
  last_p_id = bug_id ;
  last_p_pos = bug_pos ;
#ifdef DEBUG_POLY_AND_POLY
  ///////////////////
  cerr << "Caching collinear point: " << last_p_id << endl ;
  ///////////////////
#endif
  // modify other pointers accordingly
  bug = (bug+1) % poly_size ;
  ++adv ;
}
//
// "edge" is one of the segment vector
// with the x axis;
// "continue_collinear" indicates whether it is collinear
// last time (to determine whether to output points)
// "collinear_point_pos" indicates the last recorded point
// to be output in this round, it is also set during
// this time of evaluation
// "collinear_point_id" is the id of the last recorded
// collinear point
// "firstIntersect" is a flag indicating whether a normal
// segment intersection has occurred before
// "first_int_id" is the id of the first intersection point.
// "last_output_id" is the id of the last output point
// "bug_P" and "bug_Q" are the advancing pointer of polygon P and Q.
// "P_adv" and "Q_adv" are the advancing counters of polygon P and Q.
// "P_size" and "Q_size" are the total number of points defining
// the polygons P and Q
// "Ins" and "Ins_index" are the intersection points and their ids.
void
advance_collinear(const vec2d& P_cur, const vec2d& Q_cur,
                  Entity P_cur_id, Entity Q_cur_id,
                  const vec2d& edge, bool& continue_collinear,
                  vec2d& collinear_point_pos,
                  Entity& collinear_point_id,
                  bool firstIntersect,
                  Entity first_int_id,
                  Entity& last_output_id,
                  vector<vec2d>::size_type& bug_P,
                  vector<vec2d>::size_type& bug_Q,
                  vector<vec2d>::size_type& P_adv,
                  vector<vec2d>::size_type& Q_adv,
                  vector<vec2d>::size_type P_size,
                  vector<vec2d>::size_type Q_size,
                  vector<vec2d>& Ins,
                  vector<Entity>& Ins_index) {
  // first determine who to advance, we will want to advance
  // the segment who is falling behind the other one. "falling
  // behind" here is defined as if the two segments have a
  // less than PI angle with the x axis, then the one with
  // a smaller x projection is the one that is falling behind.
  // Otherwise, if they have a larger than PI degree angle with
  // the x axis, then the one with a larger x component is the
  // one to advance.
  double direction = dot(edge, vec2d(1,0)) ;
  bool advance_p = false ;
  if(direction>0) {
    // NOTE: the equality in these conditions are set deliberately.
    // the intention here is to ensure in case of tie, we will
    // always advance a segment on polygon P, this is to ensure that
    // in case of overlapped points on collinear segments, it is
    // always the one on P that gets output (this is critical to
    // maintain the correctness of the algorithm).
    if(P_cur.x <= Q_cur.x)
      advance_p = true ; // advance P
    else
      advance_p = false ; // advance Q
  } else if(direction<0) {
    if(P_cur.x >= Q_cur.x)
      advance_p = true ;
    else
      advance_p = false ;
  } else { // direction == 0
    // if direction == 0, then we want to
    // compare the y components of the two points.
    if(edge.y > 0) {// pointing upward, pick a smaller y component
      if(P_cur.y <= Q_cur.y)
        advance_p = true ;
      else
        advance_p = false ;
    } else {// pointing downward, pick a larger y component
      if(P_cur.y >= Q_cur.y)
        advance_p = true ;
      else
        advance_p = false ;
    }
  }
  if(advance_p) {
    advance_collinear_adv(P_size, continue_collinear,
                          firstIntersect, first_int_id,
                          last_output_id, collinear_point_id,
                          collinear_point_pos, bug_P, P_adv,
                          P_cur_id, P_cur, Ins, Ins_index) ;
#ifdef DEBUG_POLY_AND_POLY
    ///////////////////
    cerr << "Collinear advancing P" << endl << endl ;
    ///////////////////
#endif
  } else {
    advance_collinear_adv(Q_size, continue_collinear,
                          firstIntersect, first_int_id,
                          last_output_id, collinear_point_id,
                          collinear_point_pos, bug_Q, Q_adv,
                          Q_cur_id, Q_cur, Ins, Ins_index) ;
#ifdef DEBUG_POLY_AND_POLY
    ///////////////////
    cerr << "Collinear advancing Q" << endl << endl ;
    ///////////////////
#endif
  }
}

// this is the boundary version (uses approximate predicate)
void
advance_collinear_bf(const vec2d& P_cur, const vec2d& Q_cur,
                     Entity P_cur_id, Entity Q_cur_id,
                     const vec2d& edge, bool& continue_collinear,
                     vec2d& collinear_point_pos,
                     Entity& collinear_point_id,
                     bool firstIntersect,
                     Entity first_int_id,
                     Entity& last_output_id,
                     vector<vec2d>::size_type& bug_P,
                     vector<vec2d>::size_type& bug_Q,
                     vector<vec2d>::size_type& P_adv,
                     vector<vec2d>::size_type& Q_adv,
                     vector<vec2d>::size_type P_size,
                     vector<vec2d>::size_type Q_size,
                     vector<vec2d>& Ins,
                     vector<Entity>& Ins_index) {
  double direction = dot(edge, vec2d(1,0)) ;
  bool advance_p = false ;
  if(doubleNumEqual(direction, 0.0, ON_LINE_SEG_DOT_THRESHOLD)) {
    if(edge.y > 0) {// pointing upward, pick a smaller y component
      if(P_cur.y <= Q_cur.y)
        advance_p = true ;
      else
        advance_p = false ;
    } else {// pointing downward, pick a larger y component
      if(P_cur.y >= Q_cur.y)
        advance_p = true ;
      else
        advance_p = false ;
    }
  } else if(direction>0) {
    // NOTE: the equality in these conditions are set deliberately.
    // the intention here is to ensure in case of tie, we will
    // always advance a segment on polygon P, this is to ensure that
    // in case of overlapped points on collinear segments, it is
    // always the one on P that gets output (this is critical to
    // maintain the correctness of the algorithm).
    if(P_cur.x <= Q_cur.x)
      advance_p = true ; // advance P
    else
      advance_p = false ; // advance Q
  } else { // direction<0
    if(P_cur.x >= Q_cur.x)
      advance_p = true ;
    else
      advance_p = false ;
  }
  if(advance_p) {
    advance_collinear_adv(P_size, continue_collinear,
                          firstIntersect, first_int_id,
                          last_output_id, collinear_point_id,
                          collinear_point_pos, bug_P, P_adv,
                          P_cur_id, P_cur, Ins, Ins_index) ;
#ifdef DEBUG_POLY_AND_POLY
    ///////////////////
    cerr << "Collinear advancing P" << endl << endl ;
    ///////////////////
#endif
  } else {
    advance_collinear_adv(Q_size, continue_collinear,
                          firstIntersect, first_int_id,
                          last_output_id, collinear_point_id,
                          collinear_point_pos, bug_Q, Q_adv,
                          Q_cur_id, Q_cur, Ins, Ins_index) ;
#ifdef DEBUG_POLY_AND_POLY
    ///////////////////
    cerr << "Collinear advancing Q" << endl << endl ;
    ///////////////////
#endif
  }
}

// this function computes the intersection of two convex
// polygons. it returns false if they do not intersect,
// true if they intersect and set the intersected polygon
// and the intersected polygon node index
// it is assumed that the polygons are convex and
// their vertices are defined in counter clockwise order
bool poly_and_poly(const vector<vec2d>& P,
                   const vector<vec2d>& Q,
                   const vector<Entity>& Pid,
                   const vector<Entity>& Qid,
                   vector<vec2d>& Ins,
                   vector<Entity>& Ins_index) {
  typedef vector<vec2d>::size_type size ;
  size P_size = P.size() ;
  size Q_size = Q.size() ;
  size P_size2 = P_size*2 ;
  size Q_size2 = Q_size*2 ;
  size bug_P = 0, bug_Q = 0 ;
  size P_adv = 0, Q_adv = 0 ;
  inside_type inflag = UNKNOWN ;
  // record the first intersection point
  vec2d first_int(0.0,0.0) ;
  Entity first_int_id = -1 ;
  Entity last_output_id = -1 ;
  bool firstIntersect = true ;
  // flag indicates whether the intersection in previous iteration
  // is equal to the first intersection
  bool is_ins_last_e_first = false ;
  // data used to control output of collinear points
  bool continue_collinear = false ;
  Entity collinear_point_id = -1 ;
  vec2d collinear_point_pos(0.0, 0.0) ;

  do {
    const vec2d& P_cur = P[bug_P] ;
    const vec2d& Q_cur = Q[bug_Q] ;
    size bug_P_pre = (bug_P - 1 + P_size) % P_size ;
    size bug_Q_pre = (bug_Q - 1 + Q_size) % Q_size ;
    const vec2d& P_pre = P[bug_P_pre] ;
    const vec2d& Q_pre = Q[bug_Q_pre] ;

    Entity P_cur_id = Pid[bug_P] ;
    Entity P_pre_id = Pid[bug_P_pre] ;
    Entity Q_cur_id = Qid[bug_Q] ;
    Entity Q_pre_id = Qid[bug_Q_pre] ;

    vec2d seg_int ;
    Entity seg_int_id ;
    seg_int_type seg_int_result =
      seg_seg_int(P_cur,P_pre,Q_cur,Q_pre,
                  P_cur_id,P_pre_id,Q_cur_id,Q_pre_id,
                  seg_int,seg_int_id) ;

    int P_in_Q_plane = in_half_plane(Q_cur,Q_pre,P_cur,
                                     Q_cur_id,Q_pre_id,P_cur_id) ;
    int Q_in_P_plane = in_half_plane(P_cur,P_pre,Q_cur,
                                     P_cur_id,P_pre_id,Q_cur_id) ;

#ifdef DEBUG_POLY_AND_POLY
    ////////////////////////
    cerr.precision(32) ;
    cerr << "P_cur: " << P_cur << " P_pre: " << P_pre << endl ;
    cerr << "Q_cur: " << Q_cur << " Q_pre: " << Q_pre << endl ;
    cerr << "P_cur_id: " << P_cur_id << " P_pre_id: " << P_pre_id << endl ;
    cerr << "Q_cur_id: " << Q_cur_id << " Q_pre_id: " << Q_pre_id << endl ;
    cerr << "P_in_Q_plane: " << P_in_Q_plane << endl ;
    cerr << "Q_in_P_plane: " << Q_in_P_plane << endl ;

    if(seg_int_result == NOINT)
      cerr << "No intersection!" <<endl ;
    else if(seg_int_result == NORMALINT) {
      cerr << "Normal intersection!" << endl ;
      cerr << "intersection point: " << seg_int << endl ;
      cerr << "intersection id: " << seg_int_id << endl ;
      cerr << "last_output_id: " << last_output_id << endl ;
      cerr << "first_int_id: " << first_int_id << endl ;
      cerr << "is_ins_last_e_first: " << is_ins_last_e_first << endl ;
    }
    else if(seg_int_result == COLLINEAR)
      cerr << "Collinear!" << endl ;
    else
      cerr << "Unknown result!" << endl ;
    ///////////////////////
#endif
    if(seg_int_result == NORMALINT) {
      // check whether we could quit
      if(!is_ins_last_e_first && !firstIntersect) {
        if(seg_int_id == first_int_id)
          break ;
      }
      // if first intersection
      if( (inflag == UNKNOWN) && firstIntersect) {
        P_adv = Q_adv = 0 ;
        firstIntersect = false ;
        first_int = seg_int ;
        first_int_id = seg_int_id ;
      }
      if(seg_int_id == first_int_id)
        is_ins_last_e_first = true ;
      else
        is_ins_last_e_first = false ;

      if(seg_int_id != last_output_id) {
        Ins.push_back(seg_int) ;
        Ins_index.push_back(seg_int_id) ;
        last_output_id = seg_int_id ;
      }

      if(P_in_Q_plane > 0)
        inflag = PIN ;
      else if(Q_in_P_plane > 0)
        inflag = QIN ;
    } else {
      is_ins_last_e_first = false ;
    }
    // advance P or Q
    vec2d P_edge = P_cur - P_pre ;
    vec2d Q_edge = Q_cur - Q_pre ;
    int cross_sign = cross2d_exact(Q_cur,Q_pre,P_cur,P_pre) ;
#ifdef DEBUG_POLY_AND_POLY
    ///////////////////////////
    if(inflag == UNKNOWN)
      cerr << "inflag == UNKNOWN" << endl ;
    else if(inflag == PIN)
      cerr << "inflag == PIN" << endl ;
    else if(inflag == QIN)
      cerr << "inflag == QIN" << endl ;
    else
      cerr << "UNKNOWN inflag" << endl ;

    cerr << "cross_sign: " << cross_sign << endl ;
    cerr << "dot(P_edge,Q_edge): " << dot(P_edge,Q_edge) << endl ;
    ///////////////////////////
#endif
    if( (seg_int_result==COLLINEAR) && (dot(P_edge,Q_edge) < 0)) {
      // two edges are colinear and in oppisite position
      // this means the two polygons only share an edge
      // we therefore could exit the loop now.
      // NOTE: we could ouput the shared edges, but since
      // we do not regard them as valid intersection here,
      // we just return
      Ins.clear() ;
      Ins_index.clear() ;
      return false ;
    } else if( (seg_int_result==PARALLEL) &&
               (P_in_Q_plane < 0) && (Q_in_P_plane < 0)) {
      // two edges are parallel and separated, the two
      // polygons are therefore disjoint
      Ins.clear() ;
      Ins_index.clear() ;
      return false ;
    } /*else if(doubleNumEqual(cross_result, 0.0) &&
        (P_in_Q_plane==0) && (Q_in_P_plane==0)) {*/
    else if( (seg_int_result==COLLINEAR) &&
             (P_in_Q_plane==0) && (Q_in_P_plane==0)) {
      // special case A&B collinear and are the same direction,
      // we call the special handling routine
      advance_collinear(P_cur, Q_cur, P_cur_id, Q_cur_id, P_edge,
                        continue_collinear, collinear_point_pos,
                        collinear_point_id,
                        firstIntersect, first_int_id,
                        last_output_id, bug_P, bug_Q,
                        P_adv, Q_adv, P_size, Q_size, Ins, Ins_index) ;
//       if(inflag == PIN)
//         advance_Q(inflag,Q_size,Q_cur,Q_cur_id,first_int_id,
//                   last_output_id,bug_Q,Q_adv,Ins,Ins_index) ;
//       else
//         advance_P(inflag,P_size,P_cur,P_cur_id,first_int_id,
//                   last_output_id,bug_P,P_adv,Ins,Ins_index) ;
    } else if(cross_sign >= 0) {
      continue_collinear = false ;
      if(P_in_Q_plane > 0)
        advance_Q(inflag,Q_size,Q_cur,Q_cur_id,first_int_id,
                  last_output_id,bug_Q,Q_adv,Ins,Ins_index) ;
      else
        advance_P(inflag,P_size,P_cur,P_cur_id,first_int_id,
                  last_output_id,bug_P,P_adv,Ins,Ins_index) ;
    } else {
      continue_collinear = false ;
      if(Q_in_P_plane > 0)
        advance_P(inflag,P_size,P_cur,P_cur_id,first_int_id,
                  last_output_id,bug_P,P_adv,Ins,Ins_index) ;
      else
        advance_Q(inflag,Q_size,Q_cur,Q_cur_id,first_int_id,
                  last_output_id,bug_Q,Q_adv,Ins,Ins_index) ;
    }
  } while( ((P_adv < P_size) || (Q_adv < Q_size)) &&
           (P_adv < P_size2) && (Q_adv < Q_size2)) ;

#ifdef DEBUG_POLY_AND_POLY
  /////////////////
  if(Ins.size() != Ins_index.size()) {
    cerr << "ERROR: Ins and Ins_index have different size!!!" << endl ;
  }
  cerr << "intersection: " << endl ;
  for(vector<vec2d>::const_iterator vi=Ins.begin();
      vi!=Ins.end();++vi)
    cerr << *vi << endl ;
  cerr << "intersection index: " << endl ;
  for(vector<Entity>::const_iterator vi=Ins_index.begin();
      vi!=Ins_index.end();++vi)
    cerr << *vi << endl ;
  cerr << endl ;
  /////////////////
#endif
  if(!Ins.empty()) {
    if(!collinear(Ins,Ins_index)) {
      return true ;
    }
    // we do not regard an edge or a vertex
    // as a valid intersection
    else {
      Ins.clear() ;
      Ins_index.clear() ;
#ifdef DEBUG_POLY_AND_POLY
      /////////////////
      cerr << "No intersection!" << endl ;
      cerr << endl ;
      /////////////////
#endif
      return false ;
    }
  }

  if(poly_vertex_in(P,Q)) {
    // P in Q completely
    Ins = P ;
    Ins_index = Pid ;
#ifdef DEBUG_POLY_AND_POLY
    /////////////////
    cerr << "Intersection is P!" << endl ;
    cerr << endl ;
    /////////////////
#endif
    return true ;
  } else if(poly_vertex_in(Q,P)) {
    // Q in P completely
    Ins = Q ;
    Ins_index = Qid ;
#ifdef DEBUG_POLY_AND_POLY
    /////////////////
    cerr << "Intersection is Q!" << endl ;
    cerr << endl ;
    /////////////////
#endif
    return true ;
  } else {
    // no intersection
    Ins.clear() ;
    Ins_index.clear() ;
#ifdef DEBUG_POLY_AND_POLY
    /////////////////
    cerr << "No intersection!" << endl ;
    cerr << endl ;
    /////////////////
#endif
    return false ;
  }
}

/**
************************************************************
* Followings are special version for polygon intersection
* at edges (the outer contours) of two boundaries.
**/
// first we need to redo some basic geometric predicates
// that are specialized for points and segments at the
// edges of the boundaries. We combine table lookup with
// approximate version of the geometric predicates.
// All predicates in this section have a _bf as suffix.

inline bool collinear3_bf(const vec2d& a, const vec2d& b, const vec2d& c,
                          Entity aid, Entity bid, Entity cid) {
  // first lookup the table
  if(global_pbetween_ltable.has_record(cid,aid,bid))
    return true ;
  // otherwise call collinear3_approx
  return collinear3_approx(a,b,c) ;
}
inline bool between_bf(const vec2d& a, const vec2d& b, const vec2d& c,
                       Entity aid, Entity bid, Entity cid) {
  // first lookup the table
  if(global_pbetween_ltable.has_record(cid,aid,bid))
    return true ;
  // otherwise call pointOnEdge_approx as an
  // approximation to the between predicate
  if(pointOnEdge_approx(c,a,b)) {
    global_pbetween_ltable.add_record(cid,aid,bid) ;
    return true ;
  }else
    return false ;
}
inline bool pointOnEdge_bf(const vec2d& p, const vec2d& lb, const vec2d& le,
                           Entity pid, Entity lb_id, Entity le_id) {
  if(global_pbetween_ltable.has_record(pid,lb_id,le_id))
    return true ;
  if(pointOnEdge_approx(p,lb,le)) {
    global_pbetween_ltable.add_record(pid,lb_id,le_id) ;
    return true ;
  }else
    return false ;
}
inline seg_seg_pos_type
seg_seg_pos_bf(const vec2d& l1s, const vec2d& l1e,
               const vec2d& l2s, const vec2d& l2e,
               Entity l1s_id, Entity l1e_id,
               Entity l2s_id, Entity l2e_id) {
  bool l1s_collinear = collinear3_bf(l1s,l2s,l2e,l1s_id,l2s_id,l2e_id) ;
  bool l1e_collinear = collinear3_bf(l1e,l2s,l2e,l1e_id,l2s_id,l2e_id) ;
  
  if(l1s_collinear && l1e_collinear) {
    return POS_COLLINEAR ;
  } else if(!l1s_collinear && l1e_collinear) {
    return POS_NOTPAR ;
  } else if(l1s_collinear && !l1e_collinear) {
    return POS_NOTPAR ;
  } else {
    // neither l1s nor l1e collinear with segment (l2s,l2e)
    // therefore there are two possibilities: they are parallel
    // or they are not parallel, thus we compute the cross2d
    // if it is zero, then they are parallel, otherwise, they
    // are not parallel
    if(cross2d_exact(l1s,l1e,l2s,l2e) == 0)
      return POS_PAR ;
    else
      return POS_NOTPAR ;
  }
}
// NOTE: this function is intended for segments whose end points
// are on the mesh boundary, i.e., at least three points of
// (l1s,l1e,l2s,l2e) should be on the boundary edge (the outer contours).
seg_int_type
seg_seg_int_bf(const vec2d& l1_start, const vec2d& l1_end,
               const vec2d& l2_start, const vec2d& l2_end,
               Entity l1s, Entity l1e, Entity l2s, Entity l2e,
               vec2d& intersection, Entity& int_id) {
#ifdef DEBUG_POLY_AND_POLY
  ////////////////////////
  cerr << "This is boundary version of seg_seg_int" << endl ;
  ////////////////////////
#endif
  // first we need to lookup the global table
  ssi_result si ;
  if(global_ssint_ltable.find_record(l1s,l1e,l2s,l2e,si)) {
    intersection = si.intersect_pos ;
    int_id = si.intersect_id ;
    return si.intersect_type ;
  }
  // we first use seg_seg_pos to determine the position
  // of the line segments. NOTE: here we are using
  // the table lookup version of seg_seg_pos predicate
  seg_seg_pos_type pos_type =
    seg_seg_pos_bf(l1_start,l1_end,l2_start,l2_end,
                   l1s,l1e,l2s,l2e) ;

  // if two segments are parallel, then they don't have
  // intersection
  if(pos_type == POS_PAR) {
    // add record to global table
    si.intersect_type = PARALLEL ;
    global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;

    return PARALLEL ;
  }
  // if they are collinear, then they may have valid
  // intersection point
  if(pos_type == POS_COLLINEAR) {
#ifdef UNSTABLE_SEG_SEG_INT_CODE
    // NOTE: the code in this ifdef is not stable
    // and is only experimental, and probably will be
    // removed in the future.
    
    // first we want to eliminate some special cases,
    // if two end points on two segments are the same,
    // and the other two end points are at the opposite
    // site, then we regard such a case normal intersection.
    // e.g., l1s<--------(l1e,l2s)<----------l2e
    // then we say the intersection point is l1e.

    // first case: l1e---------->(l1s,l2s)<----------l2e
    if(l1_start == l2_start) {
      if(pointOnEdge_approx(l1_start, l1_end, l2_end)) {
        intersection = l1_start ;
        int_id = l1s ;
        // add record
        si.intersect_type = NORMALINT ;
        si.intersect_pos = intersection ;
        si.intersect_id = int_id ;
        
        global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;
        
        return NORMALINT ;
      }
    }

    // second case: l1e--------->(l1s,l2e)---------->l2s
    if(l1_start == l2_end) {
      if(pointOnEdge_approx(l1_start, l1_end, l2_start)) {
        intersection = l1_start ;
        int_id = l1s ;
        // add record
        si.intersect_type = NORMALINT ;
        si.intersect_pos = intersection ;
        si.intersect_id = int_id ;
        
        global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;
        
        return NORMALINT ;
      }
    }

    // third case: l1s<----------(l1e,l2s)<----------l2e
    if(l1_end == l2_start) {
      if(pointOnEdge_approx(l1_end, l1_start, l2_end)) {
        intersection = l1_end ;
        int_id = l1e ;
        // add record
        si.intersect_type = NORMALINT ;
        si.intersect_pos = intersection ;
        si.intersect_id = int_id ;
        
        global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;
        
        return NORMALINT ;
      }
    }

    // last case: l1s<-----------(l1e,l2e)---------->l2s
    if(l1_end == l2_start) {
      if(pointOnEdge_approx(l1_end, l1_start, l2_start)) {
        intersection = l1_end ;
        int_id = l1e ;
        // add record
        si.intersect_type = NORMALINT ;
        si.intersect_pos = intersection ;
        si.intersect_id = int_id ;
        
        global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;
        
        return NORMALINT ;
      }
    }
#endif
    
    // intersections are l2_start & l2_end
    if(between_bf(l1_start,l1_end,l2_start,l1s,l1e,l2s) &&
       between_bf(l1_start,l1_end,l2_end,l1s,l1e,l2e)) {
      intersection = l2_start ;
      int_id = l2s ;
      // add record
      si.intersect_type = COLLINEAR ;
      si.intersect_pos = intersection ;
      si.intersect_id = int_id ;
      
      global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;
    
      return COLLINEAR ;
    }
    // intersections are l1_start & l1_end
    else if(between_bf(l2_start,l2_end,l1_start,l2s,l2e,l1s) &&
            between_bf(l2_start,l2_end,l1_end,l2s,l2e,l1e)) {
      intersection = l1_start ;
      int_id = l1s ;
      // add record
      si.intersect_type = COLLINEAR ;
      si.intersect_pos = intersection ;
      si.intersect_id = int_id ;
      
      global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;

      return COLLINEAR ;
    }
    // intersections are l2_start & l1_end
    else if(between_bf(l1_start,l1_end,l2_start,l1s,l1e,l2s) &&
            between_bf(l2_start,l2_end,l1_end,l2s,l2e,l1e)) {
      intersection = l2_start ;
      int_id = l2s ;
      // add record
      si.intersect_type = COLLINEAR ;
      si.intersect_pos = intersection ;
      si.intersect_id = int_id ;
      
      global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;

      return COLLINEAR ;
    }
    // intersections are l2_start & l1_start
    else if(between_bf(l1_start,l1_end,l2_start,l1s,l1e,l2s) &&
            between_bf(l2_start,l2_end,l1_start,l2s,l2e,l1s)) {
      intersection = l1_start ;
      int_id = l1s ;
      // add record
      si.intersect_type = COLLINEAR ;
      si.intersect_pos = intersection ;
      si.intersect_id = int_id ;
      
      global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;

      return COLLINEAR ;
    }
    // intersections are l2_end & l1_end
    else if(between_bf(l1_start,l1_end,l2_end,l1s,l1e,l2e) &&
            between_bf(l2_start,l2_end,l1_end,l2s,l2e,l1e)) {
      intersection = l1_end ;
      int_id = l1e ;
      // add record
      si.intersect_type = COLLINEAR ;
      si.intersect_pos = intersection ;
      si.intersect_id = int_id ;
      
      global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;

      return COLLINEAR ;
    }
    // intersections are l2_end & l1_start
    else if(between_bf(l1_start,l1_end,l2_end,l1s,l1e,l2e) &&
            between_bf(l2_start,l2_end,l1_start,l2s,l2e,l1s)) {
      intersection = l2_end ;
      int_id = l2e ;
      // add record
      si.intersect_type = COLLINEAR ;
      si.intersect_pos = intersection ;
      si.intersect_id = int_id ;
      
      global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;

      return COLLINEAR ;
    }
    // no intersections
    else {
      // add record
      si.intersect_type = NOINT ;
      global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;

      return NOINT ;
    }
  }

  // end points check
  if( (l1_start == l2_start) ||
      (l1_start == l2_end)) {
    intersection = l1_start ;
    int_id = l1s ;
    // add record
    si.intersect_type = NORMALINT ;
    si.intersect_pos = intersection ;
    si.intersect_id = int_id ;

    global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;
    
    return NORMALINT ;
  }
  if( (l1_end == l2_start) ||
      (l1_end == l2_end)) {
    intersection = l1_end ;
    int_id = l1e ;
    // add record
    si.intersect_type = NORMALINT ;
    si.intersect_pos = intersection ;
    si.intersect_id = int_id ;

    global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;
    
    return NORMALINT ;
  }

  // then we will need to check whether one end point
  // of one segment is on another segment
  // NOTE: we are using the table lookup version
  // of the point on edge predicate
  if(pointOnEdge_bf(l1_start,l2_start,l2_end,l1s,l2s,l2e)) {
    intersection = l1_start ;
    int_id = l1s ;
    // add record
    si.intersect_type = NORMALINT ;
    si.intersect_pos = intersection ;
    si.intersect_id = int_id ;

    global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;
    global_pbetween_ltable.add_record(int_id,l2s,l2e) ;

    return NORMALINT ;    
  }
  
  if(pointOnEdge_bf(l1_end,l2_start,l2_end,l1e,l2s,l2e)) {
    intersection = l1_end ;
    int_id = l1e ;
    // add record
    si.intersect_type = NORMALINT ;
    si.intersect_pos = intersection ;
    si.intersect_id = int_id ;

    global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;
    global_pbetween_ltable.add_record(int_id,l2s,l2e) ;

    return NORMALINT ;    
  }
  
  if(pointOnEdge_bf(l2_start,l1_start,l1_end,l2s,l1s,l1e)) {
    intersection = l2_start ;
    int_id = l2s ;
    // add record
    si.intersect_type = NORMALINT ;
    si.intersect_pos = intersection ;
    si.intersect_id = int_id ;

    global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;
    global_pbetween_ltable.add_record(int_id,l1s,l1e) ;

    return NORMALINT ;    
  }
  
  if(pointOnEdge_bf(l2_end,l1_start,l1_end,l2e,l1s,l1e)) {
    intersection = l2_end ;
    int_id = l2e ;
    // add record
    si.intersect_type = NORMALINT ;
    si.intersect_pos = intersection ;
    si.intersect_id = int_id ;

    global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;
    global_pbetween_ltable.add_record(int_id,l1s,l1e) ;

    return NORMALINT ;    
  }

  // OK by now, we have tested all the improper intersections
  // We need to test is there a proper intersection?
  // We first perform a boolean testing on the possibility
  // of intersection. This is based on the method described
  // in Section 1.5.4. of Joseph O'Rourke's book:
  // Computational Geometry in C (second edition)
  // Here is a quote from the book: 
  // if two segments ab and cd intersect in
  // their interiors, then c and d are split by the line
  // containing ab: c is to one side and d to the other.
  // And likewise, a and b are split by the line containing
  // cd. Both conditions must be true.
  // Here we use the exact predicate to do testing
  if(!(xor_bool(left_exact(l1_start,l1_end,l2_start),
                left_exact(l1_start,l1_end,l2_end)) &&
       xor_bool(left_exact(l2_start,l2_end,l1_start),
                left_exact(l2_start,l2_end,l1_end)))) {
    // no intersection
    si.intersect_type = NOINT ;
    global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;

    return NOINT ;
  }
  // NOTE: if the precondition of this function is met,
  // i.e., at least three points of (l1s,l1e,l2s,l2e)
  // are on the boundary edge. Then a proper intersection
  // is impossible unless the grid is extremely simple
  // and small. In that a segment could "cross" the grid,
  // i.e., from one boundary edge to another boundary edge.
  
  // then we compute the intersection
  // and this is done in INexact floating-point computation
  double denominator =
    (l2_end.y - l2_start.y) * (l1_end.x - l1_start.x) -
    (l2_end.x - l2_start.x) * (l1_end.y - l1_start.y) ;

  double numerator1 =
    (l2_end.x - l2_start.x) * (l1_start.y - l2_start.y) -
    (l2_end.y - l2_start.y) * (l1_start.x - l2_start.x) ;
  /*
  double numerator2 =
    (l1_end.x - l1_start.x) * (l1_start.y - l2_start.y) -
    (l1_end.y - l1_start.y) * (l1_start.x - l2_start.x) ;
  */
  double u1 = numerator1 / denominator ;
  //double u2 = numerator2 / denominator ;

  intersection = l1_start + u1 * (l1_end - l1_start) ;
  // Now, we have a true valid NEW intersection point
  
  // Since this is a new node, we need to create a
  // new index for this node. And we need to increment
  // the mesh_new_node_index too.
  int_id = mesh_new_node_index++ ;
  
  // add record
  si.intersect_type = NORMALINT ;
  si.intersect_pos = intersection ;
  si.intersect_id = int_id ;
  
  global_ssint_ltable(l1s,l1e,l2s,l2e) = si ;
  //global_ssint_ltable(l2s,l2e,l1s,l1e) = si ;
  // Since this is a proper intersection, the intersection
  // point p lies on both segments. Therefore, we also
  // need to put the record into the global_pbetween_ltable
  global_pbetween_ltable.add_record(int_id,l1s,l1e) ;
  global_pbetween_ltable.add_record(int_id,l2s,l2e) ;
  
  return NORMALINT ;
}
// a special version of poly_and_poly. This one is used
// for polygons that lie on the mesh boundary outer edges
bool poly_and_poly_bf(const vector<vec2d>& P,
                      const vector<vec2d>& Q,
                      const vector<Entity>& Pid,
                      const vector<Entity>& Qid,
                      vector<vec2d>& Ins,
                      vector<Entity>& Ins_index) {
#ifdef DEBUG_POLY_AND_POLY
  ////////////////////////
  cerr << "This is boundary version of poly_and_poly" << endl ;
  ////////////////////////
#endif
  typedef vector<vec2d>::size_type size ;
  size P_size = P.size() ;
  size Q_size = Q.size() ;
  size P_size2 = P_size*2 ;
  size Q_size2 = Q_size*2 ;
  size bug_P = 0, bug_Q = 0 ;
  size P_adv = 0, Q_adv = 0 ;
  inside_type inflag = UNKNOWN ;
  // record the first intersection point
  vec2d first_int(0.0,0.0) ;
  Entity first_int_id = -1 ;
  Entity last_output_id = -1 ;
  bool firstIntersect = true ;
  // flag indicates whether the intersection in previous iteration
  // is equal to the first intersection
  bool is_ins_last_e_first = false ;
  // data used to control output of collinear points
  bool continue_collinear = false ;
  Entity collinear_point_id = -1 ;
  vec2d collinear_point_pos(0.0, 0.0) ;

  do {
    const vec2d& P_cur = P[bug_P] ;
    const vec2d& Q_cur = Q[bug_Q] ;
    size bug_P_pre = (bug_P - 1 + P_size) % P_size ;
    size bug_Q_pre = (bug_Q - 1 + Q_size) % Q_size ;
    const vec2d& P_pre = P[bug_P_pre] ;
    const vec2d& Q_pre = Q[bug_Q_pre] ;

    Entity P_cur_id = Pid[bug_P] ;
    Entity P_pre_id = Pid[bug_P_pre] ;
    Entity Q_cur_id = Qid[bug_Q] ;
    Entity Q_pre_id = Qid[bug_Q_pre] ;

    // we will need to determine which version of
    // seg_seg_int should we use? We check the number of
    // points in (P_cur_id,P_pre_id,Q_cur_id,Q_pre_id)
    // that are on the boundary edge. If if is >= 3,
    // then we should call the special version seg_seg_int_bf
    vec2d seg_int ;
    Entity seg_int_id ;
    seg_int_type seg_int_result ;
    int bn_num = 0 ;
    if(masterbc_edge_nodes.inSet(P_cur_id))
      ++bn_num ;
    if(masterbc_edge_nodes.inSet(P_pre_id))
      ++bn_num ;
    if(slavebc_edge_nodes.inSet(Q_cur_id))
      ++bn_num ;
    if(slavebc_edge_nodes.inSet(Q_pre_id))
      ++bn_num ;

    if(bn_num >= 3)
      seg_int_result =
        seg_seg_int_bf(P_cur,P_pre,Q_cur,Q_pre,
                       P_cur_id,P_pre_id,Q_cur_id,Q_pre_id,
                       seg_int,seg_int_id) ;
    else
      seg_int_result =
        seg_seg_int(P_cur,P_pre,Q_cur,Q_pre,
                    P_cur_id,P_pre_id,Q_cur_id,Q_pre_id,
                    seg_int,seg_int_id) ;

    int P_in_Q_plane = in_half_plane(Q_cur,Q_pre,P_cur,
                                     Q_cur_id,Q_pre_id,P_cur_id) ;
    int Q_in_P_plane = in_half_plane(P_cur,P_pre,Q_cur,
                                     P_cur_id,P_pre_id,Q_cur_id) ;

#ifdef DEBUG_POLY_AND_POLY
    ////////////////////////
    cerr.precision(32) ;
    cerr << "P_cur: " << P_cur << " P_pre: " << P_pre << endl ;
    cerr << "Q_cur: " << Q_cur << " Q_pre: " << Q_pre << endl ;
    cerr << "P_cur_id: " << P_cur_id << " P_pre_id: " << P_pre_id << endl ;
    cerr << "Q_cur_id: " << Q_cur_id << " Q_pre_id: " << Q_pre_id << endl ;
    cerr << "P_in_Q_plane: " << P_in_Q_plane << endl ;
    cerr << "Q_in_P_plane: " << Q_in_P_plane << endl ;

    if(seg_int_result == NOINT)
      cerr << "No intersection!" <<endl ;
    else if(seg_int_result == NORMALINT) {
      cerr << "Normal intersection!" << endl ;
      cerr << "intersection point: " << seg_int << endl ;
      cerr << "intersection id: " << seg_int_id << endl ;
      cerr << "last_output_id: " << last_output_id << endl ;
      cerr << "first_int_id: " << first_int_id << endl ;
      cerr << "is_ins_last_e_first: " << is_ins_last_e_first << endl ;
    }
    else if(seg_int_result == COLLINEAR)
      cerr << "Collinear!" << endl ;
    else
      cerr << "Unknown result!" << endl ;
    ///////////////////////
#endif
    if(seg_int_result == NORMALINT) {
      // check whether we could quit
      if(!is_ins_last_e_first && !firstIntersect) {
        if(seg_int_id == first_int_id)
          break ;
      }
      // if first intersection
      if( (inflag == UNKNOWN) && firstIntersect) {
        P_adv = Q_adv = 0 ;
        firstIntersect = false ;
        first_int = seg_int ;
        first_int_id = seg_int_id ;
      }
      if(seg_int_id == first_int_id)
        is_ins_last_e_first = true ;
      else
        is_ins_last_e_first = false ;

      if(seg_int_id != last_output_id) {
        Ins.push_back(seg_int) ;
        Ins_index.push_back(seg_int_id) ;
        last_output_id = seg_int_id ;
      }

      if(P_in_Q_plane > 0)
        inflag = PIN ;
      else if(Q_in_P_plane > 0)
        inflag = QIN ;
    } else {
      is_ins_last_e_first = false ;
    }
    // advance P or Q
    vec2d P_edge = P_cur - P_pre ;
    vec2d Q_edge = Q_cur - Q_pre ;
    int cross_sign = cross2d_exact(Q_cur,Q_pre,P_cur,P_pre) ;
#ifdef DEBUG_POLY_AND_POLY
    ///////////////////////////
    if(inflag == UNKNOWN)
      cerr << "inflag == UNKNOWN" << endl ;
    else if(inflag == PIN)
      cerr << "inflag == PIN" << endl ;
    else if(inflag == QIN)
      cerr << "inflag == QIN" << endl ;
    else
      cerr << "UNKNOWN inflag" << endl ;

    cerr << "cross_sign: " << cross_sign << endl ;
    cerr << "dot(P_edge,Q_edge): " << dot(P_edge,Q_edge) << endl ;
    ///////////////////////////
#endif
    if( (seg_int_result==COLLINEAR) && (dot(P_edge,Q_edge) < 0)) {
      // two edges are colinear and in oppisite position
      // this means the two polygons only share an edge
      // we therefore could exit the loop now.
      // NOTE: we could ouput the shared edges, but since
      // we do not regard them as valid intersection here,
      // we just return
      Ins.clear() ;
      Ins_index.clear() ;
      return false ;
    } else if( (seg_int_result==PARALLEL) &&
               (P_in_Q_plane < 0) && (Q_in_P_plane < 0)) {
      // two edges are parallel and separated, the two
      // polygons are therefore disjoint
      Ins.clear() ;
      Ins_index.clear() ;
      return false ;
    } else if(seg_int_result==COLLINEAR) {
      // special case A&B collinear and have the same direction
      advance_collinear_bf(P_cur, Q_cur, P_cur_id, Q_cur_id, P_edge,
                           continue_collinear, collinear_point_pos,
                           collinear_point_id,
                           firstIntersect, first_int_id,
                           last_output_id, bug_P, bug_Q,
                           P_adv, Q_adv, P_size, Q_size, Ins, Ins_index) ;
    } else if(cross_sign >= 0) {
      continue_collinear = false ;
      if(P_in_Q_plane > 0)
        advance_Q(inflag,Q_size,Q_cur,Q_cur_id,first_int_id,
                  last_output_id,bug_Q,Q_adv,Ins,Ins_index) ;
      else
        advance_P(inflag,P_size,P_cur,P_cur_id,first_int_id,
                  last_output_id,bug_P,P_adv,Ins,Ins_index) ;
    } else {
      continue_collinear = false ;
      if(Q_in_P_plane > 0)
        advance_P(inflag,P_size,P_cur,P_cur_id,first_int_id,
                  last_output_id,bug_P,P_adv,Ins,Ins_index) ;
      else
        advance_Q(inflag,Q_size,Q_cur,Q_cur_id,first_int_id,
                  last_output_id,bug_Q,Q_adv,Ins,Ins_index) ;
    }
  } while( ((P_adv < P_size) || (Q_adv < Q_size)) &&
           (P_adv < P_size2) && (Q_adv < Q_size2)) ;

#ifdef DEBUG_POLY_AND_POLY
  /////////////////
  if(Ins.size() != Ins_index.size()) {
    cerr << "ERROR: Ins and Ins_index have different size!!!" << endl ;
  }
  cerr << "intersection: " << endl ;
  for(vector<vec2d>::const_iterator vi=Ins.begin();
      vi!=Ins.end();++vi)
    cerr << *vi << endl ;
  cerr << "intersection index: " << endl ;
  for(vector<Entity>::const_iterator vi=Ins_index.begin();
      vi!=Ins_index.end();++vi)
    cerr << *vi << endl ;
  cerr << endl ;
  /////////////////
#endif
  if(!Ins.empty()) {
    if(!collinear(Ins,Ins_index)) {
      return true ;
    }
    // we do not regard an edge or a vertex
    // as a valid intersection
    else {
      Ins.clear() ;
      Ins_index.clear() ;
#ifdef DEBUG_POLY_AND_POLY
      /////////////////
      cerr << "No intersection!" << endl ;
      cerr << endl ;
      /////////////////
#endif
      return false ;
    }
  }

  if(poly_vertex_in(P,Q)) {
    // P in Q completely
    Ins = P ;
    Ins_index = Pid ;
#ifdef DEBUG_POLY_AND_POLY
    /////////////////
    cerr << "Intersection is P!" << endl ;
    cerr << endl ;
    /////////////////
#endif
    return true ;
  } else if(poly_vertex_in(Q,P)) {
    // Q in P completely
    Ins = Q ;
    Ins_index = Qid ;
#ifdef DEBUG_POLY_AND_POLY
    /////////////////
    cerr << "Intersection is Q!" << endl ;
    cerr << endl ;
    /////////////////
#endif
    return true ;
  } else {
    // no intersection
    Ins.clear() ;
    Ins_index.clear() ;
#ifdef DEBUG_POLY_AND_POLY
    /////////////////
    cerr << "No intersection!" << endl ;
    cerr << endl ;
    /////////////////
#endif
    return false ;
  }
}
/************* end of special polygon intersection *************/

// check counter clockwise order for convex polygon
bool check_ccw(const vector<vec2d>& P) {
  typedef vector<vec2d>::size_type size ;
  size p_size = P.size() ;
  size i ;
  for(i=0;i<p_size;++i) {
    int x = (i-1+p_size)%p_size ;
    int y = i ;
    int z = (i+1)%p_size ;
    const vec2d& p_pre = P[x] ;
    const vec2d& p_cur = P[y] ;
    const vec2d& p_nxt = P[z] ;

    double a[2] = {p_pre.x, p_pre.y} ;
    double b[2] = {p_cur.x, p_cur.y} ;
    double c[2] = {p_nxt.x, p_nxt.y} ;
    double o = orient2d(a,b,c) ;
    if(o < 0.0)
      return false ;
  }
  return true ;
}
// the approximate version of check_ccw
// created for boundary faces
// the difference is that this one performs
// approximate collinear testing
bool check_ccw_approx(const vector<vec2d>& P) {
  typedef vector<vec2d>::size_type size ;
  size p_size = P.size() ;
  size i ;
  for(i=0;i<p_size;++i) {
    int x = (i-1+p_size)%p_size ;
    int y = i ;
    int z = (i+1)%p_size ;
    const vec2d& p_pre = P[x] ;
    const vec2d& p_cur = P[y] ;
    const vec2d& p_nxt = P[z] ;
    // first we perform approximate collinear testing
    if(collinear3_approx(p_pre,p_cur,p_nxt))
      continue ;
    // otherwise, do the normal testing
    double a[2] = {p_pre.x, p_pre.y} ;
    double b[2] = {p_cur.x, p_cur.y} ;
    double c[2] = {p_nxt.x, p_nxt.y} ;
    double o = orient2d(a,b,c) ;
    if(o < 0.0)
      return false ;
  }
  return true ;
}
// check for convexity of a polygon
bool check_convex(const vector<vec2d>& P) {
  typedef vector<vec2d>::size_type size ;
  size p_size = P.size() ;
  size i ;
  int flag = 0 ;
  for(i=0;i<p_size;++i) {
    int x = (i-1+p_size)%p_size ;
    int y = i ;
    int z = (i+1)%p_size ;
    const vec2d& p_pre = P[x] ;
    const vec2d& p_cur = P[y] ;
    const vec2d& p_nxt = P[z] ;
    
    double a[2] = {p_pre.x, p_pre.y} ;
    double b[2] = {p_cur.x, p_cur.y} ;
    double c[2] = {p_nxt.x, p_nxt.y} ;

    double o = orient2d(a,b,c) ;
    if(o > 0.0) {
      if(flag == -1)
        return false ;
      flag = 1 ;
    } else if(o < 0.0) {
      if(flag == 1)
        return false ;
      flag = -1 ;
    }
  }
  return true ;
}
// the approximate version of check_convex
// created for boundary faces
// the difference is that this one performs
// approximate collinear testing
bool check_convex_approx(const vector<vec2d>& P) {
  typedef vector<vec2d>::size_type size ;
  size p_size = P.size() ;
  size i ;
  int flag = 0 ;
  for(i=0;i<p_size;++i) {
    int x = (i-1+p_size)%p_size ;
    int y = i ;
    int z = (i+1)%p_size ;
    const vec2d& p_pre = P[x] ;
    const vec2d& p_cur = P[y] ;
    const vec2d& p_nxt = P[z] ;
    // first we perform approximate collinear testing
    if(collinear3_approx(p_pre,p_cur,p_nxt))
      continue ;
    // otherwise, do the normal testing
    double a[2] = {p_pre.x, p_pre.y} ;
    double b[2] = {p_cur.x, p_cur.y} ;
    double c[2] = {p_nxt.x, p_nxt.y} ;

    double o = orient2d(a,b,c) ;
    if(o > 0.0) {
      if(flag == -1)
        return false ;
      flag = 1 ;
    } else if(o < 0.0) {
      if(flag == 1)
        return false ;
      flag = -1 ;
    }
  }
  return true ;
}
// wrapper of poly_and_poly function
// this is the actual function that is
// called in the face splitting
bool poly_and_poly(const vector<Entity>& P_nodes,
                   const store<vec2d>& P_pos,
                   const vector<Entity>& Q_nodes,
                   const store<vec2d>& Q_pos,
                   vector<vec2d>& Ins,
                   vector<Entity>& Ins_index) {
  // check check validity of polygon
  if( (P_nodes.size() < 3) || (Q_nodes.size() < 3)) {
    cerr << "ERROR: polygon in poly_and_poly ill defined!" << endl ;
    Loci::Abort() ;
  }
  // NOTE: we will need to check which version of
  // poly_and_poly we will need to call.
  // We will need to call poly_and_poly_bf if two
  // polygons are both and boundary edge polygon.
  // A polygon is considered a boundary edge polygon
  // if any of its nodes is on the boundary edge
  bool P_bf = false ;
  bool Q_bf = false ;
  vector<vec2d> p_poly, q_poly ;
  vector<Entity> pid, qid ;
  for(vector<Entity>::const_iterator vi=P_nodes.begin();
      vi!=P_nodes.end();++vi) {
    if(masterbc_edge_nodes.inSet(*vi))
      P_bf = true ;
    p_poly.push_back(P_pos[*vi]) ;
  }
  for(vector<Entity>::const_iterator vi=Q_nodes.begin();
      vi!=Q_nodes.end();++vi) {
    if(slavebc_edge_nodes.inSet(*vi))
      Q_bf = true ;
    q_poly.push_back(Q_pos[*vi]) ;
  }
  // Here we first perform a quick bounding box
  // test, i.e., if the bounding box of P and Q
  // does not intersect, then P and Q are surely
  // not intersected
  rectangle p_box, q_box ;
  nodes_bbox(p_poly,p_box) ;
  nodes_bbox(q_poly,q_box) ;
  if(!p_box.rectangle_int(q_box))
    return false ;

  bool p_convex, q_convex ;
  if(P_bf)
    p_convex = check_convex_approx(p_poly) ;
  else
    p_convex = check_convex(p_poly) ;

  if(Q_bf)
    q_convex = check_convex_approx(q_poly) ;
  else
    q_convex = check_convex(q_poly) ;
  
#ifdef POLY_AND_POLY_CONVEXITY_CHECK
  if(!p_convex) {
    cerr << endl ;
    cerr.precision(16) ;
    cerr << "non-convex polygon : " << endl ;
    for(size_t i=0;i!=P_nodes.size();++i) {
      Entity n = P_nodes[i] ;
      cerr << n <<  ": ("
           << masterbc_edge_nodes.inSet(n)
           << ") = " << p_poly[i] << endl ;
    }
  }
  if(!q_convex) {
    cerr << endl ;
    cerr.precision(16) ;
    cerr << "non-convex polygon : " << endl ;
    for(size_t i=0;i!=Q_nodes.size();++i) {
      Entity n = Q_nodes[i] ;
      cerr << n <<  ": ("
           << masterbc_edge_nodes.inSet(n)
           << ") = " << q_poly[i] << endl ;
    }
  }
#endif
  
  // check if convex
  if(!p_convex || !q_convex) {
    cerr << "ERROR: polygon(s) are not convex in poly_and_poly!" << endl ;
    Loci::Abort() ;
  }

  bool p_ccw, q_ccw ;
  if(P_bf)
    p_ccw = check_ccw_approx(p_poly) ;
  else
    p_ccw = check_ccw(p_poly) ;

  if(Q_bf)
    q_ccw = check_ccw_approx(q_poly) ;
  else
    q_ccw = check_ccw(q_poly) ;
  
  pid = P_nodes ; qid = Q_nodes ;
  // check for ccw
  if(!p_ccw) {
    //cout << "Warning: polygon not defined counter clockwise, reversed!"
    //   << endl ;
    reverse(p_poly.begin(),p_poly.end()) ;
    reverse(pid.begin(),pid.end()) ;
  }
  if(!q_ccw) {
    //cout << "Warning: polygon not defined counter clockwise, reversed!"
    //   << endl ;
    reverse(q_poly.begin(),q_poly.end()) ;
    reverse(qid.begin(),qid.end()) ;
  }
  // now ready for intersection
  if(P_bf && Q_bf)
    return poly_and_poly_bf(p_poly,q_poly,pid,qid,Ins,Ins_index) ;
  else
    return poly_and_poly(p_poly,q_poly,pid,qid,Ins,Ins_index) ;
}
/**
 ************************************************************
 *                   FaceSet QuadTree codes                 *
 ************************************************************
 */
class QuadTree {
  bool is_empty ; // whether the tree is empty
  class qtnode {
  public:
    int depth ;
    int max_depth ;
    int face_set_size ;
    // indicate whether this node contains data
    // or just contains pointers to four childrens
    bool data_node ;
    rectangle bbox ; // bounding box for this node
    entitySet faces ; // face list if "data_node" is true
    // four of its children
    qtnode* top_left ;
    qtnode* top_right ;
    qtnode* bot_left ;
    qtnode* bot_right ;
    // bounding box for four children
    rectangle tlb, trb, blb, brb ;
    
    qtnode() ;
    qtnode(int md, int fs, const rectangle& b) ;
    // set the bounding box of its four children
    void set_children_bbox() {
      double half_x = (bbox.l.x + bbox.r.x) / 2 ;
      double half_y = (bbox.l.y + bbox.r.y) / 2 ;
      tlb = rectangle(bbox.l, vec2d(half_x,half_y)) ;
      trb = rectangle(vec2d(half_x,bbox.l.y),
                      vec2d(bbox.r.x,half_y)) ;
      blb = rectangle(vec2d(bbox.l.x,half_y),
                      vec2d(half_x,bbox.r.y)) ;
      brb = rectangle(tlb.r, bbox.r) ;
    }
    void erase_node() ;    
    void copy_node(const qtnode& n) ;
    // the copy constructor, assignment operator, and destructor
    qtnode(const qtnode& n) {
      copy_node(n) ;
    }
    qtnode& operator=(const qtnode& qtn) {
      if(&qtn != this) {
        erase_node() ;
        copy_node(qtn) ;
      }
      return *this ;
    }
    ~qtnode() {
      erase_node() ;
    }
    // function that builds the node
    void build_node(const entitySet& allfaces,
                    const multiMap& face2node,
                    const store<vec2d>& pos,
                    int current_depth) ;
    // locate a given point in the tree, return
    // the potential faces that contain the point.
    // Normally only one face will contain the point,
    // but if the point lies on a vertex or an edge
    // then multiple faces could all contain the point
    entitySet locate_point(const vec2d& p) const ;
    // return the maximum depth of the tree
    int tree_depth() const ;
    // return the maximum number of face set size of any subtree
    int max_face_size() const ;
    // return the number of nonempty subtrees
    int leafs_num() const ;
    // return all the data contained in this node
    entitySet leafs() const ;
    // return the sum of size in all leaves
    int all_leafs_size() const ;
    // locate the given rectangle in the tree. return
    // all the faces that contained in the subspaces that
    // touched by the rectangle
    entitySet locate_rectangle(const rectangle& r) const ;
  } ;
  
  qtnode* root ;
  
public:
  //  QuadTree() {is_empty = true ; root = new qtnode ;}
  QuadTree(int md, int fs, const rectangle& b)
  {is_empty = true ; root = new qtnode(md,fs,b) ;}
  // copy constructor, assignment operator, and the destructor
  QuadTree(const QuadTree& t) {
    is_empty = t.is_empty ;
    root->erase_node() ;
    root->copy_node(*(t.root)) ;
  }
  QuadTree& operator=(const QuadTree& t) {
    if(&t != this) {
      is_empty = t.is_empty ;
      root->erase_node() ;
      root->copy_node(*(t.root)) ;
    }
    return *this ;
  }
  ~QuadTree()
  {root->erase_node() ; delete root ; is_empty=true ;}
  // public utility methods
  bool empty() const {return is_empty ;}
  void build_tree(const entitySet& allfaces,
                  const multiMap& face2node,
                  const store<vec2d>& pos) {
    root->build_node(allfaces,face2node,pos,0) ;
    is_empty=false ;
  }
  entitySet locate_point(const vec2d& p) const
  {return root->locate_point(p) ;}
  int depth() const
  {return root->tree_depth() ;}
  int max_face_size() const
  {return root->max_face_size() ;}
  int leafs_num() const
  {return root->leafs_num() ;}
  entitySet locate_rectangle(const rectangle& r) const {
    // check validity of the rectangle
    r.check_valid() ;
    return root->locate_rectangle(r) ;
  }
  // locating a given polygon
  entitySet locate_polygon(const vector<int>& nodes,
                           const store<vec2d>& pos) const ;
  // return all the data in the tree
  entitySet all_leafs() const
  {return root->leafs() ;}
  // return the sum of all leaf size
  int all_leafs_size() const
  {return root->all_leafs_size() ;}
} ;

/////////////////// QuadTree class implementation ////////////////////
QuadTree::qtnode::qtnode() {
  depth = 0 ;
  max_depth = 20 ;
  face_set_size = 10 ;
  data_node = false ;
  faces = EMPTY ;
  top_left = 0 ;
  top_right = 0 ;
  bot_left = 0 ;
  bot_right = 0 ;
}
QuadTree::qtnode::qtnode(int md, int fs, const rectangle& b) {
  depth = 0 ;
  max_depth = md ;
  face_set_size = fs ;
  data_node = false ;
  faces = EMPTY ;
  top_left = top_right = bot_left = bot_right = 0 ;
  bbox = b ;
  set_children_bbox() ;
}
void QuadTree::qtnode::erase_node() {
  if(top_left != 0) {
    if(top_left->data_node) {
      delete top_left ;
      top_left = 0 ;
    }
    else
      top_left->erase_node() ;
  }
  
  if(top_right != 0) {
    if(top_right->data_node) {
      delete top_right ;
      top_right = 0 ;
    }
    else
      top_right->erase_node() ;
  }

  if(bot_left != 0) {
    if(bot_left->data_node) {
      delete bot_left ;
      bot_left = 0 ;
    }
    else
      bot_left->erase_node() ;
  }

  if(bot_right != 0) {
    if(bot_right->data_node) {
      delete bot_right ;
      bot_right = 0 ;
    }
    else
      bot_right->erase_node() ;
  }
}
void QuadTree::qtnode::copy_node(const qtnode& n) {
  erase_node() ;
  
  depth = n.depth ;
  max_depth = n.max_depth ;
  face_set_size = n.face_set_size ;
  data_node = n.data_node ;
  bbox = n.bbox ;
  faces = n.faces ;
  tlb = n.tlb ;
  trb = n.trb ;
  blb = n.blb ;
  brb = n.brb ;
  
  if(!data_node) {
    top_left = new qtnode ;
    top_right = new qtnode ;
    bot_left = new qtnode ;
    bot_right = new qtnode ;
    top_left->copy_node(*(n.top_left)) ;
    top_right->copy_node(*(n.top_right)) ;
    bot_left->copy_node(*(n.bot_left)) ;
    bot_right->copy_node(*(n.bot_right)) ;
  }else {
    top_left = 0;
    top_right = 0 ;
    bot_left = 0 ;
    bot_right = 0 ;
  }
}

void QuadTree::qtnode::build_node(const entitySet& allfaces,
                                  const multiMap& face2node,
                                  const store<vec2d>& pos,
                                  int current_depth) {
  int allfaces_size = allfaces.size() ;
  if(allfaces_size == 0) {
    cerr << "ERROR: empty faces in build_node! This should not occur!"
         << endl ;
    Loci::Abort() ;
  }
  depth = current_depth ;
  data_node = true ;
  if(depth >= max_depth) {
    faces = allfaces ;
    return ;
  }
  if(allfaces_size <= face_set_size) {
    faces = allfaces ;    
    return ;
  }

  faces = EMPTY ;
  data_node = false ;
  entitySet tl,tr,bl,br ;
  for(entitySet::const_iterator ei=allfaces.begin();
      ei!=allfaces.end();++ei) {
    vector<int> nodes = get_face_nodes(face2node,*ei) ;
    // test where this face falls into
    if(tlb.is_convex_polygon_int(nodes,pos))
      tl+=*ei ;
    if(trb.is_convex_polygon_int(nodes,pos))
      tr+=*ei ;
    if(blb.is_convex_polygon_int(nodes,pos))
      bl+=*ei ;
    if(brb.is_convex_polygon_int(nodes,pos))
      br+=*ei ;
  }
  int tl_size = tl.size() ;
  int tr_size = tr.size() ;
  int bl_size = bl.size() ;
  int br_size = br.size() ;
  // build subnodes
  if(tl_size > 0) {
    top_left = new qtnode(max_depth,face_set_size,tlb) ;
    top_left->build_node(tl,face2node,pos,depth+1) ;
  } else
    top_left = 0 ;
  
  if(tr_size > 0) {
    top_right = new qtnode(max_depth,face_set_size,trb) ;
    top_right->build_node(tr,face2node,pos,depth+1) ;
  } else
    top_right = 0 ;
  
  if(bl_size > 0) {
    bot_left = new qtnode(max_depth,face_set_size,blb) ;
    bot_left->build_node(bl,face2node,pos,depth+1) ;
  } else
    bot_left = 0 ;
  
  if(br_size > 0) {
    bot_right = new qtnode(max_depth,face_set_size,brb) ;
    bot_right->build_node(br,face2node,pos,depth+1) ;
  } else
    bot_right = 0 ;
}

entitySet QuadTree::qtnode::locate_point(const vec2d& p) const {
  entitySet ret_faces ;
  if(bbox.is_point_in(p)) {
    if(data_node)
      ret_faces += faces ;
    else {
      // recursively testing subnodes
      if(top_left != 0)
        if(tlb.is_point_in(p))
          ret_faces += top_left->locate_point(p) ;

      if(top_right != 0)
        if(trb.is_point_in(p))
          ret_faces += top_right->locate_point(p) ;

      if(bot_left != 0)
        if(blb.is_point_in(p))
          ret_faces += bot_left->locate_point(p) ;

      if(bot_right != 0)
        if(brb.is_point_in(p))
          ret_faces += bot_right->locate_point(p) ;
    }
  }

  return ret_faces ;
}
int QuadTree::qtnode::tree_depth() const {
  if(data_node)
    return depth ;
  else {
    int tld=0, trd=0, bld=0, brd=0 ;
    if(top_left != 0)
      tld = top_left->tree_depth() ;
    if(top_right != 0)
      trd = top_right->tree_depth() ;
    if(bot_left != 0)
      bld = bot_left->tree_depth() ;
    if(bot_right != 0)
      brd = bot_right->tree_depth() ;
    
    return max(max(tld,trd),max(bld,brd)) ;
  }
}
int QuadTree::qtnode::max_face_size() const {
  if(data_node)
    return faces.size() ;
  else {
    int tlf=0, trf=0, blf=0, brf=0 ;
    if(top_left != 0)
      tlf = top_left->max_face_size() ;
    if(top_right != 0)
      trf = top_right->max_face_size() ;
    if(bot_left != 0)
      blf = bot_left->max_face_size() ;
    if(bot_right != 0)
      brf = bot_right->max_face_size() ;
    
    return max(max(tlf,trf),max(blf,brf)) ;
  }
}
int QuadTree::qtnode::leafs_num() const {
  if(data_node)
    return 1 ;
  else {
    int tle=0, tre=0, ble=0, bre=0 ;
    if(top_left != 0)
      tle = top_left->leafs_num() ;
    if(top_right != 0)
      tre = top_right->leafs_num() ;
    if(bot_left != 0)
      ble = bot_left->leafs_num() ;
    if(bot_right != 0)
      bre = bot_right->leafs_num() ;
    
    return tle+tre+ble+bre ;
  }
}
entitySet QuadTree::qtnode::leafs() const {
  entitySet data ;
  if(data_node)
    data += faces ;
  else {
    // recursively collect subnodes data
    if(top_left != 0)
      data += top_left->leafs() ;
    
    if(top_right != 0)
      data += top_right->leafs() ;
    
    if(bot_left != 0)
      data += bot_left->leafs() ;
    
    if(bot_right != 0)
      data += bot_right->leafs() ;
  }
  return data ;
}
int QuadTree::qtnode::all_leafs_size() const {
  int s = 0 ;
  if(data_node)
    return faces.size() ;
  else {
    if(top_left != 0)
      s += top_left->all_leafs_size() ;
    if(top_right != 0)
      s += top_right->all_leafs_size() ;
    if(bot_left != 0)
      s += bot_left->all_leafs_size() ;
    if(bot_right != 0)
      s += bot_right->all_leafs_size() ;
  }
  return s ;
}
entitySet QuadTree::qtnode::
locate_rectangle(const rectangle& r) const {
  entitySet data ;
  if(data_node)
    data += faces ;
  else {
    // splitting the passed in rectangle
    rectangle rsub ;
    if(tlb.rectangle_int(r,rsub))
      if(top_left != 0)
        data += top_left->locate_rectangle(rsub) ;

    if(trb.rectangle_int(r,rsub))
      if(top_right != 0)
        data += top_right->locate_rectangle(rsub) ;

    if(blb.rectangle_int(r,rsub))
      if(bot_left != 0)
        data += bot_left->locate_rectangle(rsub) ;

    if(brb.rectangle_int(r,rsub))
      if(bot_right != 0)
        data += bot_right->locate_rectangle(rsub) ;
  }
  return data ;
}
// we first compute a bounding box for the given polygon
// then call locate_rectangle to locate the bounding box
entitySet QuadTree::locate_polygon(const vector<int>& nodes,
                                   const store<vec2d>& pos) const {
  if(nodes.empty()) {
    cerr << "Error: empty polygon in locate_polygon!" << endl ;
    Loci::Abort() ;
  }
//   // first we construct a store of position
//   vector<vec2d> nodes_pos ;
//   for(vector<int>::const_iterator ni=nodes.begin();
//       ni!=nodes.end();++ni)
//     nodes_pos.push_back(pos[*ni]) ;
//   rectangle b ;
//   nodes_bbox(nodes_pos,b) ;
  rectangle b ;
  nodes_bbox(nodes,pos,b) ;
  return locate_rectangle(b) ;
}
///////////end of FaceSet QuadTree implementation////////////

/**
 ************************************************************
 *                  PointSet QuadTree codes                 *
 ************************************************************
 */
class PointSet_QuadTree {
  bool is_empty ; // whether the tree is empty
  class qtnode {
  public:
    int depth ;
    // indicate whether this node contains data
    // or just contains pointers to four childrens
    bool data_node ;
    rectangle bbox ; // bounding box for this node
    vec2d point ; // the data (point position)
    entitySet id ; // the point's id, it could be a set that
                   // actually indicates there are multiple
                   // points with the same position
    // four of its children
    qtnode* top_left ;
    qtnode* top_right ;
    qtnode* bot_left ;
    qtnode* bot_right ;
    // bounding box for four children
    rectangle tlb, trb, blb, brb ;
    
    qtnode() ;
    qtnode(const rectangle& b) ;
    // set the bounding box of its four children
    void set_children_bbox() {
      double half_x = (bbox.l.x + bbox.r.x) / 2 ;
      double half_y = (bbox.l.y + bbox.r.y) / 2 ;
      tlb = rectangle(bbox.l, vec2d(half_x,half_y)) ;
      trb = rectangle(vec2d(half_x,bbox.l.y),
                      vec2d(bbox.r.x,half_y)) ;
      blb = rectangle(vec2d(bbox.l.x,half_y),
                      vec2d(half_x,bbox.r.y)) ;
      brb = rectangle(tlb.r, bbox.r) ;
    }
    void erase_node() ;    
    void copy_node(const qtnode& n) ;
    // the copy constructor, assignment operator, and destructor
    qtnode(const qtnode& n) {
      copy_node(n) ;
    }
    qtnode& operator=(const qtnode& qtn) {
      if(&qtn != this) {
        erase_node() ;
        copy_node(qtn) ;
      }
      return *this ;
    }
    ~qtnode() {
      erase_node() ;
    }
    // function that builds the node
    template<class STOREOFVEC2D>
    void build_node(const entitySet& allpoints,
                    const STOREOFVEC2D& pos,
                    int current_depth) ;
    // return the maximum depth of the tree
    int tree_depth() const ;
    // return the number of nonempty subtrees
    int leafs_num() const ;
    // return all the data contained within this subspace
    entitySet leafs() const ;
    // locate the given rectangle in the tree. return
    // all the points that contained in the subspaces that
    // touched by the rectangle
    entitySet locate_rectangle(const rectangle& r) const ;
    // locate the given circle (disk) in the tree. return
    // all the points taht contained in the subspaces that
    // touched by the disk. This effectively returns all
    // candidate points that are within the range "radius"
    // of the circle's "center"
    entitySet locate_circle(const vec2d& center, double radius) const ;
    // insert a point into the tree
    void insert_point(const vec2d& p, Entity identity) ;
  } ;
  
  qtnode* root ;
  
public:
  PointSet_QuadTree(const rectangle& b)
  {is_empty = true ; root = new qtnode(b) ;}
  // copy constructor, assignment operator, and the destructor
  PointSet_QuadTree(const PointSet_QuadTree& t) {
    is_empty = t.is_empty ;
    root->erase_node() ;
    root->copy_node(*(t.root)) ;
  }
  PointSet_QuadTree& operator=(const PointSet_QuadTree& t) {
    if(&t != this) {
      is_empty = t.is_empty ;
      root->erase_node() ;
      root->copy_node(*(t.root)) ;
    }
    return *this ;
  }
  ~PointSet_QuadTree()
  {root->erase_node() ; delete root ; is_empty=true ;}
  // public utility methods
  bool empty() const {return is_empty ;}
  void build_tree(const entitySet& allpoints,
                  const store<vec2d>& pos) {
    root->build_node(allpoints,pos,0) ;
    is_empty=false ;
  }
  int depth() const
  {return root->tree_depth() ;}
  int leafs_num() const
  {return root->leafs_num() ;}
  entitySet locate_rectangle(const rectangle& r) const {
    // check validity of the rectangle
    r.check_valid() ;
    return root->locate_rectangle(r) ;
  }
  entitySet locate_circle(const vec2d& center, double radius) const
  {return root->locate_circle(center,radius) ;}
  // locating a given polygon
  entitySet locate_polygon(const vector<int>& nodes,
                           const store<vec2d>& pos) const ;
  // return all the data in the tree
  entitySet all_leafs() const
  {return root->leafs() ;}
  void insert_point(const vec2d& p, Entity identity)
  {root->insert_point(p,identity) ;}
} ;

/////////////////// PointSet_QuadTree class implementation ////////////////////
PointSet_QuadTree::qtnode::qtnode() {
  depth = 0 ;
  data_node = false ;
  point = vec2d(0.0, 0.0) ;
  id = EMPTY ;
  top_left = 0 ;
  top_right = 0 ;
  bot_left = 0 ;
  bot_right = 0 ;
}
PointSet_QuadTree::qtnode::qtnode(const rectangle& b) {
  depth = 0 ;
  data_node = false ;
  point = vec2d(0.0, 0.0) ;
  id = EMPTY ;
  top_left = top_right = bot_left = bot_right = 0 ;
  bbox = b ;
  set_children_bbox() ;
}
void PointSet_QuadTree::qtnode::erase_node() {
  if(top_left != 0) {
    if(top_left->data_node) {
      delete top_left ;
      top_left = 0 ;
    }
    else
      top_left->erase_node() ;
  }
  
  if(top_right != 0) {
    if(top_right->data_node) {
      delete top_right ;
      top_right = 0 ;
    }
    else
      top_right->erase_node() ;
  }

  if(bot_left != 0) {
    if(bot_left->data_node) {
      delete bot_left ;
      bot_left = 0 ;
    }
    else
      bot_left->erase_node() ;
  }

  if(bot_right != 0) {
    if(bot_right->data_node) {
      delete bot_right ;
      bot_right = 0 ;
    }
    else
      bot_right->erase_node() ;
  }
}
void PointSet_QuadTree::qtnode::copy_node(const qtnode& n) {
  erase_node() ;
  
  depth = n.depth ;
  data_node = n.data_node ;
  bbox = n.bbox ;
  point = n.point ;
  id = n.id ;
  tlb = n.tlb ;
  trb = n.trb ;
  blb = n.blb ;
  brb = n.brb ;
  
  if(!data_node) {
    top_left = new qtnode ;
    top_right = new qtnode ;
    bot_left = new qtnode ;
    bot_right = new qtnode ;
    top_left->copy_node(*(n.top_left)) ;
    top_right->copy_node(*(n.top_right)) ;
    bot_left->copy_node(*(n.bot_left)) ;
    bot_right->copy_node(*(n.bot_right)) ;
  }else {
    top_left = 0;
    top_right = 0 ;
    bot_left = 0 ;
    bot_right = 0 ;
  }
}
template<class STOREOFVEC2D>
bool is_allpoints_equal(const entitySet& allpoints,
                        const STOREOFVEC2D& pos) {
  entitySet::const_iterator ei = allpoints.begin() ;
  const vec2d& first = pos[*ei] ; ++ei ;
  for(;ei!=allpoints.end();++ei)
    if(first != pos[*ei])
      return false ;
  return true ;
}
template<class STOREOFVEC2D>
void PointSet_QuadTree::qtnode::build_node(const entitySet& allpoints,
                                           const STOREOFVEC2D& pos,
                                           int current_depth) {
  if(allpoints.size() == 0) {
    cerr << "ERROR: empty points in build_node! This should not occur!"
         << endl ;
    Loci::Abort() ;
  }
  depth = current_depth ;
  data_node = true ;
  if(allpoints.size() == 1) {
    Entity i = *(allpoints.begin()) ;
    id += i ;
    point = pos[i] ;
    return ;
  }

  if(is_allpoints_equal(allpoints,pos)) {
    id += allpoints ;
    point = pos[*(allpoints.begin())] ;
    return ;
  }

  data_node = false ;
  entitySet tl,tr,bl,br ;
  for(entitySet::const_iterator ei=allpoints.begin();
      ei!=allpoints.end();++ei) {
    const vec2d& p = pos[*ei] ;
    // test where this face falls into
    if(tlb.is_point_in(p))
      tl+=*ei ;
    if(trb.is_point_in(p))
      tr+=*ei ;
    if(blb.is_point_in(p))
      bl+=*ei ;
    if(brb.is_point_in(p))
      br+=*ei ;
  }
  // build subnodes
  if(tl.size() > 0) {
    top_left = new qtnode(tlb) ;
    top_left->build_node(tl,pos,depth+1) ;
  } else
    top_left = 0 ;
  
  if(tr.size() > 0) {
    top_right = new qtnode(trb) ;
    top_right->build_node(tr,pos,depth+1) ;
  } else
    top_right = 0 ;
  
  if(bl.size() > 0) {
    bot_left = new qtnode(blb) ;
    bot_left->build_node(bl,pos,depth+1) ;
  } else
    bot_left = 0 ;
  
  if(br.size() > 0) {
    bot_right = new qtnode(brb) ;
    bot_right->build_node(br,pos,depth+1) ;
  } else
    bot_right = 0 ;
}

int PointSet_QuadTree::qtnode::tree_depth() const {
  if(data_node)
    return depth ;
  else {
    int tld=0, trd=0, bld=0, brd=0 ;
    if(top_left != 0)
      tld = top_left->tree_depth() ;
    if(top_right != 0)
      trd = top_right->tree_depth() ;
    if(bot_left != 0)
      bld = bot_left->tree_depth() ;
    if(bot_right != 0)
      brd = bot_right->tree_depth() ;
    
    return max(max(tld,trd),max(bld,brd)) ;
  }
}

int PointSet_QuadTree::qtnode::leafs_num() const {
  if(data_node)
    return 1 ;
  else {
    int tle=0, tre=0, ble=0, bre=0 ;
    if(top_left != 0)
      tle = top_left->leafs_num() ;
    if(top_right != 0)
      tre = top_right->leafs_num() ;
    if(bot_left != 0)
      ble = bot_left->leafs_num() ;
    if(bot_right != 0)
      bre = bot_right->leafs_num() ;
    
    return tle+tre+ble+bre ;
  }
}

entitySet PointSet_QuadTree::qtnode::leafs() const {
  entitySet data ;
  if(data_node)
    data += id ;
  else {
    // recursively collect subnodes data
    if(top_left != 0)
      data += top_left->leafs() ;
    
    if(top_right != 0)
      data += top_right->leafs() ;
    
    if(bot_left != 0)
      data += bot_left->leafs() ;
    
    if(bot_right != 0)
      data += bot_right->leafs() ;
  }
  return data ;
}

entitySet PointSet_QuadTree::qtnode::
locate_rectangle(const rectangle& r) const {
  entitySet data ;
  if(data_node) {
    data += id ;
  } else {
    // splitting the passed in rectangle
    rectangle rsub ;
    if(tlb.rectangle_int(r,rsub))
      if(top_left != 0)
        data += top_left->locate_rectangle(rsub) ;

    if(trb.rectangle_int(r,rsub))
      if(top_right != 0)
        data += top_right->locate_rectangle(rsub) ;

    if(blb.rectangle_int(r,rsub))
      if(bot_left != 0)
        data += bot_left->locate_rectangle(rsub) ;

    if(brb.rectangle_int(r,rsub))
      if(bot_right != 0)
        data += bot_right->locate_rectangle(rsub) ;
  }
  return data ;
}

entitySet PointSet_QuadTree::qtnode::
locate_circle(const vec2d& center, double radius) const {
  entitySet data ;
  if(data_node) {
    data += id ;
  } else {
    if(tlb.is_circle_int(center,radius))
      if(top_left != 0)
        data += top_left->locate_circle(center,radius) ;

    if(trb.is_circle_int(center,radius))
      if(top_right != 0)
        data += top_right->locate_circle(center,radius) ;

    if(blb.is_circle_int(center,radius))
      if(bot_left != 0)
        data += bot_left->locate_circle(center,radius) ;

    if(brb.is_circle_int(center,radius))
      if(bot_right != 0)
        data += bot_right->locate_circle(center,radius) ;
  }
  return data ;
}

void PointSet_QuadTree::qtnode::
insert_point(const vec2d& p, Entity identity) {
  if(bbox.is_point_in(p)) {
    dstore<vec2d> pos ;
    pos[identity] = p ;
    entitySet ps ; ps += identity ;
    
    if(data_node) {
      // we need to rebuild the current quadrant
      for(entitySet::const_iterator ei=id.begin();
          ei!=id.end();++ei)
        pos[*ei] = point ;
      ps += id ;
      data_node = false ;
      this->build_node(ps,pos,depth) ;
    } else {

      if(tlb.is_point_in(p)) {
        if(top_left == 0) {
          top_left = new qtnode(tlb) ;
          top_left->build_node(ps,pos,depth+1) ;
        }else
          top_left->insert_point(p,identity) ;
      }
      if(trb.is_point_in(p)) {
        if(top_right == 0) {
          top_right = new qtnode(trb) ;
          top_right->build_node(ps,pos,depth+1) ;
        }else
          top_right->insert_point(p,identity) ;
      }
      if(blb.is_point_in(p)) {
        if(bot_left == 0) {
          bot_left = new qtnode(blb) ;
          bot_left->build_node(ps,pos,depth+1) ;
        }else
          bot_left->insert_point(p,identity) ;
      }
      if(brb.is_point_in(p)) {
        if(bot_right == 0) {
          bot_right = new qtnode(brb) ;
          bot_right->build_node(ps,pos,depth+1) ;
        }else
          bot_right->insert_point(p,identity) ;
      }
    }
  }
}

// we first compute a bounding box for the given polygon
// then call locate_rectangle to locate the bounding box
entitySet
PointSet_QuadTree::locate_polygon(const vector<int>& nodes,
                                  const store<vec2d>& pos) const {
  if(nodes.empty()) {
    cerr << "Error: empty polygon in locate_polygon!" << endl ;
    Loci::Abort() ;
  }
  // first we construct a store of position
//   vector<vec2d> nodes_pos ;
//   for(vector<int>::const_iterator ni=nodes.begin();
//       ni!=nodes.end();++ni)
//     nodes_pos.push_back(pos[*ni]) ;
//   rectangle b ;
//   nodes_bbox(nodes_pos,b) ;
  rectangle b ;
  nodes_bbox(nodes,pos,b) ;
  return locate_rectangle(b) ;
}
/////////////end of PointSet QuadTree implementation///////////////////


///////////////////////////////////////////////////////////////////////
// some other utility functions (structure preserving point moving)  //
///////////////////////////////////////////////////////////////////////
// this function computes the edge nodes set based on the relations
// this edge nodes is to the 2d case, i.e., all the nodes that lie
// on the edge (the outer contour line or boundary line) of the 2d grid
entitySet get_edge_nodes(const dMap& N1, const dMap& N2,
                         const dMap& Er) {
  entitySet edge_nodes ;
  entitySet edge_domain = Er.domain() ;
  for(entitySet::const_iterator ei=edge_domain.begin();
      ei!=edge_domain.end();++ei)
    if(Er[*ei] < 0) {
      edge_nodes += N1[*ei] ;
      edge_nodes += N2[*ei] ;
    }
  return edge_nodes ;
}

// return the closest point (id and distance)
// in a pointset to a given point
template<class STOREOFVEC2D>
void closest_point(const vec2d& p,
                   const entitySet& pset,
                   const STOREOFVEC2D& pset_pos,
                   int& id, double& d) {
  double min_d = std::numeric_limits<double>::max() ;
  id = -1 ;
  d = min_d ;
  for(entitySet::const_iterator ei=pset.begin();
      ei!=pset.end();++ei) {
    const vec2d& ppos = pset_pos[*ei] ;
    d = dist(p,ppos) ;
    if(d < min_d) {
      id = *ei ;
      min_d = d ;
    }
  }
  d = min_d ;
}

// given a point, compute the smallest length of
// all the edges that connect to it
double smallest_edge_length(int node_id,
                            const store<vec2d>& pos,
                            const multiMap& face2node,
                            const multiMap& node2face,
                            const entitySet& bc_faces) {
  // first get all the neighbor nodes
  entitySet faces = get_multiMap_elems(node2face,node_id) ;
  // get only the boundary faces
  faces &= bc_faces ;
  entitySet neighbor ;
  for(entitySet::const_iterator ei=faces.begin();
      ei!=faces.end();++ei)
    neighbor += get_multiMap_elems(face2node,*ei) ;
  neighbor -= node_id ;
  const vec2d& p = pos[node_id] ;
  double min_d = std::numeric_limits<double>::max() ;
  for(entitySet::const_iterator ei=neighbor.begin();
      ei!=neighbor.end();++ei) {
    const vec2d& ppos = pos[*ei] ;
    double d = dist(p,ppos) ;
    if(d < min_d)
      min_d = d ;
  }
  return min_d ;
}
// given a point, compute the smallest distance between
// the point and all the bounding edges
double smallest_bedge_dist(Entity node_id,
                           const store<vec2d>& pos,
                           const multiMap& face2node,
                           const multiMap& node2face,
                           const entitySet& bc_faces) {
  entitySet faces = get_multiMap_elems(node2face,node_id) ;
  // get only the boundary faces
  faces &= bc_faces ;

  const vec2d& p = pos[node_id] ;
  double min_d = std::numeric_limits<double>::max() ;

  for(entitySet::const_iterator ei=faces.begin();
      ei!=faces.end();++ei) {
    // here we need to keep the order of nodes
    vector<int> nodes ;
    int num_elems = face2node.num_elems(*ei) ;
    for(int i=0;i<num_elems;++i)
      nodes.push_back(face2node[*ei][i]) ;
    
    vector<int>::size_type i, j, size=nodes.size() ;
    for(i=0,j=size-1;i<size;j=i++) {
      // we don't compute the distance of the point
      // and the edges formed by the point itself
      if( (nodes[j] == node_id) || (nodes[i] == node_id))
        continue ;
      const vec2d& b = pos[nodes[j]] ;
      const vec2d& e = pos[nodes[i]] ;
      // we also need to test if point p lies collinear
      // to the segment, if it is, we don't compute the
      // distance. NOTE: this is possible because two
      // edges of a polygon could be collinear
      // and also note we use the approximate version
      // of the collinearity predicate. Because we
      // don't need an exact computing here. We use
      // the approximate version to help to eliminate
      // extremely small distance.
      if(collinear3_approx(p,b,e))
        continue ;
      double d = point2lineDist(p,b,e) ;
      if(d < min_d)
        min_d = d ;
    }
  }
  
  return min_d ;
}
// this function moves the points on BC2 to the closest
// one on BC1 while still preserving BC2's topological structure
// it returns all the points moved on BC2 and modifies
// corresponding BC2 nodes
// it also set the corresponding maps for corresponding
// identical points on both boundaries.
//
// This is an improved version that uses pointset quadTree
// (hence the name is called shift_points2)
// It should be faster and more reliable
// The original code is cleaned
// (but we remain the name for historical reason)
entitySet shift_points2(const entitySet& bc1_faces,
                        const entitySet& bc2_faces,
                        const entitySet& bc1_nodes,
                        const store<vec2d>& bc1_npos,
                        const entitySet& bc2_nodes,
                        store<vec2d>& bc2_npos,
                        const multiMap& face2node,
                        const multiMap& node2face,
                        dMap& bc1p_2_bc2p,
                        dMap& bc2p_2_bc1p) {
#ifdef POINT_SHIFT_CONVEXITY_CHECK
  {
    // check for bc2 face convexity
    int non_convex = 0 ;
    for(entitySet::const_iterator
          fi=bc2_faces.begin();fi!=bc2_faces.end();++fi) {
      // first get face nodes
      vector<int> nodes ;
      vector<vec2d> nodes_pos ;
      int bc_nodes = 0 ;
      int num_elems = face2node.num_elems(*fi) ;
      for(int i=0;i<num_elems;++i) {
        int n = face2node[*fi][i] ;
        if(slavebc_edge_nodes.inSet(n))
          ++bc_nodes ;
        nodes.push_back(n) ;
        nodes_pos.push_back(bc2_npos[n]) ;
      }
      if(bc_nodes >= 3) {
        if(!check_convex_approx(nodes_pos))
          ++non_convex ;
      } else {
        if(!check_convex(nodes_pos))
          ++non_convex ;
      }
    }
    cerr << "Before SHIFT, total NON convex polygon: "
         << non_convex << endl ;
  }
#endif
  
  entitySet shifted ;
  // Here we do point shifting, that is we choose BC1 as
  // the master boundary, and we shift BC2 points to match
  // BC1 points in the given range.

  // First we construct a PointSet QuadTree for BC1 nodes
  gettimeofday(&time_shift_point_qt_start,NULL) ;
  rectangle bounding_box ;
  nodes_bbox(bc1_npos,bounding_box) ;
  // make a slightly larger bounding box
  bounding_box.l -= vec2d(0.01,-0.01) ;
  bounding_box.r += vec2d(0.01,-0.01) ;
  PointSet_QuadTree qtree(bounding_box) ;
  qtree.build_tree(bc1_nodes,bc1_npos) ;
  gettimeofday(&time_shift_point_qt_end,NULL) ;

  gettimeofday(&time_shift_point_start,NULL) ;
  // then we shift bc2 points
  for(entitySet::const_iterator ei=bc2_nodes.begin();
      ei!=bc2_nodes.end();++ei) {
    // first we need to compute the shifting threshold
    // for this bc2 node. It is the smallest distance
    // to its bounding polygon
    double threshold =
      smallest_bedge_dist(*ei,bc2_npos,face2node,node2face,bc2_faces) ;
    threshold *= SHIFT_RANGE_THRESHOLD ;

    // then we search in the bc1 pointset quadtree and see
    // if there are any bc1 points lie in the range
    entitySet neighbor
      = qtree.locate_circle(bc2_npos[*ei],threshold) ;

    if(neighbor == EMPTY) // we don't have points in range
      continue ;
    
    int p_id ; double min_d ;
    closest_point(bc2_npos[*ei],neighbor,bc1_npos,p_id,min_d) ;
    
    if(min_d <= threshold) {
      // then this point can be shifted
      shifted += *ei ;
      // edit position
      bc2_npos[*ei] = bc1_npos[p_id] ;
      // and we also set up the relation map
      bc1p_2_bc2p[p_id] = *ei ;
      bc2p_2_bc1p[*ei] = p_id ;
    }
  }
  gettimeofday(&time_shift_point_end,NULL) ;

#ifdef POINT_SHIFT_CONVEXITY_CHECK
  {
    // check for bc2 face convexity
    int non_convex = 0 ;
    for(entitySet::const_iterator
          fi=bc2_faces.begin();fi!=bc2_faces.end();++fi) {
      // first get face nodes
      vector<int> nodes ;
      vector<vec2d> nodes_pos ;
      int bc_nodes = 0 ;
      int num_elems = face2node.num_elems(*fi) ;
      for(int i=0;i<num_elems;++i) {
        int n = face2node[*fi][i] ;
        if(slavebc_edge_nodes.inSet(n))
          ++bc_nodes ;
        nodes.push_back(n) ;
        nodes_pos.push_back(bc2_npos[n]) ;
      }
      if(bc_nodes >= 3) {
        if(!check_convex_approx(nodes_pos))
          ++non_convex ;
      } else {
        if(!check_convex(nodes_pos))
          ++non_convex ;
      }
    }
    cerr << "SHIFT done, total Non convex polygon: "
         << non_convex << endl ;
  }
#endif

  return shifted ;
}

/**
***********************************************************
*  face splitting routines, split bc1 faces according to  *
*  the overlap of the two boundaries.                     *
***********************************************************
**/
template<typename T1, typename T2>
std::ostream& print_map(const map<T1,T2>& m, std::ostream& s) {
  for(typename map<T1,T2>::const_iterator mi=m.begin();
      mi!=m.end();++mi)
    s << mi->first << " --> " << mi->second << endl ;
  return s ;
}

// return removable faces for bc1 (face that do not
// need to be splitted
entitySet
bc1_removable_faces(const entitySet& bc1_faces,
                    const entitySet& bc2_faces,
                    const multiMap& face2node,
                    const multiMap& node2face,
                    const dMap& bc1p_2_bc2p) {
  entitySet removable ;
  
  entitySet bc1_matched = bc1p_2_bc2p.domain() ;
  for(entitySet::const_iterator ei=bc1_faces.begin();
      ei!=bc1_faces.end();++ei) {
    // first get all the face nodes
    entitySet ns = get_multiMap_elems(face2node,*ei) ;
    // see if all the face nodes are matched
    if( (ns - bc1_matched) != EMPTY)
      continue ;
    // now we test are the corresponding matched points
    // on bc2 also form a valid bc2 face? If it is so,
    // then this bc1 face are matched with a bc2 face
    // and can be removed from splitting
    entitySet bc2_correspond = bc1p_2_bc2p.image(ns) ;
    entitySet bc2_face_cand = interval(Loci::UNIVERSE_MIN,
                                       Loci::UNIVERSE_MAX) ;
    for(entitySet::const_iterator fi=bc2_correspond.begin();
        fi!=bc2_correspond.end();++fi)
      bc2_face_cand &= get_multiMap_elems(node2face,*fi) ;
    bc2_face_cand &= bc2_faces ;
    if(bc2_face_cand.size() == 1) {
      // then we'll need to compare the face nodes number
      // they should be equal if two faces are equal
      entitySet bc2_face_nodes =
        get_multiMap_elems(face2node,*(bc2_face_cand.begin())) ;
      if(ns.size() == bc2_face_nodes.size())
        removable += *ei ;
    }
  }
  
  return removable ;
}

// this function converts a dstore to store
template<typename T>
void dstore2store(const dstore<T>& ds, store<T>& s) {
  // clear the store first
  s.allocate(EMPTY) ;
  entitySet ds_dom = ds.domain() ;
  s.allocate(ds_dom) ;
  for(entitySet::const_iterator ei=ds_dom.begin();
      ei!=ds_dom.end();++ei)
    s[*ei] = ds[*ei] ;
}
// converts an STL map to Loci multiMap
void map2multiMap(const map<Entity,vector<Entity> >& stlm,
                  multiMap& locim) {
  // clear the multiMap first
  locim.allocate(EMPTY) ;
  entitySet count_domain ;
  for(map<Entity,vector<Entity> >::const_iterator mi=stlm.begin();
      mi!=stlm.end();++mi)
    count_domain += mi->first ;
  store<int> count ;
  count.allocate(count_domain) ;
  for(map<Entity,vector<Entity> >::const_iterator mi=stlm.begin();
      mi!=stlm.end();++mi)
    count[mi->first] = (mi->second).size() ;

  locim.allocate(count) ;
  for(map<Entity,vector<Entity> >::const_iterator mi=stlm.begin();
      mi!=stlm.end();++mi) {
    const vector<Entity> ve = mi->second ;
    for(vector<Entity>::size_type i=0;i!=ve.size();++i)
      locim[mi->first][i] = ve[i] ;
  }
}
// converts an STL map to Loci multiMap
void map2multiMap(const map<Entity,entitySet>& stlm,
                  multiMap& locim) {
  // clear the multiMap first
  locim.allocate(EMPTY) ;
  entitySet count_domain ;
  for(map<Entity,entitySet>::const_iterator mi=stlm.begin();
      mi!=stlm.end();++mi)
    count_domain += mi->first ;
  store<int> count ;
  count.allocate(count_domain) ;
  for(map<Entity,entitySet>::const_iterator mi=stlm.begin();
      mi!=stlm.end();++mi)
    count[mi->first] = (mi->second).size() ;

  locim.allocate(count) ;
  for(map<Entity,entitySet>::const_iterator mi=stlm.begin();
      mi!=stlm.end();++mi) {
    const entitySet es = mi->second ;
    int i=0 ;
    for(entitySet::const_iterator ei=es.begin();ei!=es.end();++ei,++i)
      locim[mi->first][i] = *ei ;
  }
}

void print_face_node_pos(const vector<int>& nodes,
                         const store<vec2d>& pos) {
  for(vector<int>::const_iterator vi=nodes.begin();
      vi!=nodes.end();++vi)
    cout << pos[*vi] << endl ;
  cout << "++++++++" << endl ;
}

// this function returns the information that needed
// for inverse projection of 2d points back to 3d space
struct alphaBeta {
  // P = p0 + a(p1 - p0) + b(p2 - p0)
  Entity p0, p1, p2 ;
  double a, b ;
  alphaBeta():p0(-1),p1(-1),p2(-1),a(0.0),b(0.0) {}
  alphaBeta(Entity i, Entity j, Entity k, double d1, double d2):
    p0(i),p1(j),p2(k),a(d1),b(d2) {}
} ;
namespace Loci {
  template<>
  struct  data_schema_traits<alphaBeta> {
  public:
    typedef IDENTITY_CONVERTER Schema_Converter;
    static DatatypeP get_type() {
      alphaBeta t ;
      CompoundDatatypeP ct = CompoundFactory(t) ;
      LOCI_INSERT_TYPE(ct,alphaBeta,p0) ;
      LOCI_INSERT_TYPE(ct,alphaBeta,p1) ;
      LOCI_INSERT_TYPE(ct,alphaBeta,p2) ;
      LOCI_INSERT_TYPE(ct,alphaBeta,a) ;
      LOCI_INSERT_TYPE(ct,alphaBeta,b) ;
      return DatatypeP(ct) ;
    }
  };
}

// triangulate a given convex polygon
// it is assumed that the passed in nodes
// define a convex polygon in either clockwise
// or counter-clockwise order
//
// This function will output a triangulation that includes
// every node in the polygon. i.e., the triangulation will
// not skip any node in the polygon, even if the original
// polygon contains several nodes that are collinear.
// Of course, the triangulation will not contain any
// degenerated triangles formed by collinear nodes.
//
// We first pick a triangle that is not degenerated (p0, p1, p2).
// Then we use the p0 as the fixed point and move p1, p2 to form
// new triangles until we detect a degenerated triangle (if any).
// Then we use the current p1 as the starting point p0 and moving
// p1 and p2 in the remaining points to form triangles.
// This is a linear time algorithm. 
vector<vector<Entity> >
triangulate_convex_polygon(const vector<Entity>& poly_nodes,
                           const store<vec2d>& pos) {
  vector<vector<Entity> > triangulation ;
  vector<Entity>::size_type size = poly_nodes.size() ;
  if(size < 3) {
    cout << "ERROR: polygon ill defined in triangulate_convex_polygon!"
         << endl ;
    Loci::Abort() ;
  }
  if(size == 3) {
    // trivial case
    triangulation.push_back(poly_nodes) ;
    return triangulation ;
  }
  vector<Entity>::size_type p0_index = 0 ;
  vector<Entity>::size_type p1_index = 2 ;
  vector<Entity>::size_type p2_index = 1 ;
  bool polygon_degenerated = false ;
  vector<Entity>::size_type start_p0 = p0_index ;
  vector<Entity>::size_type start_p1 = p1_index ;
  vector<Entity>::size_type start_p2 = p2_index ;
  // we first need to find a non-degenerated triangle
  while(true) {
    Entity p0 = poly_nodes[p0_index] ;
    Entity p1 = poly_nodes[p1_index] ;
    Entity p2 = poly_nodes[p2_index] ;
    const vec2d& p0_pos = pos[p0] ;
    const vec2d& p1_pos = pos[p1] ;
    const vec2d& p2_pos = pos[p2] ;
    // we need to confirm that this triangle is
    // not degenerated to an edge
    if(!collinear3_approx(p0_pos,p1_pos,p2_pos)) {
      start_p0 = p0_index ;
      start_p1 = p1_index ;
      start_p2 = p2_index ;
      break ;
    }else {
      p0_index = (p0_index+1)%size ;
      p1_index = (p1_index+1)%size ;
      p2_index = (p2_index+1)%size ;
      if(p0_index == 0) {
        // we looped back to the start, hence the
        // polygon itself is degenerated and this is
        // an error
        polygon_degenerated = true ;
        break ;
      }
    }
  }
  if(polygon_degenerated) {
    cout << "ERROR: passed in polygon degenerated in "
         << "triangulate_convex_polygon!" << endl ;
    Loci::Abort() ;
  }

  // now we triangle the first part of the polygon
  p0_index = start_p0 ;
  p1_index = start_p1 ;
  p2_index = start_p2 ;
  bool triangulation_complete = false ;
  while(true) {
    vector<Entity> triangle ;
    Entity p0 = poly_nodes[p0_index] ;
    Entity p1 = poly_nodes[p1_index] ;
    Entity p2 = poly_nodes[p2_index] ;
    const vec2d& p0_pos = pos[p0] ;
    const vec2d& p1_pos = pos[p1] ;
    const vec2d& p2_pos = pos[p2] ;
    if(!collinear3_approx(p0_pos,p1_pos,p2_pos)) {
      triangle.push_back(p0) ;
      triangle.push_back(p1) ;
      triangle.push_back(p2) ;
      triangulation.push_back(triangle) ;
      p1_index = (p1_index+1)%size ;
      p2_index = (p2_index+1)%size ;
      if(p1_index == start_p0) {
        // by now, we've completed the triangulation
        // and did not found any degenerated triangles
        triangulation_complete = true ;
        break ;
      }
    }else {
      break ;
    }
  }
  if(!triangulation_complete) {
    // first we need to remove the last triangle from
    // triangulation because it is not a valid one
    triangulation.pop_back() ;
    // then setup the new starting node p0
    p0_index = (p2_index-1+size)%size ;
    p1_index = (p0_index+2)%size ;
    p2_index = (p0_index+1)%size ;
    while(p1_index != start_p2) {
      vector<Entity> triangle ;
      Entity p0 = poly_nodes[p0_index] ;
      Entity p1 = poly_nodes[p1_index] ;
      Entity p2 = poly_nodes[p2_index] ;
      // we don't need to test for degenerated triangle
      // here because there will be no such triangles
      // from now on.
      triangle.push_back(p0) ;
      triangle.push_back(p1) ;
      triangle.push_back(p2) ;
      triangulation.push_back(triangle) ;
      p1_index = (p1_index+1)%size ;
      p2_index = (p2_index+1)%size ;
    }
  }

  return triangulation ;
}

alphaBeta get_inverse_proj_info(const vector<vector<Entity> >& tri,
                                const store<vec2d>& pos,
                                const vec2d& p, Entity pid) {
  // we need to first determine where the point p lies
  // in the convex polygon (i.e., in which triangle)
  bool found = false ;
  Entity p0=-1,p1=-1,p2=-1 ;
  vector<vector<Entity> >::const_iterator vi ;
  for(vi=tri.begin();vi!=tri.end();++vi) {
    p0 = (*vi)[0] ;
    p1 = (*vi)[1] ;
    p2 = (*vi)[2] ;
    // We first query the table if p is on any edge
    // of the triangle (p0,p1,p2) because p might be
    // computed inexactly
    if(global_pbetween_ltable.has_record(pid,p0,p1)) {
      found = true ;
      break ;
    }    
    if(global_pbetween_ltable.has_record(pid,p1,p2)) {
      found = true ;
      break ;
    }
    if(global_pbetween_ltable.has_record(pid,p0,p2)) {
      found = true ;
      break ;
    }
    vec2d p0_pos = pos[p0] ;
    vec2d p1_pos = pos[p1] ;
    vec2d p2_pos = pos[p2] ;
    // We then test if p is inside the triangle or not
    if(!left_exact(p2_pos,p1_pos,p0_pos)) {
      swap(p2,p1) ;
      swap(p2_pos,p1_pos) ;
    }
    if(leftOn_exact(p0_pos,p2_pos,p) &&
       leftOn_exact(p1_pos,p0_pos,p) &&
       leftOn_exact(p2_pos,p1_pos,p) ) {
      found = true ;
      break ;
    }
    // If it failed by now, it is also possible
    // that this node is on boundary edge so that
    // it value is not precise and its relation
    // to the triangle (p0,p1,p2) has not been
    // recorded. So we give it a last chance
    // using approximate predicate testing
    bool p0_bn = false ; bool p1_bn = false ; bool p2_bn = false ;
    if(!slavebc_edge_nodes.inSet(pid))
      continue ;
    if(masterbc_edge_nodes.inSet(p0))
      p0_bn = true ;
    if(masterbc_edge_nodes.inSet(p1))
      p1_bn = true ;
    if(masterbc_edge_nodes.inSet(p2))
      p2_bn = true ;

    if(p0_bn && p2_bn)
      if(pointOnEdge_approx(p,p0_pos,p2_pos)) {
        global_pbetween_ltable.add_record(pid,p0,p2) ;
        found = true ;
        break ;
      }
    if(p1_bn && p0_bn)
      if(pointOnEdge_approx(p,p1_pos,p0_pos)) {
        global_pbetween_ltable.add_record(pid,p1,p2) ;
        found = true ;
        break ;
      }
    if(p2_bn && p1_bn)
      if(pointOnEdge_approx(p,p2_pos,p1_pos)) {
        global_pbetween_ltable.add_record(pid,p2,p1) ;
        found = true ;
        break ;
      }
  }
  if(!found) {
    cerr << "ERROR: in triangle testing failed in computing alphaBeta info!"
         << endl ;
    /////////////////////
    cerr.precision(32) ;
    cerr << "p: " << p << endl ;
    cerr << "pid: " << pid << endl ;
    cerr << "p0: " << p0 << " p1: " << p1 << " p2: " << p2 << endl ;
    cerr << "tri.size: " << tri.size() << endl ;
    cerr << "p0: " << pos[p0] << endl ;
    cerr << "p1: " << pos[p1] << endl ;
    cerr << "p2: " << pos[p2] << endl ;
    /////////////////////
    Loci::Abort() ;
  }
  // now we know the triangle that contains p is(p0,p1,p2) ;
  vec2d v1 = pos[p1] - pos[p0] ;
  vec2d v2 = pos[p2] - pos[p0] ;
  vec2d P = p - pos[p0] ;
  double a = cross(P,v2) / cross(v1,v2) ;
  double b = cross(P,v1) / cross(v2,v1) ;

  return alphaBeta(p0,p1,p2,a,b) ;
}

// struct used for sorting added nodes on boundary
struct bc_node_info {
  Entity id ;
  double ab_value ; // the alpha or beta value
  bc_node_info(Entity i, double a):id(i),ab_value(a) {}
} ;
// comparison functions
inline bool ascend_ab(const bc_node_info& n1, const bc_node_info& n2)
{return n1.ab_value < n2.ab_value ;}

// structure for computing new face2node relation for
// the face on the other boundary that share the boundary
// edges with BC1
struct new_interior_face_info {
  bool valid ;
  Entity p0, p1 ;
  entitySet faces ;
  vector<bc_node_info> new_nodes ;
  new_interior_face_info():valid(false),p0(-1),p1(-1) {}
} ;
namespace Loci {
  // this is actually not correctly defined for Loci,
  // but we just want the program to be able to compile
  // so we are sloppy here
  template<>
  struct  data_schema_traits<new_interior_face_info> {
  public:
    typedef IDENTITY_CONVERTER Schema_Converter;
    static DatatypeP get_type() {
      new_interior_face_info t ;
      CompoundDatatypeP ct = CompoundFactory(t) ;
      LOCI_INSERT_TYPE(ct,new_interior_face_info,valid) ;
      LOCI_INSERT_TYPE(ct,new_interior_face_info,p0) ;
      LOCI_INSERT_TYPE(ct,new_interior_face_info,p1) ;
      //LOCI_INSERT_TYPE(ct,new_interior_face_info,faces) ;
      //LOCI_INSERT_TYPE(ct,new_interior_face_info,new_nodes) ;
      return DatatypeP(ct) ;
    }
  };
}

// this function checks for any points added to the interior face
void
check_interior_node_insertion(const alphaBeta& abinfo,
                              const multiMap& node2edge,
                              const multiMap& edge2iface,
                              Entity new_node_id,
                              Entity real_id,
                              dstore<new_interior_face_info>& nifi) {
  Entity p0=0, p1=0 ;
  double ab_value=0.0 ;
  // We first need to check whether the new node lies on
  // an edge of the triangle (p0,p1,p2) ;
  // first we check if it lies on (p0,p2)
  if(global_pbetween_ltable.has_record(new_node_id,
                                       abinfo.p0, abinfo.p2)) {
    // then the point lies on (p0 p2)
    p0 = abinfo.p0 ;
    p1 = abinfo.p2 ;
    ab_value = abinfo.b ;
  } else if(global_pbetween_ltable.has_record(new_node_id,
                                              abinfo.p0, abinfo.p1)) {
    // then the point lies on (p0,p1)
    p0 = abinfo.p0 ;
    p1 = abinfo.p1 ;
    ab_value = abinfo.a ;
  } else if(global_pbetween_ltable.has_record(new_node_id,
                                              abinfo.p1, abinfo.p2)) {
    // then the point lies on (p1 p2)
    p0 = abinfo.p1 ;
    p1 = abinfo.p2 ;
    ab_value = abinfo.b ;
  } else
    return ;

  // then get the edge defined by p0 and p1
  entitySet n1_edges = get_multiMap_elems(node2edge,p0) ;
  entitySet n2_edges = get_multiMap_elems(node2edge,p1) ;
  entitySet eint = n1_edges & n2_edges ;
  if(eint.size() == 0)
    return ; // this is not an original edge in the mesh
  if(eint.size() != 1) {
    cerr << "ERROR: data inconsistent in check_bc_node_insertion!"
         << endl ;
    Loci::Abort() ;
  }
  Entity e = *(eint.begin()) ;

  if(!nifi.domain().inSet(e))
    nifi[e] = new_interior_face_info() ;

  new_interior_face_info& nifi_elem = nifi[e] ;

  if(!nifi_elem.valid) {
    nifi_elem.valid = true ;
    nifi_elem.p0 = p0 ;
    nifi_elem.p1 = p1 ;
    nifi_elem.faces = get_multiMap_elems(edge2iface,e) ;
  }else {
    if( (p0 != nifi_elem.p0) || (p1 != nifi_elem.p1)) {
      cerr << "ERROR: data inconsistent in check_interior_node_insertion!"
           << " edge does not agree!" << endl ;
      Loci::Abort() ;
    }
    entitySet faces = get_multiMap_elems(edge2iface,e) ;
    if(faces != nifi_elem.faces) {
      cerr << "ERROR: data inconsistent in check_interior_node_insertion!"
           << " face does not agree!" << endl ;
      Loci::Abort() ;
    }
  }

  nifi_elem.new_nodes.push_back(bc_node_info(real_id,ab_value)) ;
}
// this function reorders the passed in vector so that
// the first two elements in the vector are the same
// as p0 and p1. e.g., p0=79,p1=78,v=78,171,79
// then v will be reordered to: v=79,78,171
// e.g., p0=66,p1=65,v=66,152,65 then v will be
// reordered to: v=65,66,152 and p0 and p1 will be swapped
// NOTE: the reordering of vector does not change the relative
// order of elements in it, it merely just is a rotation
// of the vector
// This function returns if p0 and p1 will need to be swapped,
// but it does not really swap p0 and p1. As this flag is
// just an indication of a reverse of the new_nodes vector
// in the new_interior_face_info struct
bool reorder_face_node(vector<Entity>& v,
                       Entity p0, Entity p1) {
  if(v.size() < 3) {
    cerr << "ERROR: vector size less than three in reorder_face_node!"
         << endl ;
    Loci::Abort() ;
  }
  vector<Entity>::size_type i, size = v.size() ;
  bool match = false ;
  bool p0p1swap = false ;
  for(i=0;i<size;++i) {
    if( (v[0] == p0)&&(v[1] == p1)) {
      match = true ;
      break ;
    }
    vector<Entity>::iterator next = v.begin()+1 ;
    rotate(v.begin(),next,v.end()) ;
  }
  if(!match) {
    p0p1swap = true ;
    swap(p0,p1) ;
    for(i=0;i<size;++i) {
      if( (v[0] == p0)&&(v[1] == p1)) {
        match = true ;
        break ;
      }
      vector<Entity>::iterator next = v.begin()+1 ;
      rotate(v.begin(),next,v.end()) ;
    }
  }
  if(!match) {
    // If we are here, that means there are something in between
    // p0 and p1 in the vector. This is not possible if the polygon
    // splitting and new nodes lookup worked correctly. This error
    // is irrecoverable and we should give warning and stop the
    // program.
    cerr << "ERROR: data inconsistent in reorder_face_node!" << endl ;
    Loci::Abort() ;
  }
  return p0p1swap ;
}
// get new face2node relation for interior faces
void
add_new_interior_face2node(dstore<new_interior_face_info>& nifi,
                           const multiMap& face2node,
                           map<Entity,vector<Entity> >& m) {
  entitySet edges = nifi.domain() ;
  for(entitySet::const_iterator ei=edges.begin();ei!=edges.end();++ei) {
    new_interior_face_info& nifi_elem = nifi[*ei] ;
    // first sort the nifi_elem new nodes
    sort(nifi_elem.new_nodes.begin(),nifi_elem.new_nodes.end(),ascend_ab) ;

    entitySet interior_faces = nifi_elem.faces ;
    for(entitySet::const_iterator fi=interior_faces.begin();
        fi!=interior_faces.end();++fi) {
      vector<Entity> orig_face_nodes ;
      // here we need to get the original face definition.
      // it is possible during the face splitting, the face
      // definition has already been changed. Hence we need
      // to first query the STL map to see whether the face
      // definition is there. If it is not there in the STL map,
      // then we retrieve the face definition from the original
      // face2node multiMap.
      map<Entity,vector<Entity> >::const_iterator mi = m.find(*fi) ;
      if(mi != m.end())
        orig_face_nodes = mi->second ;
      else
        orig_face_nodes = get_face_nodes(face2node,*fi) ;
      
      Entity p0 = nifi_elem.p0, p1 = nifi_elem.p1 ;
      bool reverse = false ;
      if(reorder_face_node(orig_face_nodes, p0, p1)) {
        reverse = true ;
        swap(p0,p1) ;
      }

      vector<Entity> new_face_nodes ;
      new_face_nodes.push_back(p0) ;
      if(reverse) {
        // push back in the reverse order (from the end -> begin)
        for(vector<bc_node_info>::reverse_iterator
              vi=(nifi_elem.new_nodes).rbegin();
            vi!=(nifi_elem.new_nodes).rend();++vi)
          new_face_nodes.push_back(vi->id) ;
      }else {
        for(vector<bc_node_info>::const_iterator
              vi=(nifi_elem.new_nodes).begin();
            vi!=(nifi_elem.new_nodes).end();++vi)
          new_face_nodes.push_back(vi->id) ;
      }
      new_face_nodes.push_back(p1) ;  
  
      for(vector<Entity>::size_type i=2;i<orig_face_nodes.size();++i)
        new_face_nodes.push_back(orig_face_nodes[i]) ;

      // finally we put this in record
      m[*fi] = new_face_nodes ;
    } // end of for (interior face)
  } // end of for (edges)
}

// utility function for debugging purpose
// visualize the given faceset in the "showme" program
void write_out_face_poly(const entitySet& faces,
                         const multiMap& face2node,
                         const store<vec2d>& pos,
                         const char* fname) {
  // first computes node number
  int node_number = 0 ;
  for(entitySet::const_iterator ei=faces.begin();
      ei!=faces.end();++ei)
    node_number += face2node.num_elems(*ei) ;
  
  // output to poly file
  std::ofstream out(fname,ios::out) ;
  out.precision(16) ;
  // output nodes
  out <<  node_number << ' ' << "2 1 0" << endl ;
  int cnt = 1 ;
  for(entitySet::const_iterator ei=faces.begin();
      ei!=faces.end();++ei) {
    const vector<Entity>& vn = get_face_nodes(face2node,*ei) ;
    for(vector<Entity>::const_iterator vi=vn.begin();
        vi!=vn.end();++vi,++cnt)
      out << cnt << ' ' << pos[*vi].x << ' ' << pos[*vi].y << endl ;
    out << endl ;
  }
  // output edges
  out << node_number << " 0" << endl ;
  cnt = 1 ;
  
  for(entitySet::const_iterator ei=faces.begin();
      ei!=faces.end();++ei) {
    int f_size = face2node.num_elems(*ei) ;
    int i ;
    int begin = cnt ;
    for(i=0;i<f_size-1;++i,++cnt) {
      out << cnt << ' ' << i+begin << ' ' << i+begin+1 << endl ;
    }
    out << cnt++ << ' ' << i+begin << ' ' <<  begin << ' ' << endl ;
    out << endl ;
  }
  // 0 holes
  out << "0" << endl ;
  out.close() ;
}

void write_out_face_poly(const vector<vector<vec2d> >& mesh,
                         const char* poly_name) {
  // computes node number
  int all_nodes_number = 0 ;
  for(vector<vector<vec2d> >::const_iterator vi=mesh.begin();
      vi!=mesh.end();++vi)
    all_nodes_number += vi->size() ;
  
  // output to poly file
  std::ofstream out(poly_name,ios::out) ;
  out.precision(16) ;
  // output nodes
  out <<  all_nodes_number << ' ' << "2 1 0" << endl ;
  int cnt = 1 ;
  for(vector<vector<vec2d> >::const_iterator vi=mesh.begin();
      vi!=mesh.end();++vi) {
    for(vector<vec2d>::const_iterator vi2=vi->begin();
        vi2!=vi->end();++vi2,++cnt)
      out << cnt << ' ' << vi2->x << ' ' << vi2->y << endl ;
    out << endl ;
  }
  // output edges
  out << all_nodes_number << " 0" << endl ;
  cnt = 1 ;
  typedef vector<vec2d>::size_type size ;

  for(vector<vector<vec2d> >::const_iterator vi=mesh.begin();
      vi!=mesh.end();++vi) {
    size p_size = vi->size() ;
    size i ;
    size begin = cnt ;
    for(i=0;i<p_size-1;++i,++cnt) {
      out << cnt << ' ' << i+begin << ' ' << i+begin+1 << endl ;
    }
    out << cnt++ << ' ' << i+begin << ' ' <<  begin << ' ' << endl ;
    out << endl ;
  }
  // 0 holes
  out << "0" << endl ;
  out.close() ;
}

// a small object for vec2d comparison (operator<)
struct vec2dLess {
  bool operator()(const vec2d& v1, const vec2d& v2) const {
    if(v1.x == v2.x)
      return v1.y < v2.y ;
    else
      return v1.x < v2.x ;
  }
} ;
// This function remaps the polygon intersection index
// to the global new node numbering
vector<Entity> remap_intersection_index
(const vector<Entity>& index,
 const entitySet& master_nodes,
 const entitySet& slave_nodes,
 const vector<Entity>& master_poly,
 const store<vec2d>& master_pos,
 const store<vec2d>& slave_pos) {

  entitySet Q_remap_dom = global_Q_id_remap.domain() ;
  // first we construct a map for the original
  // polygon on the master boundary for later lookup
  map<vec2d,Entity,vec2dLess> pmap ;
  for(vector<Entity>::const_iterator vi=master_poly.begin();
      vi!=master_poly.end();++vi)
    pmap[master_pos[*vi]] = *vi ;
  // then we generate new index
  // we will generate a new index for any node that
  // is NOT already on the master boundary
  vector<Entity> real_index ;
  map<vec2d,Entity,vec2dLess>::const_iterator mi ;
  for(vector<Entity>::const_iterator vi=index.begin();
      vi!=index.end();++vi) {
    if(master_nodes.inSet(*vi)) {
      real_index.push_back(*vi) ;
    } else if(Q_remap_dom.inSet(*vi)) {
      // first we query the Q remap,
      // the Q remap is a map used to record
      // node id maps for nodes in the intersection polygon
      // that are NOT on the master boundary.
      real_index.push_back(global_Q_id_remap[*vi]) ;
    } else if(slave_nodes.inSet(*vi) &&
              (mi=pmap.find(slave_pos[*vi])) != pmap.end()) {
      // if not in and originally a node on the slave boundary,
      // and the node overlaps with one of the master boundary
      // node, then we'd like to use the master node id
      real_index.push_back(mi->second) ;
      global_Q_id_remap[*vi] = mi->second ;
    } else {
      // otherwise, generate a new index
      real_index.push_back(mesh_new_node_index) ;
      global_Q_id_remap[*vi] = mesh_new_node_index ;
      ++mesh_new_node_index ;
    }
  }

  return real_index ;
}

// this is a small struct for recording face splitting results
struct split_unit {
  // each split consists of the nodes and their ids
  vector<vec2d> nodes ;
  vector<Entity> ids ;
  split_unit():nodes(vector<vec2d>()),ids(vector<Entity>()) {}
  split_unit(const vector<vec2d>& n, const vector<Entity>& i)
    :nodes(n),ids(i) {}
} ;
// a face split consists of a vector of such units
map<Entity, vector<split_unit> > global_bc1_split_records ;
map<Entity, vector<split_unit> > global_bc2_split_records ;

// a utility function to write out a collection of polygons
// in the format of the 2d viz program "showme"
void polygons2showme(const string& fname,
                     const vector<vector<vec2d> >& polys) {
  ofstream output(fname.c_str(), std::ios::out) ;
  output.precision(16) ;
    
  int total_nodes = 0 ;
  for(size_t i=0;i!=polys.size();++i)
    total_nodes += polys[i].size() ;

  // header
  output << total_nodes << " 2 0 0" << endl << endl ;
  // polygon node position
  int index = 0 ;
  for(size_t i=0;i!=polys.size();++i) {
    for(size_t j=0;j!=polys[i].size();++j,++index) {
      output << index << " " << polys[i][j] << endl ;
    }
    output << endl ;
  }
  output << total_nodes << " 0" << endl << endl ;
  // polygon segment definition
  index = 0 ;
  for(size_t i=0;i!=polys.size();++i) {
    int start = index ;
    for(size_t j=0;j!=polys[i].size();++j,++index) {
      if(j == polys[i].size() - 1)
        output << index << " " << index << " " << start << endl ;
      else
        output << index << " " << index  << " " << index+1 << endl ;
    }
    output << endl ;
  }
  // no holes
  output << "0" << endl ;
  output.close() ;
}

// a small function that creates a string represents a status bar
// i.e., [.....             ],
// passed in s is the total points, n is the displayed points
string
get_status_bar(size_t s, size_t n) {
  string bar = "[" ;
  bar.append(string(n,'>')) ;
  bar.append(string(s-n,'_')) ;
  bar.append("]") ;
  return bar ;
}

#ifdef FACE_SPLIT_SELF_CHECK
inline bool eq_vec2d(const vec2d& v1, const vec2d& v2) {
  return v1==v2 ;
}
#endif
 
// the main face splitting function,
// this function only creates the global_face_split_records,
// later functions will build up new mesh topology from that
// record
// returns the total number of pairs of faces that are split
int face_split(const entitySet& bc1_faces,
               const entitySet& bc2_faces,
               const entitySet& bc1_nodes,
               const store<vec2d>& bc1_npos,
               const entitySet& bc2_nodes,
               const store<vec2d>& bc2_npos,
               const multiMap& face2node,
               const multiMap& node2face,
               const entitySet& bc1_edge_nodes,
               const entitySet& bc2_edge_nodes,
               const dMap& bc1p_2_bc2p) {
  int total_split_pairs = 0 ;
  // init the other global variables in special
  // routines for boundary faces
  masterbc_edge_nodes = bc1_edge_nodes ;
  slavebc_edge_nodes = bc2_edge_nodes ;
  // first we need to remove any bc1 faces that coincide
  // with one bc2 face, because we do not need to split
  // those bc1 faces
  if(global_verbose) {
    cout << "  Computing overlapped faces... " ;
    cout.flush() ;
  }

  gettimeofday(&time_face_split_frm_start,NULL) ;

  entitySet bc1_face_remove =
    bc1_removable_faces(bc1_faces,bc2_faces,face2node,
                        node2face,bc1p_2_bc2p) ;

  gettimeofday(&time_face_split_frm_end,NULL) ;
  
  if(global_verbose) {
    cout << "[" << bc1_face_remove.size()
         << " faces overlapped]" << endl ;
    cout << "  constructing quad tree... " ;
    cout.flush() ;
  }
  /******
   *  now we start to split the faces
   **/
  gettimeofday(&time_face_split_qt_start,NULL) ;
  
  rectangle bounding_box_bc2 ;
  nodes_bbox(bc2_npos,bounding_box_bc2) ;
  // make a slightly larger bounding box
  bounding_box_bc2.l -= vec2d(0.01,-0.01) ;
  bounding_box_bc2.r += vec2d(0.01,-0.01) ;

  entitySet remaining_bc1_faces = bc1_faces - bc1_face_remove ;
  // we first build a quadtree for all bc2 faces
  QuadTree bc2_qtree(15/*max_depth*/,45/*max_face_set size*/,
                    bounding_box_bc2) ;
  bc2_qtree.build_tree(bc2_faces,face2node,bc2_npos) ;
  
  gettimeofday(&time_face_split_qt_end,NULL) ;

#ifdef SHOW_QT_PROPERTY
  cout << endl << "QuadTree build time: "
       << time_face_split_qt_end - time_face_split_qt_start << endl ;
  cout << "QuadTree depth: " << bc2_qtree.depth() << endl ;
  cout << "QuadTree max face size: " << bc2_qtree.max_face_size() << endl ;
  cout << "QuadTree leaf number: " << bc2_qtree.leafs_num() << endl ;
  cout << "QuadTree total leaf size: " << bc2_qtree.all_leafs_size() << endl ;
  cout << "BC2 faces size: " << bc2_faces.size() << endl ;
#endif
  
  int finished_bc1_face = 1 ;

  float total_work_load = (float)remaining_bc1_faces.size() ;
  int total_work_points = 40 ;
  int previous_work_progress = 0 ;
  
  if(global_verbose) {
    cout << "Done" << endl ;
    cout << "  splitting "
         << get_status_bar(total_work_points,0) ;
    cout.flush() ;
  }

  gettimeofday(&time_face_split_start,NULL) ;
  for(entitySet::const_iterator ei=remaining_bc1_faces.begin();
      ei!=remaining_bc1_faces.end();++ei,++finished_bc1_face) {
    // reference to the records of bc1 face *ei
    vector<split_unit>& ei_record = global_bc1_split_records[*ei] ;
    // get all the nodes for this bc1 face
    vector<Entity> bc1_fnodes = get_face_nodes(face2node,*ei) ;
    // query the quad tree that contains bc2 faces to find
    // out the possible bc2 faces for intersection test with
    // this bc1 face
    entitySet intersect_face_cand =
      bc2_qtree.locate_polygon(bc1_fnodes,bc1_npos) ;

#ifdef FACE_SPLIT_SELF_CHECK
    vector<double> split_area_sum ;
    vector<vec2d> bc1_fpos ;
    std::set<vec2d,vec2dLess> origin_nodes ;
    std::set<vec2d,vec2dLess> split_nodes ;
    for(size_t i=0;i!=bc1_fnodes.size();++i) {
      vec2d p = bc1_npos[bc1_fnodes[i]] ;
      bc1_fpos.push_back(p) ;
      origin_nodes.insert(p) ;
    }
    double bc1_farea = polygon_area(bc1_fpos) ;
#endif

#ifdef DEBUG_SPLIT_SHOWME
    vector<vector<vec2d> > faces ;
    vector<vector<vec2d> > splits ;
    cerr << "face " << *ei << ": " << endl ;
    cerr.precision(16) ;
    
    vector<vec2d> bc1_fpos2 ;
    for(size_t i=0;i!=bc1_fnodes.size();++i) {
      vec2d p = bc1_npos[bc1_fnodes[i]] ;
      bc1_fpos2.push_back(p) ;
      cerr << bc1_fnodes[i] << ": " << p << endl ;
    }
    faces.push_back(bc1_fpos2) ;
#endif

    // test intersection for each bc2 faces found in the quad tree
    for(entitySet::const_iterator fi=intersect_face_cand.begin();
        fi!=intersect_face_cand.end();++fi) {
      // get this bc2 face nodes
      vector<int> bc2_fnodes = get_face_nodes(face2node,*fi) ;
      // compute intersection
      vector<vec2d> intersect ;
      vector<Entity> intersect_index ;

#ifdef DEBUG_SPLIT_SHOWME
      cerr << "intersecting with face " << *fi << ":" << endl ;
      vector<vec2d> bc2_fpos ;
      for(size_t i=0;i!=bc2_fnodes.size();++i) {
        vec2d p = bc2_npos[bc2_fnodes[i]] ;
        bc2_fpos.push_back(p) ;
        cerr << p << endl ;
      }
      faces.push_back(bc2_fpos) ;
#endif
      
      if(poly_and_poly(bc1_fnodes,bc1_npos,
                       bc2_fnodes,bc2_npos,
                       intersect,intersect_index)) {
        ++total_split_pairs ;
        // add split record
        // reference to the records of bc2 face *fi
        vector<split_unit>& fi_record = global_bc2_split_records[*fi] ;
        split_unit spu(intersect, intersect_index) ;
        ei_record.push_back(spu) ;
        fi_record.push_back(spu) ;

#ifdef FACE_SPLIT_SELF_CHECK
        cerr << "***Intersect***" << endl ;
        split_area_sum.push_back(polygon_area(intersect)) ;
        // check for collapsed nodes
        vector<vec2d> intersect2 = intersect ;
        std::sort(intersect2.begin(), intersect2.end(), vec2dLess()) ;
        vector<vec2d>::iterator new_end =
          std::unique(intersect2.begin(), intersect2.end(), eq_vec2d) ;
        size_t unique_nodes = new_end - intersect2.begin() ;
        if(unique_nodes != intersect.size()) {
          cerr << "Warning: face split (" << *ei << ", " << *fi
               << ") has nodes problem: " << endl ;
          cerr << "   there are " << unique_nodes
               << " unique nodes" << endl ;
          cerr << "         and " << intersect.size()
               << " total nodes" << endl ;

        }
        split_nodes.insert(intersect.begin(), intersect.end()) ;
#endif
#ifdef DEBUG_SPLIT_SHOWME
        splits.push_back(intersect) ;
#endif
      } // end of if(poly_and_poly)
    } // end of for(intersect_face_cand)

#ifdef DEBUG_SPLIT_SHOWME
    polygons2showme("split-overview.poly", faces) ;
    polygons2showme("split.poly", splits) ;
#endif

#ifdef FACE_SPLIT_SELF_CHECK
    double split_total = double_summation(split_area_sum) ;
    double area_diff = abs(bc1_farea - split_total) ;
    if(area_diff > bc1_farea/100000) {
      cerr << "Warning: split for face (" << *ei
           << ") might be problematic:" << endl ;
      cerr << "    original area:    " << bc1_farea << endl ;
      cerr << "    split total area: " << split_total << endl ;
      cerr << "    area sum diff:    " << area_diff << endl ;
    }
    vector<vec2d> set_diff ;
    std::back_insert_iterator<vector<vec2d> > diff_ii(set_diff) ;
    std::set_difference(origin_nodes.begin(), origin_nodes.end(),
                        split_nodes.begin(), split_nodes.end(),
                        diff_ii, vec2dLess()) ;
    if(!set_diff.empty()) {
      cerr << "Warning: split for face " << *ei
           << " has lost " << set_diff.size()
           << " original nodes" << endl ;
    }
#endif
    
    if(global_verbose) {
      int work_progress =
        (int)(40 * finished_bc1_face / total_work_load) ;
      if(work_progress > previous_work_progress) {
        previous_work_progress = work_progress ;
        cout << '\r' << "  splitting "
             << get_status_bar(total_work_points,work_progress) ;
        cout.flush() ;
      }
    }
  } // end of for(remaining_bc1_faces)
  gettimeofday(&time_face_split_end,NULL) ;
  if(global_verbose) {
    cout << '\r' << "  splitting... Done"
         << string(total_work_points-5,' ') << endl ;
  }
  return total_split_pairs ;
} // end of function (face_split)

// this function takes the split result and generates a new
// topology for the specified boundary mesh
// returns the set of faces processed
entitySet gen_boundary_topo
(const map<Entity, vector<split_unit> >& split_result,
 const entitySet& master_nodes,
 const store<vec2d>& master_npos,
 const entitySet& slave_nodes,
 const store<vec2d>& slave_npos,
 const multiMap& face2node,
 const multiMap& node2edge,// must contain master_node2edge
 const multiMap& edge2iface,//must contain master_edge2iface
 const Map& cl, const Map& cr,
 store<vec2d>& new_nodes,
 multiMap& new_face2node,
 multiMap& new_interior_face2node,
 dMap& new_cl, dMap& new_cr,
 dstore<alphaBeta>& inverse_info) {
  
  entitySet processed_faces = EMPTY ;
  
  dstore<vec2d> new_nodes_dstore ;
  map<Entity,vector<Entity> > new_faces ;
  map<Entity,vector<Entity> > new_interior_face2node_STL ;
  entitySet processed_nodes = master_nodes ;

  for(map<Entity,vector<split_unit> >::const_iterator
        mi=split_result.begin();mi!=split_result.end();++mi) {
    Entity master_face = mi->first ;
    processed_faces += master_face ;
    // get all the nodes for this master face
    vector<Entity> mf_nodes = get_face_nodes(face2node,master_face) ;
    // next, we triangulate this bc1 face
    vector<vector<Entity> > triangulated_mf_nodes =
      triangulate_convex_polygon(mf_nodes,master_npos) ;
    // variable to store new interior face info.
    dstore<new_interior_face_info> nifi ;
    
    const vector<split_unit>& vsu = mi->second ;
    for(vector<split_unit>::const_iterator
          vi=vsu.begin();vi!=vsu.end();++vi) {
      // first we need to remap the split polygon
      vector<Entity> remapped_index =
        remap_intersection_index(vi->ids,
                                 master_nodes,
                                 slave_nodes,
                                 mf_nodes,
                                 master_npos,
                                 slave_npos) ;
      // then we get the reverse projection info.
      for(vector<Entity>::size_type ii=0;
          ii<remapped_index.size();++ii) {
        // an index(remapped) not in the current known nodes set
        // is a new node produced in the face splitting

        // ii is the index to remapped_index vector
        // jj is the value of the ii th elements in vi->ids (Entity)
        // rr is the value of the ii th elements
        //                           in remapped_index (Entity)
        // kk is the value of the ii th elements in vi->nodes (vec2d)
        Entity jj = (vi->ids)[ii] ;
        Entity rr = remapped_index[ii] ;
        if(!processed_nodes.inSet(rr)) {
          const vec2d& kk = (vi->nodes)[ii] ;
          new_nodes_dstore[rr] = kk ;
          // compute inverse projection information
          // NOTE in computing the inverse projection info
          // we use elements in vi->ids index, NOT from the
          // remapped intersect index because all our previous
          // stored intersection info is in the original
          // nodes index as in the intersected index
          alphaBeta abinfo =
            get_inverse_proj_info(triangulated_mf_nodes,
                                  master_npos,kk,jj) ;
          inverse_info[rr] = abinfo ;
          // edit the data structure for new nodes on interior face
          check_interior_node_insertion(abinfo,node2edge,edge2iface,
                                        jj,rr,nifi) ;
          processed_nodes += rr ;
        }
      } // end of for(remapped_index)
      // set up new face information
      new_faces[mesh_new_face_index] = remapped_index ;
      new_cl[mesh_new_face_index] = cl[master_face] ;
      new_cr[mesh_new_face_index] = cr[master_face] ;
      ++mesh_new_face_index ;
    } // end of for(split_unit)
    // edit interior_face2node map
    add_new_interior_face2node(nifi,face2node,
                               new_interior_face2node_STL) ;
  } // end of for(split_result)
  // finally converts the dstore and map into store and multiMap
  dstore2store(new_nodes_dstore,new_nodes) ;
  map2multiMap(new_faces,new_face2node) ;
  map2multiMap(new_interior_face2node_STL,new_interior_face2node) ;

  // it is important to clear the global_Q_id_remap
  global_Q_id_remap.allocate(EMPTY) ;

  return processed_faces ;
}

// output the splitting to 2dgv and showme for visualization
void viz_split_showme_2dgv(const entitySet& bc1_nodes,
                           const store<vec2d> &bc1_npos,
                           const store<vec2d> &new_nodes,
                           const multiMap& new_face2node,
                           const char* poly_name,
                           const char* Tdgv_name) {
  // poly output visualize
  entitySet old_node_dom = bc1_nodes ;
  entitySet new_node_dom = new_nodes.domain() ;
  entitySet new_faces_dom = new_face2node.domain() ;
  vector<vector<vec2d> > mesh ;
  int all_nodes_number = 0 ;
  for(entitySet::const_iterator ei=new_faces_dom.begin();
      ei!=new_faces_dom.end();++ei) {
    const vector<Entity> nodes = get_face_nodes(new_face2node,*ei) ;
    vector<vec2d> polygon ;
    for(vector<Entity>::const_iterator vi=nodes.begin();
        vi!=nodes.end();++vi) {
      if(old_node_dom.inSet(*vi))
        polygon.push_back(bc1_npos[*vi]) ;
      else if(new_node_dom.inSet(*vi))
        polygon.push_back(new_nodes[*vi]) ;
      else {
        cerr << "ERROR: cannot find new node!" << endl ;
        Loci::Abort() ;
      }
    }
    all_nodes_number += nodes.size() ;
    mesh.push_back(polygon) ;
  }

  // output to poly file
  std::ofstream out(poly_name,ios::out) ;
  out.precision(16) ;
  // output nodes
  out <<  all_nodes_number << ' ' << "2 1 0" << endl ;
  int cnt = 1 ;
  for(vector<vector<vec2d> >::const_iterator vi=mesh.begin();
      vi!=mesh.end();++vi) {
    for(vector<vec2d>::const_iterator vi2=vi->begin();
        vi2!=vi->end();++vi2,++cnt)
      out << cnt << ' ' << vi2->x << ' ' << vi2->y << endl ;
    out << endl ;
  }
  // output edges
  out << all_nodes_number << " 0" << endl ;
  cnt = 1 ;
  typedef vector<vec2d>::size_type size ;

  for(vector<vector<vec2d> >::const_iterator vi=mesh.begin();
      vi!=mesh.end();++vi) {
    size p_size = vi->size() ;
    size i ;
    size begin = cnt ;
    for(i=0;i<p_size-1;++i,++cnt) {
      out << cnt << ' ' << i+begin << ' ' << i+begin+1 << endl ;
    }
    out << cnt++ << ' ' << i+begin << ' ' <<  begin << ' ' << endl ;
    out << endl ;
  }
  // 0 holes
  out << "0" << endl ;
  out.close() ;

  // output to 2dgv file
  out.open(Tdgv_name) ;
  out << "general" << endl ;
  out << all_nodes_number << ' ' << '1' << endl ;
  for(vector<vector<vec2d> >::const_iterator vi=mesh.begin();
      vi!=mesh.end();++vi) {
    for(vector<vec2d>::const_iterator vi2=vi->begin();
        vi2!=vi->end();++vi2,++cnt)
      out << vi2->x << ' ' << vi2->y << endl ;
    out << endl ;
  }
  // then output edges
  out <<  all_nodes_number << " 1 "  <<  new_faces_dom.size()
      << " 1" << endl ;
  cnt = 1 ;
  for(vector<vector<vec2d> >::const_iterator vi=mesh.begin();
      vi!=mesh.end();++vi) {
    size p_size = vi->size() ;
    size i ;
    size begin = cnt ;
    for(i=0;i<p_size-1;++i,++cnt) {
      out << i+begin << ' ' << i+begin+1 << " 1 -1" << endl ;
    }
    out << i+begin << ' ' <<  begin << ' ' << " 1 -1" << endl ;
    ++cnt ;
    out << endl ;
  }
  timeval t ;
  gettimeofday(&t,NULL) ;
  srand48(t.tv_usec) ;
  for(int i=0;i<all_nodes_number;++i)
    out << drand48() << endl ;
  out.close() ;
}

/*************************************************************
 *        Several face splitting checking functions          *
 ************************************************************/
// this function returns the number of the nodes that are only
// used by a single faces. Normally this number should not
// exceed the corner number in the original boundary mesh
void get_node2face(const multiMap&, multiMap&) ;
multiMap merge_multiMap(const multiMap&, const multiMap&) ;
int face_split_polygon_node_check(const multiMap& old_face2node,
                                  const multiMap& new_face2node,
                                  const entitySet& bc_faces,
                                  entitySet& corner_nodes) {
  int count = 0 ;
  // we first merge these two map together
  multiMap all_face2node = merge_multiMap(old_face2node,
                                          new_face2node) ;
  // we first compute an inverse map all_node2face
  multiMap all_node2face ;
  get_node2face(all_face2node,all_node2face) ;
  entitySet dom = all_node2face.domain() ;
  for(entitySet::const_iterator ei=dom.begin();ei!=dom.end();++ei) {
    entitySet faces = get_multiMap_elems(all_node2face,*ei) ;
    faces &= bc_faces ;
    int num = faces.size() ;
    if(num == 1) {
      corner_nodes += *ei ;
      ++count ;
    }
  }
  return count ;
}
// this function checks the split face (new polygons)
// definition. it detects any duplicated vertices
// definition in the polygons
// returns true to indicate test passed, false if failed (there
// are duplicated vertices in polygon definition)
bool face_split_polygon_check(const multiMap& new_face2node) {
  entitySet new_faces = new_face2node.domain() ;
  for(entitySet::const_iterator ei=new_faces.begin();
      ei!=new_faces.end();++ei) {
    vector<Entity> nodes = get_face_nodes(new_face2node,*ei) ;
    if(detect_dup(nodes)) {
      return false ;
    }
  }
  return true ;
}
// this function checks the area sum after face splitting
// it should be the same than the face area before splitting
// returns true if area is the same, false if different
bool face_split_areaSum_check(const entitySet& pre_faces,
                              const store<vec2d>& pre_pos,
                              const multiMap& pre_face2node,
                              const store<vec2d>& new_pos,
                              const multiMap& new_face2node) {
  entitySet new_faces = new_face2node.domain() ;
  entitySet new_nodes = new_pos.domain() ;
  entitySet pre_nodes = pre_pos.domain() ;
  double pre_area=0.0, new_area=0.0 ;
  vector<double> pre_area_vector, new_area_vector ;
  for(entitySet::const_iterator ei=pre_faces.begin();
      ei!=pre_faces.end();++ei) {
    vector<Entity> nodes = get_face_nodes(pre_face2node,*ei) ;
    pre_area_vector.push_back(polygon_area(nodes,pre_pos)) ;
  }
  for(entitySet::const_iterator ei=new_faces.begin();
      ei!=new_faces.end();++ei) {
    vector<Entity> nodes = get_face_nodes(new_face2node,*ei) ;
    vector<vec2d> polygon ;
    for(vector<Entity>::const_iterator vi=nodes.begin();
        vi!=nodes.end();++vi) {
      if(new_nodes.inSet(*vi))
        polygon.push_back(new_pos[*vi]) ;
      else if(pre_nodes.inSet(*vi))
        polygon.push_back(pre_pos[*vi]) ;
      else {
        cerr << "ERROR: cannot locate new nodes record!" << endl ;
        Loci::Abort() ;
      }
    }
    new_area_vector.push_back(polygon_area(polygon)) ;
  }
  pre_area = double_summation(pre_area_vector) ;
  new_area = double_summation(new_area_vector) ;

  if(abs(pre_area - new_area) <=
     pre_area / std::max(static_cast<double>(pre_faces.size()),
                         AREA_SUM_CHECK_THRESHOLD))
    return true ;
  else
    return false ;
}
// this function checks the edge map after face splitting
// returns true if passed, false if failed
bool face_split_edgeMap_check(const multiMap& new_face2node,
                              const multiMap& new_node2face) {
  entitySet new_faces = new_face2node.domain() ;
  
  entitySet processed_nodes ;
  for(entitySet::const_iterator ei=new_faces.begin();
      ei!=new_faces.end();++ei) {
    // get the number of nodes on this face
    int nf = new_face2node.end(*ei) - new_face2node.begin(*ei) ;
    // loop over all the edges on this face
    for(int i=0;i<nf;++i) {
      // get the nodes sequentially
      int n1 = new_face2node[*ei][(i==0)?(nf-1):(i-1)] ;
      int n2 = new_face2node[*ei][i] ;
      if( (processed_nodes.inSet(n1)) &&
          (processed_nodes.inSet(n2)))
        continue ;
      processed_nodes += n1 ;
      processed_nodes += n2 ;
      // get the two faces that share the edge (n1 n2)
      entitySet n1_faces = get_multiMap_elems(new_node2face,n1) ;
      entitySet n2_faces = get_multiMap_elems(new_node2face,n2) ;
      // intersection will give us the faces
      entitySet neighbor_faces = n1_faces & n2_faces ;
      neighbor_faces &= new_faces ;
      int nfsize = neighbor_faces.size() ;
      // there should only have one or two faces
      if( (nfsize == 0) || (nfsize > 2)) {
        return false ;
      }
    } // end of for(edges)
  } // end of for(new_faces)

  return true ;
}
// this function checks the convexity of faces and determines
// the smallest face area in the mesh
// returns true if all faces are convex, false otherwise,
// also sets the smallest face area and the number of faces
// whose area is less than 1e-12.
bool
face_split_convexity_area_scan
(const map<Entity,vector<split_unit> >& split_record,
 double& smallest_area, int& small_faces) {
  smallest_area = std::numeric_limits<double>::max() ;
  small_faces = 0 ;
  bool all_convex = true ;
  for(map<Entity,vector<split_unit> >::const_iterator
        mi=split_record.begin();mi!=split_record.end();++mi) {
    const vector<split_unit>& faces = mi->second ;
    for(size_t i=0;i!=faces.size();++i) {
      double a = polygon_area(faces[i].nodes) ;
      if(a < smallest_area)
        smallest_area = a ;
      if(a < 1e-12)
        ++small_faces ;
      if(all_convex)
        if(!check_convex_approx(faces[i].nodes))
          all_convex = false ;
    }
  }
  return all_convex ;
}
                               

/************************************************************
 * Several Map and Store operations for obtaining necessary *
 * new Maps and Stores for the grid geometry                *
 ***********************************************************/
// merge two multiMap to get a new one,
// any thing in m2 overides the same thing in m1
multiMap merge_multiMap(const multiMap& m1, const multiMap& m2) {
  multiMap ret ;
  store<int> alloc ;
  entitySet m1dom = m1.domain() ;
  entitySet m2dom = m2.domain() ;
  alloc.allocate(m1dom | m2dom) ;
  for(entitySet::const_iterator ei=m1dom.begin();
      ei!=m1dom.end();++ei) {
    if(!m2dom.inSet(*ei))
      alloc[*ei] = m1.num_elems(*ei) ;
  }
  for(entitySet::const_iterator ei=m2dom.begin();
      ei!=m2dom.end();++ei)
    alloc[*ei] = m2.num_elems(*ei) ;

  ret.allocate(alloc) ;
  for(entitySet::const_iterator ei=m1dom.begin();
      ei!=m1dom.end();++ei) {
    if(!m2dom.inSet(*ei)) {
      int num = m1.num_elems(*ei) ;
      for(int i=0;i<num;++i)
        ret[*ei][i] = m1[*ei][i] ;
    }
  }
  for(entitySet::const_iterator ei=m2dom.begin();
      ei!=m2dom.end();++ei) {
    int num = m2.num_elems(*ei) ;
    for(int i=0;i<num;++i)
      ret[*ei][i] = m2[*ei][i] ;
  }
  return ret ;
}
// get the boundary edge information
void get_edge2interiorFace(const dMap& N1, const dMap& N2,
                           const dMap& El, const dMap& Er,
                           const multiMap& node2face,
                           multiMap& edge2iface) {
  map<Entity,entitySet> STL_edge2iface ;
  entitySet dom = N1.domain() ;
  for(entitySet::const_iterator ei=dom.begin();
      ei!=dom.end();++ei) {
    entitySet local_face ; local_face += El[*ei] ;
    if(Er[*ei] >= 0)
      local_face += Er[*ei] ;
    entitySet n1_faces = get_multiMap_elems(node2face,N1[*ei]) ;
    entitySet n2_faces = get_multiMap_elems(node2face,N2[*ei]) ;
    entitySet interior_faces = n1_faces & n2_faces ;
    interior_faces -= local_face ;
    // now we have found all the interior faces
    // and we have to build the map
    /*
    for(entitySet::const_iterator fi=interior_faces.begin();
        fi!=interior_faces.end();++fi) {
      STL_edge2iface[*ei] += *fi ;
    }
    */
    STL_edge2iface[*ei] += interior_faces ;
  }
  // converts to multiMap
  map2multiMap(STL_edge2iface,edge2iface) ;
}

// this function inverses the N1, N2 relation to get
// the new node2edge map
void get_node2edge(const dMap& N1, const dMap& N2,
                   multiMap& node2edge) {
  // we first contruct a STL map
  map<Entity,entitySet> STL_node2edge ;
  entitySet dom = N1.domain() ;
  for(entitySet::const_iterator ei=dom.begin();
      ei!=dom.end();++ei)
    STL_node2edge[N1[*ei]] += *ei ;
  dom = N2.domain() ;
  for(entitySet::const_iterator ei=dom.begin();
      ei!=dom.end();++ei)
    STL_node2edge[N2[*ei]] += *ei ;

  map2multiMap(STL_node2edge,node2edge) ;
}

// this function computes the node2node relation.
// i.e., given a node, through the node2node map,
// we can find all other nodes connected by edges
void get_node2node(const entitySet& nodes,
                   const dMap& N1, const dMap& N2,
                   const multiMap& node2edge,
                   multiMap& node2node) {
  // we first construct an STL map
  map<Entity,entitySet> STL_node2node ;
  for(entitySet::const_iterator ei=nodes.begin();
      ei!=nodes.end();++ei) {
    // find all the edges
    entitySet others ;
    entitySet edges = get_multiMap_elems(node2edge,*ei) ;
    for(entitySet::const_iterator edi=edges.begin();
        edi!=edges.end();++edi) {
      Entity n1 = N1[*edi] ;
      Entity n2 = N2[*edi] ;
      if(n1 != *ei)
        others += n1 ;
      else
        others += n2 ;
    }
    STL_node2node[*ei] = others ;
  }

  map2multiMap(STL_node2node,node2node) ;
}
// this function inverses the face2node relation to get
// the new node2face map
void get_node2face(const multiMap& face2node,
                   multiMap& node2face) {
  // we first contruct a STL map
  map<Entity,entitySet> STL_node2face ;
  entitySet dom = face2node.domain() ;
  for(entitySet::const_iterator ei=dom.begin();
      ei!=dom.end();++ei) {
    vector<Entity> nodes = get_face_nodes(face2node,*ei) ;
    for(vector<Entity>::const_iterator vi=nodes.begin();
        vi!=nodes.end();++vi)
      STL_node2face[*vi] += *ei ;
  }
  map2multiMap(STL_node2face,node2face) ;
}

/************************************************************
 *  3D points projection: point to point match, inverse     *
 *  project one boundary, and transform to another boundary *
 *  The functions in this section are used to get the       *
 *  3D positions for all nodes in two split boundaries      *
 ************************************************************/
// inverse projection of 2D points back to 3D space
// according to the inverse projection information
void inverse_project(const dstore<alphaBeta>& info,
                     const store<vec3d>& pos3d,
                     store<vec3d>& new_pos) {
  entitySet dom = info.domain() ;
  new_pos.allocate(EMPTY) ;
  new_pos.allocate(dom) ;
  for(entitySet::const_iterator ei=dom.begin();
      ei!=dom.end();++ei) {
    const alphaBeta& i = info[*ei] ;
    const vec3d& P0 = pos3d[i.p0] ;
    const vec3d& P1 = pos3d[i.p1] ;
    const vec3d& P2 = pos3d[i.p2] ;
    new_pos[*ei] = P0 + i.a * (P1 - P0) + i.b * (P2 - P0) ;
  }
}
/**
***************************************
* The 3D space transformation objects *
**/
inline vec3d dot(ten3d t, const vec3d v) {
  return vec3d(dot(t.x,v),dot(t.y,v),dot(t.z,v)) ;
}

struct rigid_transform {
  vec3d t1,t2 ;
  ten3d R,Rinv ;
  rigid_transform(vec3d center, vec3d v, double angle, vec3d translate) {
    t1 = -1.*center ;
    t2 = center + translate ;
    R.x = vec3d(1,0,0) ;
    R.y = vec3d(0,1,0) ;
    R.z = vec3d(0,0,1) ;
    Rinv = R ;
    if(angle == 0)
      return ;
    double s = sin(angle) ;
    double c = cos(angle) ;
    double C = 1.-c ;
    double x = v.x ;
    double y = v.y ;
    double z = v.z ;
    R.x += vec3d(C*(x*x-1.),C*x*y-z*s,C*x*z+y*s) ;
    R.y += vec3d(C*x*y+z*s,C*(y*y-1.),C*y*z-x*s) ;
    R.z += vec3d(C*x*z-y*s,C*y*z+x*s,C*(z*z-1.)) ;
    s = -s ;
    Rinv.x += vec3d(C*(x*x-1.),C*x*y-z*s,C*x*z+y*s) ;
    Rinv.y += vec3d(C*x*y+z*s,C*(y*y-1.),C*y*z-x*s) ;
    Rinv.z += vec3d(C*x*z-y*s,C*y*z+x*s,C*(z*z-1.)) ;
  }
  vec3d rotate_vec(vec3d v) const {
    return dot(R,v) ;
  }
  vec3d transform(vec3d v) const {
    return rotate_vec(v+t1)+t2 ;
  }
} ;
// This function is not used in this section. It is left
// from previous code development. However since this is a
// general function that computes facecenter for arbitrary
// faces, it might be useful. We therefore keep it here.
vec2d get_facecenter(const vector<vec2d>& nodes) {
  vector<vec2d>::size_type size = nodes.size() ;
  if(size < 3) {
    cerr << "ERROR: face ill defined in get_facecenter!" << endl ;
    Loci::Abort() ;
  }
  vector<vec2d>::size_type i, j ;
  double atmp = 0., xtmp = 0., ytmp = 0. ;
  for(i=size-1,j=0;j<size;i=j,++j) {
    double ai = nodes[i].x * nodes[j].y - nodes[j].x * nodes[i].y ;
    atmp += ai ;
    xtmp += (nodes[j].x + nodes[i].x) * ai ;
    ytmp += (nodes[j].y + nodes[i].y) * ai ;
  }
  // atmp is the polygon's area * 2
  if(atmp == 0) {
    cerr << "ERROR: face degenerated in get_facecenter!" << endl ;
    Loci::Abort() ;
  }
  return vec2d(xtmp/(3*atmp), ytmp/(3*atmp)) ;
}
// this function computes the minimal edge length of the passed in face
double min_face_edge_len(const vector<vec2d>& nodes) {
  double min_len = std::numeric_limits<double>::max() ;
  vector<vec2d>::size_type size = nodes.size() ;
  if(size < 2)
    return 0.0 ;
  vector<vec2d>::size_type i, j ;
  for(i=size-1,j=0;j<size;i=j,++j) {
    double len = dist(nodes[i], nodes[j]) ;
    if(len < min_len)
      min_len = len ;
  }
  return min_len ;
}
// this function matches corresponding points on two boundaries
void twoD_point_match(const store<vec2d>& bc1_nodes,
                      const store<vec2d>& bc2_nodes,
                      const store<vec2d>& bc1_new_nodes,
                      const store<vec2d>& bc2_new_nodes,
                      const entitySet& bc1_all_faces,
                      const multiMap& face2node,
                      const multiMap& bc1_new_face2node,
                      dMap& bc2p_2_bc1p,
                      entitySet& bc1_not_matched,
                      entitySet& bc2_not_matched) {
  // some entitySet useful
  entitySet bc1_orig_nodes_dom = bc1_nodes.domain() ;
  entitySet bc2_orig_nodes_dom = bc2_nodes.domain() ;
  entitySet bc1_new_nodes_dom = bc1_new_nodes.domain() ;
  entitySet bc2_new_nodes_dom = bc2_new_nodes.domain() ;
  entitySet orig_faces_dom = face2node.domain() ;
  entitySet bc1_new_faces_dom = bc1_new_face2node.domain() ;
  // this records the minimal face edge length on bc1
  double min_len = std::numeric_limits<double>::max() ;
  for(entitySet::const_iterator ei=bc1_all_faces.begin();
      ei!=bc1_all_faces.end();++ei) {
    vector<int> nodes_id ;
    vector<vec2d> nodes ;
    // find out all the nodes for this faces
    if(orig_faces_dom.inSet(*ei)) {
      nodes_id = get_face_nodes(face2node,*ei) ;
    }else if(bc1_new_faces_dom.inSet(*ei)) {
      nodes_id = get_face_nodes(bc1_new_face2node,*ei) ;
    }else {
      cerr << "ERROR: face definition error in twoD_point_match!"
           << endl ;
      Loci::Abort() ;
    }
    for(vector<int>::const_iterator vi=nodes_id.begin();
        vi!=nodes_id.end();++vi) {
      if(bc1_orig_nodes_dom.inSet(*vi))
        nodes.push_back(bc1_nodes[*vi]) ;
      else if(bc1_new_nodes_dom.inSet(*vi))
        nodes.push_back(bc1_new_nodes[*vi]) ;
      else {
        cerr << "ERROR: read node position error in twoD_point_match"
             << endl ;
        Loci::Abort() ;
      }
    }
    // compute the minimal edge length for this face
    double len = min_face_edge_len(nodes) ;
    if(len < min_len)
      min_len = len ;
  }
  // union all node records
  store<vec2d> bc1n, bc2n ;
  entitySet bc1ndom = bc1_orig_nodes_dom | bc1_new_nodes_dom ;
  entitySet bc2ndom = bc2_orig_nodes_dom | bc2_new_nodes_dom ;
  // checking...
  if(bc1ndom.size() != bc2ndom.size()) {
    cout << endl ;
    cout << "Error: nodes number not equal on two boundaries" << endl ;
    cout << "       please report this problem!" << endl ;
    Loci::Abort() ;
  }
  bc1n.allocate(bc1ndom) ;
  bc2n.allocate(bc2ndom) ;
  for(entitySet::const_iterator ei=bc1_orig_nodes_dom.begin();
      ei!=bc1_orig_nodes_dom.end();++ei)
    bc1n[*ei] = bc1_nodes[*ei] ;
  for(entitySet::const_iterator ei=bc1_new_nodes_dom.begin();
      ei!=bc1_new_nodes_dom.end();++ei)
    bc1n[*ei] = bc1_new_nodes[*ei] ;
  for(entitySet::const_iterator ei=bc2_orig_nodes_dom.begin();
      ei!=bc2_orig_nodes_dom.end();++ei)
    bc2n[*ei] = bc2_nodes[*ei] ;
  for(entitySet::const_iterator ei=bc2_new_nodes_dom.begin();
      ei!=bc2_new_nodes_dom.end();++ei)
    bc2n[*ei] = bc2_new_nodes[*ei] ;
  // now we need to construct a PointSet QuadTree for bc1 nodes
  rectangle bounding_box ;
  nodes_bbox(bc1n, bounding_box) ;
  // make it slightly larger
  bounding_box.l -= vec2d(0.01, -0.01) ;
  bounding_box.r += vec2d(0.01, -0.01) ;
  PointSet_QuadTree bc1_qtree(bounding_box) ;
  bc1_qtree.build_tree(bc1ndom, bc1n) ;
  // matching...
  for(entitySet::const_iterator ei=bc2ndom.begin();
      ei!=bc2ndom.end();++ei) {
    const vec2d& node = bc2n[*ei] ;
    // we search through the quadtree to match this node
    entitySet cand = bc1_qtree.locate_circle(node,min_len) ;
    if(cand == EMPTY) {
      // then this node does not match
      continue ;
    }
    // then we compute if this center matches
    int partner = -1 ;
    double min_dist = std::numeric_limits<double>::max() ;
    for(entitySet::const_iterator ei2=cand.begin();
        ei2!=cand.end();++ei2) {
      double d = dist(node,bc1n[*ei2]) ;
      if(d < min_dist) {
        partner = *ei2 ;
        min_dist = d ;
      }
    }
    // fill in the match records
    bc2p_2_bc1p[*ei] = partner ;
  }
  // then we will check the matching results to make
  // sure nothing wrong happened
  // first all bc2 nodes must have a record in bc2p_2_bc1p
  entitySet matching_dom = bc2p_2_bc1p.domain() ;
  bc2_not_matched = bc2ndom - matching_dom ;
  bc1_not_matched = bc1ndom - bc2p_2_bc1p.image(matching_dom) ;
}
// this function projects back the 2d boundaries into their 3d space.
// bc1 is inverse projected by interpolation (using the info obtained
// from face splitting), bc2 is projected by transforming bc1 points
// directly. This is to minimize the possible floating point errors.
void get_3D_nodes_pos(const dstore<alphaBeta>& bc1_inverse_proj_info,
                      const dMap& bc2p_2_bc1p,
                      const rigid_transform& rt,
                      const entitySet& bc1_orig_nodes,
                      const entitySet& bc1_new_nodes,
                      const entitySet& bc2_orig_nodes,
                      const entitySet& bc2_new_nodes,
                      store<vec3d>& orig_3d_pos,
                      store<vec3d>& bc1_new_3d_pos,
                      store<vec3d>& bc2_new_3d_pos) {
  // first we inverse project bc1 nodes to their 3d pos
  inverse_project(bc1_inverse_proj_info,orig_3d_pos,bc1_new_3d_pos) ;
  // then we need to obtain the bc2 nodes' 3d position
  bc2_new_3d_pos.allocate(bc2_new_nodes) ;
  entitySet transform_domain = bc2p_2_bc1p.domain() ;
  for(entitySet::const_iterator ei=transform_domain.begin();
      ei!=transform_domain.end();++ei) {
    if(bc2_orig_nodes.inSet(*ei)) {
      int partner = bc2p_2_bc1p[*ei] ;
      if(bc1_orig_nodes.inSet(partner))
        orig_3d_pos[*ei] = rt.transform(orig_3d_pos[partner]) ;
      else if(bc1_new_nodes.inSet(partner))
        orig_3d_pos[*ei] = rt.transform(bc1_new_3d_pos[partner]) ;
      else {
        cout << "Error in obtaining 3D position for boundary nodes" << endl ;
        cout << "Please report this problem." << endl ;
      }
    }else if(bc2_new_nodes.inSet(*ei)) {
      int partner = bc2p_2_bc1p[*ei] ;
      if(bc1_orig_nodes.inSet(partner))
        bc2_new_3d_pos[*ei] = rt.transform(orig_3d_pos[partner]) ;
      else if(bc1_new_nodes.inSet(partner))
        bc2_new_3d_pos[*ei] = rt.transform(bc1_new_3d_pos[partner]) ;
      else {
        cout << "Error in obtaining 3D position for boundary nodes" << endl ;
        cout << "Please report this problem." << endl ;
      }
    }else {
      // an error
      cout << "Error in obtaining 3D position for boundary nodes!" << endl ;
      cout << "Please report this problem." << endl ;
      Loci::Abort() ;
    }
  }
}
// visualize the face set for debugging purpose
void facecenter_result_viz(const entitySet& faces,
                           const store<vec2d>& nodes1,
                           const store<vec2d>& nodes2,
                           const multiMap& face2node1,
                           const multiMap& face2node2,
                           const char* output_name) {
  // we will use the write_out_face_poly function
  // we first construct a vector<vector<vec2d> > for the faces
  entitySet n1dom = nodes1.domain() ;
  entitySet n2dom = nodes2.domain() ;
  entitySet f2n1dom = face2node1.domain() ;
  entitySet f2n2dom = face2node2.domain() ;
  
  vector<vector<vec2d> > mesh ;
  for(entitySet::const_iterator ei=faces.begin();
      ei!=faces.end();++ei) {
    vector<vec2d> nodes ;
    vector<int> nodes_id ;
    if(f2n1dom.inSet(*ei))
      nodes_id = get_face_nodes(face2node1,*ei) ;
    else
      nodes_id = get_face_nodes(face2node2,*ei) ;
    for(vector<int>::const_iterator vi=nodes_id.begin();
        vi!=nodes_id.end();++vi)
      if(n1dom.inSet(*vi))
        nodes.push_back(nodes1[*vi]) ;
      else
        nodes.push_back(nodes2[*vi]) ;
    mesh.push_back(nodes) ;
  }
  write_out_face_poly(mesh,output_name) ;
}
// debugging function that computes the face center
// for the passed in faces and returns a vector of face center
vector<vec2d> comp_face_set_center(const entitySet& faces,
                                   const store<vec2d>& nodes1,
                                   const store<vec2d>& nodes2,
                                   const multiMap& face2node1,
                                   const multiMap& face2node2) {
  entitySet n1dom = nodes1.domain() ;
  entitySet n2dom = nodes2.domain() ;
  entitySet f2n1dom = face2node1.domain() ;
  entitySet f2n2dom = face2node2.domain() ;
  
  vector<vec2d> centers ;
  for(entitySet::const_iterator ei=faces.begin();
      ei!=faces.end();++ei) {
    vector<vec2d> nodes ;
    vector<int> nodes_id ;
    if(f2n1dom.inSet(*ei))
      nodes_id = get_face_nodes(face2node1,*ei) ;
    else
      nodes_id = get_face_nodes(face2node2,*ei) ;
    for(vector<int>::const_iterator vi=nodes_id.begin();
        vi!=nodes_id.end();++vi)
      if(n1dom.inSet(*vi))
        nodes.push_back(nodes1[*vi]) ;
      else
        nodes.push_back(nodes2[*vi]) ;
    centers.push_back(get_facecenter(nodes)) ;
  }
  return centers ;
}
//////////// End of 3D point projection /////////////

// a small utility function that extracts the 3D vector definition
// from the command line string.
// return true to indicate success, false if any failure.
// the extracted string is stored in s1, s2, and s3
// the vector is store as 0,0,0 format
inline bool extract_vec3d(const string& s,
                          string& s1, string& s2, string& s3) {
  // s should at least have a length of 7 (e.g., 0,0,0)
  if(s.size() < 5)
    return false ;
  string::const_iterator b = s.begin() ;
  string::const_iterator e = s.end() ;
  string::const_iterator si,si2 ;
  // find the first comma
  si = find(b,e,',') ;
  if(si == e) // if not found, return false
    return false ;
  s1 = string(b,si) ;
  // find the second blank
  si2 = find(si+1,e,',') ;
  if(si2 == e)
    return false ;
  s2 = string(si+1,si2) ;
  s3 = string(si2+1,e) ;
  if(s1.empty() || s2.empty() || s3.empty())
    return false ;
  
  return true ;
}
std::ostream& show_brief_msg(const string& pb, std::ostream& s) {
  s << pb << " is a periodic boundary tool for the CHEM program." << endl ;
  s << "Version 0.123 (Try -h option for detailed help)" << endl ;

  return s ;
}
// show the help message
std::ostream& show_usage(const string& pb, std::ostream& s) {
  s << pb << " is a periodic boundary tool. It is mainly designed " << endl ;
  s << "to be used with the CHEM program (the vog grid format). " << endl ;
  s << "It makes exactly same boundaries for the given grid. " << endl ;
  s << "This version uses the robust geometric predicates from" << endl ;
  s << "Jonathan Richard Shewchuk." << endl ;
  s << endl << "This is version 0.123" << endl ;
  s << endl ;
  s << "Some preconditions must be met before using this tool" << endl ;
  s << "  1) The two specified boundaries must have same contours" << endl ;
  s << "  2) The boundary cannot have arbitrary complex topology," << endl ;
  s << "         planes or simple curves are the best inputs" << endl ;
  s << "  3) Boundary faces must be convex" << endl ;
  s << endl ;
  s << "Usage: " << pb << " [options] input_case" << endl ;
  s << "    These options are available:" << endl ;
  s << "    -h  Help: show this help message you are reading" << endl ;
  s << "    -o  Specify the output grid file name (default is pb_casename.vog)"
    << endl ;
  s << "    -Q  Quiet: No terminal output except errors" << endl ;
  s << "    -2  Output intermediate results for 2D visualization" << endl ;
  s << "          program \"2dgv\" and \"showme\"" << endl ;
  s << "    -b  Specify two boundaries. Followed by two boundary names" << endl ;
  s << "          separated by a space. E.g., -b periodic1 periodic2" << endl ;
  s << "    -r  Specify rotation vectors. \"center\" \"axis\" and \"angle\""
    << endl ;
  s << "          respectively. The rotation is specified from the first"
    << endl ;
  s << "          boundary to the second boundary. Vectors are in 3D space,"
    << endl ;
  s << "          and in the format of 0,0,0 angle is in degree" << endl ;
  s << "    -t  Specify the translation vector. The translation is " << endl ;
  s << "          specified from the first boundary to the second boundary."
    << endl ;
  s << "          The vector is in the format of 0,0,0" << endl ;
  s << endl ;
  s << "    Followings are some advanced options:" << endl ;
  s << "    --shiftrange  The threshold in the point shifting phase."
    << endl ;
  s << "                    During the point shifting phase, one boundary's"
    << endl ;
  s << "                    nodes are shifted to match the corresponding nodes"
    << endl ;
  s << "                    on the other boundary. The threshold for such"
    << endl ;
  s << "                    shifting is the minimal distance of the point"
    << endl ;
  s << "                    to its bounding polygon times \"shiftrange.\""
    << endl ;
  s << "                    This value is valid in range [0,1] only." << endl ;
  s << "                    To reset it, supply a valid number, e.g., "
    << endl ;
  s << "                    --shiftrange 0.05. The default is 0.05" << endl ;
  s << "    --cta         The collinear threshold angle. This value is"
    << endl ;
  s << "                    the tolerance for the inexactness occurred at"
    << endl ;
  s << "                    the mesh boundary edges. Specifically this value"
    << endl ;
  s << "                    is used to control the collinearity predicate for"
    << endl ;
  s << "                    boundary edges and nodes. If the largest internal"
    << endl ;
  s << "                    angle in the triangle formed by three nodes is"
    << endl ;
  s << "                    within this threshold to 180 degrees, they are"
    << endl ;
  s << "                    still considered collinear. To reset this value,"
    << endl ;
  s << "                    supply a valid number in degree, e.g., --cta 1e-3."
    << endl ;
  s << "                    The default is approximately 1e-3 degree, which"
    << endl ;
  s << "                    means the largest internal angle in the triangle"
    << endl ;
  s << "                    by three points is allowed to be within" << endl ;
  s << "                    [180-1e-3,180+1e-3] for the three points to be"
    << endl ;
  s << "                    considered collinear" << endl ;
  s << "    --ast         The tolerance used in performing area sum check."
    << endl ;
  s << "                    The default value is set to be 1e8. Specifically,"
    << endl ;
  s << "                    if the sum of split areas are off by " << endl ;
  s << "                    bounary area / threshold, then it is considered"
    << endl ;
  s << "                    a failure. Setting a smaller value will decrease"
    << endl ;
  s << "                    the sensitivity. To reset this value, supply"
    << endl ;
  s << "                    a valid double number, e.g., --cta 1e7" << endl ;
  s << "    --nocheck     Do not perform any of the after split checks"
    << endl ;
  s << "    NOTE: These options change the program's internal behavior and"
    << endl ;
  s << "            are NOT recommended to use unless you fully understand the"
    << endl ;
  s << "            consequence of changing the default value." << endl ;
  s << endl ;
  s << "Final note: the order of options and input_case does not matter"
    << endl ;
  s << "      input_case is the grid file name without \".vog\" suffix"
    << endl ;
  s << "      (e.g., if you have \"abc.vog\", then input_case is \"abc\")"
    << endl ;
  s << "      the options \"-b\" and either \"-r\" or \"-t\" must be"
    << endl ;
  s << "      specified. But \"-r\" and \"-t\" cannot be both specified."
    << endl ;
  s << "      this function is not currently supported." << endl ;
  s << endl ;
  s << "Examples: " << endl ;
  s << "     " << pb << " -o abc-pb -r 0,0,0 1,0,0 30 "
    << " -b periodic1 periodic2 abc" << endl ;
  s << "     " << pb << " abc -o abc-pb -b periodic1 periodic2 -t 1,0,0" << endl ;
  
  return s ;
}

namespace Loci {
  bool readGridVOG(vector<entitySet> &local_nodes,
		   vector<entitySet> &local_faces,
		   vector<entitySet> &local_cells,
		   store<vector3d<double> > &pos, Map &cl, Map &cr,
		   multiMap &face2node, 
		   store<string> &boundary_names,
		   store<string> &boundary_tags,
                   vector<pair<string,entitySet> > &volTags,
		   int max_alloc, string filename) ;
}
///////////////////////////////////////////////////////////
//                      The Main                         //
///////////////////////////////////////////////////////////
int main(int ac, char* av[]) {
  gettimeofday(&time_prog_start,NULL) ;
  
  Loci::Init(&ac, &av) ;
  // init for exact arithmetic
  exactinit() ;
  bool optimize = true;
  /**********
   * The command line parser
   **/
  // the input problem name (e.g., "a" for a.vog)
  string problem_name ;
  // flag to indicate whether problem name has been specified or not
  bool problem_name_seen = false ;
  // the output grid file name
  string output_name = "" ;
  // flag to indicate whether output name has been specified or not
  bool output_name_seen = false ;
  // flag to indicate is the boundaries rotates?
  bool is_rotate = false ;
  // or they are translated?
  bool is_translate = false ;
  // rotation information
  vec3d rotation_center(0,0,0), rotation_axis(0,0,0) ;
  double rotation_angle_in_degree = 0 ;
  // translation information
  vec3d translation_vector(0,0,0) ;
  // flags to check multiple specification of rotation and translation
  bool rotation_set = false ;
  bool translation_set = false ;
  // indicate the boundary id, they actually should be
  // "Entity" rather than "entitySet". Here we just set
  // it for later convenience

  int BC1_id=0, BC2_id=0 ;
  string BC1name, BC2name ;
  // flag to indicate whether boundaries are set or not
  bool boundary_set = false ;
  // whether to suppress output or not
  // (but error messages are NOT suppressed)
  bool verbose = true ;
  // whether to output intermediate results for
  // 2D visualization program (2dgv & showme)
  bool twoD_aux_viz = false ;
  // flag to indicate whether we have set shift_range or not
  bool shift_range_set = false ;
  // const factor multiplied in GRID_POINT_COMPARE_THRESHOLD
  // defaults to 1e-8
  double collinear_threshold_angle ;
  bool collinear_threshold_angle_set = false ;
  bool area_sum_check_threshold_set = false ;
  bool report_timing = false ;
  // whether to print help message or not
  bool print_help = false ;
  // brief description message
  bool brief_help = false ;
  // whether to perform after split check
  bool after_split_check = true ;

  string pb_name = av[0] ;
  int avi = 1 ;

  if(ac==1)
    brief_help = true ;
  
  while(avi < ac) {
    string opt = av[avi] ;
    if(opt[0] == '-') {
      if(opt == "-o") {
        if(output_name_seen) {
          cout << pb_name << ": warning: option -o appears more than once, "
               << "first setting taken" << endl ;
          avi+=2 ;
          continue ;
        }
        output_name_seen = true ;
        output_name = av[avi+1] ;
        avi+=2 ;
      } else if(opt == "-Q") {
        verbose = false ;
        global_verbose = false ;
        avi+=1 ;
      } else if(opt == "--time") {
        report_timing = true ;
        avi+=1 ;
      } else if(opt == "-2") {
        twoD_aux_viz = true ;
        avi+=1 ;
      } else if(opt == "-h") {
        print_help = true ;
        avi+=1 ;
      } else if(opt == "-b") {
        if(boundary_set) {
          cout << pb_name << ": warning: option -b appears more than once, "
               << "first setting taken" << endl ;
          avi+=3 ;
          continue ;
        }
        BC1name = string(av[avi+1]) ;
        BC2name = string(av[avi+2]) ;
        
        boundary_set = true ;
        avi+=3 ;
      } else if(opt == "-r") {
        if(rotation_set) {
          cout << pb_name << ": warning: option -r appears more than once, "
               << "first setting taken" << endl ;
          avi+=4 ;
          continue ;
        }
        // then get the "center" "axis" and "angle" respectively
        // "center" and "axis" are 3d vectors
        string center_x, center_y, center_z ;
        string axis_x, axis_y, axis_z ;
        string angle_string ;
        if(!extract_vec3d(string(av[avi+1]),center_x,center_y,center_z)) {
          cout << pb_name << ": error: \"center\" vector should be "
               << "specified as 0,0,0 format" << endl ;
          return -1 ;
        }
        if(!extract_vec3d(string(av[avi+2]),axis_x,axis_y,axis_z)) {
          cout << pb_name << ": error: \"axis\" vector should be "
               << "specified as 0,0,0 format" << endl ;
          return -1 ;
        }
        angle_string = string(av[avi+3]) ;
        if(!valid_double(center_x) ||
           !valid_double(center_y) ||
           !valid_double(center_z)) {
          cout << pb_name << ": error: \"center\" vector incorrect,"
               << " vector elements are NOT valid double numbers" << endl ;
          return -1 ;
        }
        if(!valid_double(axis_x) ||
           !valid_double(axis_y) ||
           !valid_double(axis_z)) {
          cout << pb_name << ": error: \"axis\" vector incorrect,"
               << " vector elements are NOT valid double numbers" << endl ;
          return -1 ;
        }
        if(!valid_double(angle_string)) {
          cout << pb_name << ": error: rotation angle incorrect,"
               << " it is NOT a valid double number" << endl ;
          return -1 ;
        }
        rotation_center = vec3d(str2double(center_x),
                                str2double(center_y),
                                str2double(center_z)) ;
        rotation_axis = vec3d(str2double(axis_x),
                              str2double(axis_y),
                              str2double(axis_z)) ;
        rotation_angle_in_degree = str2double(angle_string) ;

        is_rotate = true ;
        rotation_set = true ;
        avi+=4 ;
      } else if(opt == "-t") {
        if(translation_set) {
          cout << pb_name << ": warning: option -t appears more than once, "
               << "first setting taken" << endl ;
          avi+=2 ;
          continue ;
        }
        // then get the "translate" vector it is a 3d vectors
        string tran_x, tran_y, tran_z ;
        if(!extract_vec3d(string(av[avi+1]),tran_x,tran_y,tran_z)) {
          cout << pb_name << ": error: \"translate\" vector should be "
               << "specified as 0,0,0 format" << endl ;
          return -1 ;
        }
        if(!valid_double(tran_x) ||
           !valid_double(tran_y) ||
           !valid_double(tran_z)) {
          cout << pb_name << ": error: \"translate\" vector incorrect,"
               << " vector elements are NOT valid double numbers" << endl ;
          return -1 ;
        }
        translation_vector = vec3d(str2double(tran_x),
                                   str2double(tran_y),
                                   str2double(tran_z)) ;
        is_translate = true ;
        translation_set = true ;
        avi+=2 ;
      } else if(opt == "--shiftrange") {
        if(shift_range_set) {
          cout << pb_name << ": warning: option --shiftrange appears "
               << "more than once, first setting taken" << endl ;
          avi+=2 ;
          continue ;
        }
        string range_str = av[avi+1] ;
        if(!valid_double(range_str)) {
          cout << pb_name << ": error: shift range is not a valid "
               << "double number" << endl ;
          return -1 ;
        }
        double range = str2double(range_str) ;
        if( (range<0) || (range>1)) {
          cout << pb_name << ": error: shift range should be "
               << "within [0,1]" << endl ;
          return -1 ;
        }
        
        SHIFT_RANGE_THRESHOLD = range ;
        shift_range_set = true ;
        avi+=2 ;
      } else if(opt == "--cta") {
        if(collinear_threshold_angle_set) {
          cout << pb_name << ": warning: option --cta appears "
               << "more than once, first setting taken" << endl ;
          avi+=2 ;
          continue ;
        }
        string angle_str = av[avi+1] ;
        if(!valid_double(angle_str)) {
          cout << pb_name << ": error: collinear threshold angle is "
               << "not a valid double number" << endl ;
          return -1 ;
        }
        collinear_threshold_angle = str2double(angle_str) ;
        ON_LINE_SEG_DOT_THRESHOLD =
          1.0 - cos(deg2rad(collinear_threshold_angle)) ;
        collinear_threshold_angle_set = true ;
        avi+=2 ;
      } else if(opt == "--ast") {
        if(area_sum_check_threshold_set) {
          cout << pb_name << ": warning: option --ast appears "
               << "more than once, first setting taken" << endl ;
          avi+=2 ;
          continue ;
        }
        string ast_str = av[avi+1] ;
        if(!valid_double(ast_str)) {
          cout << pb_name << ": error: area sum check threshold is "
               << "not a valid double number" << endl ;
          return -1 ;
        }
        AREA_SUM_CHECK_THRESHOLD = str2double(ast_str) ;
        area_sum_check_threshold_set = true ;
        avi+=2 ;
      } else if (opt == "--nocheck") {
        after_split_check = false ;
        avi+=1 ;
      } else {
        cout << pb_name << ": warning: option " << opt
             << " not understood! ignored"
             << endl ;
        avi+=1 ;
      }
    } else {
      if(problem_name_seen) {
        cout << pb_name << ": error: multiple problem name specified! "
             << "Confused!" << endl ;
        return -1 ;
      }
      problem_name = opt ;
      problem_name_seen = true ;
      avi+=1 ;
    }
  }

  // post checking for command line parse
  if(brief_help) {
    show_brief_msg(pb_name,cout) ;
    return 0 ;
  }
  if(print_help) {
    show_usage(pb_name,cout) ;
    return 0 ;
  }
  if(is_rotate && is_translate) {
    cout << pb_name << ": error: both rotation and translation are NOT "
         << "supported currently" << endl ;
    return -1 ;
  }
  if(!is_rotate && !is_translate) {
    cout << pb_name << ": error: either rotation or translation must be "
         << "specified" << endl ;
    return -1 ;
  }
  if(!problem_name_seen) {
    cout << pb_name << ": error: no input problem specified!" << endl ;
    return -1 ;
  }
  if(!boundary_set) {
    cout << pb_name << ": error: no boundaries specified!" << endl ;
    return -1 ;
  }
  /***************************
   ** command line parser ends
   **/
#ifdef NORMAL
  /*******************************************************
   *            Actual operations begin here             *
   *******************************************************/

  //read in grid
  fact_db facts ;
  string gridfile = problem_name + string(".vog") ;
  vector<pair<int,string> > boundary_ids ;
  if(!Loci::readBCfromVOG(gridfile,boundary_ids)) {
    cerr << "unable to open grid file '" << gridfile << "'" << endl ;
    exit(-1) ;
  }
  BC1_id = -1 ;
  BC2_id = -1 ;
  if(boundary_ids.size()==0){
    BC1_id=atoi(BC1name.substr(3).c_str());
    BC2_id=atoi(BC2name.substr(3).c_str()); 
  }else{
    for(size_t i = 0;i<boundary_ids.size();++i) {
      if(BC1name == boundary_ids[i].second)
        BC1_id = boundary_ids[i].first ;
      if(BC2name == boundary_ids[i].second)
        BC2_id = boundary_ids[i].first ;
    }
  }
  if(BC1_id == -1) 
    cerr << "unable to find boundary name '" << BC1name << "' in grid." ;
  if(BC2_id == -1) 
    cerr << "unable to find boundary name '" << BC2name << "' in grid." ;
  if(BC1_id == -1 || BC2_id == -1)
    exit(-1) ;
  
  if(verbose) {
    cout << "Reading grid file: " << gridfile << "... " ;
    cout.flush() ;
  }

  gettimeofday(&time_grid_read_start,NULL) ;
 
  vector<entitySet> local_nodes,local_faces,local_cells ;
  store<vector3d<double> > pos ;
  Map cl,cr ;
  multiMap face2node ;
  int max_alloc=0 ;
  vector<pair<string,entitySet> > volTags ;
  
  store<string> boundary_names ;
  store<string> boundary_tags ;

  if(!Loci::readGridVOG(local_nodes,local_faces,local_cells,
                        pos,cl,cr,face2node,
			boundary_names,boundary_tags,volTags,
			max_alloc,gridfile)) {
                        
    if(Loci::MPI_rank == 0) {
      cerr << endl
           << "Reading grid file '" << gridfile <<"' failed in grid reader!"
           << endl ;
      Loci::Abort() ;
    }
  }
  if(verbose) {
    cout << "Done" << endl ;
  }

  gettimeofday(&time_grid_read_end,NULL) ;

  gettimeofday(&time_essential_start,NULL) ;
  
  // NOTE this has been revised but not tested after changes to
  // VOG file reader
  entitySet dom = boundary_names.domain() ;
  BC1_id = -1 ;
  BC2_id = -1 ;
  FORALL(dom,bc) {
    if(BC1name == boundary_names[bc])
      BC1_id = bc ;
    if(BC2name == boundary_names[bc])
      BC2_id = bc ;
  } ENDFORALL ;

  entitySet BC1_faces, BC2_faces;
  dom = cr.domain() ;
  FORALL(dom,fc) {
    if(cr[fc] == BC1_id)
      BC1_faces += fc ;
    if(cr[fc] == BC2_id)
      BC2_faces += fc ;
  } ENDFORALL ;
   
  
  gettimeofday(&time_data_structure_start,NULL) ;
  
  multiMap node2face ;
  Loci::inverseMap(node2face,face2node,pos.domain(),face2node.domain()) ;
  
  
  entitySet BC1_nodes = get_boundary_nodes(face2node,BC1_faces) ;
  entitySet BC2_nodes = get_boundary_nodes(face2node,BC2_faces) ;

  store<vec3d> BC1_nodes_pos(get_nodes_pos(pos,BC1_nodes)) ;
  store<vec3d> BC2_nodes_pos(get_nodes_pos(pos,BC2_nodes)) ;

  // then we get the edges in two boundary meshes
  // BC*_N1 and BC*_N2 defines the end nodes of each edge
  // on boundary BC*. BC*_El and BC*_Er defines the left and
  // right faces of the edge. BC*_edge_num hold the total
  // edge number for BC*
  dMap BC1_N1, BC1_N2 ; dMap BC1_El, BC1_Er ; int BC1_edge_num ;
  get_edges_info(face2node,node2face,BC1_faces,cr,
                 BC1_N1,BC1_N2,BC1_El,BC1_Er,BC1_edge_num) ;
  dMap BC2_N1, BC2_N2 ; dMap BC2_El, BC2_Er ; int BC2_edge_num ;
  get_edges_info(face2node,node2face,BC2_faces,cr,
                 BC2_N1,BC2_N2,BC2_El,BC2_Er,BC2_edge_num) ;

  gettimeofday(&time_data_structure_end,NULL) ;
  data_structure_time +=
    time_data_structure_end - time_data_structure_start ;

  gettimeofday(&time_essential_end,NULL) ;
  prog_essential_time +=
    time_essential_end - time_essential_start ;

  if(verbose) {
    cout << "--------" << endl ;  
    cout << BC1name <<" faces number: " << BC1_faces.size() << endl ;
    cout << BC1name <<" nodes number: " << BC1_nodes.size() << endl ;
    cout << BC1name <<" edges number: " << BC1_edge_num << endl ;
    cout << "--------" << endl ;  
    cout << BC2name <<" faces number: " << BC2_faces.size() << endl ;
    cout << BC2name <<" nodes number: " << BC2_nodes.size() << endl ;
    cout << BC2name <<" edges number: " << BC2_edge_num << endl ;
  }

  if(verbose) {
    cout << "--------" << endl ;  
    cout << "Projecting boundaries... " ;
  }
  
  gettimeofday(&time_essential_start,NULL) ;
  gettimeofday(&time_3D_to_2D_projection_start,NULL) ;
  // then we do the projection for the two boundary meshes
  store<vec2d> BC1_proj_pos, BC2_proj_pos ;
  // first we compute the projection vectors
  if(is_rotate) {
    normalize(rotation_axis) ;
    vec3d u,v ;
    orthogonal_coords(rotation_axis,u,v) ;

    axis_radius_projection(rotation_axis,u,v,rotation_center,
                           BC1_nodes_pos,BC1_proj_pos) ;
    axis_radius_projection(rotation_axis,u,v,rotation_center,
                           BC2_nodes_pos,BC2_proj_pos) ;
  }else if(is_translate) {
    vec3d translation_vector_cp = translation_vector ;
    // we cannot normalize translation_vector because
    // we still need it in the following code for 2D->3D projection
    normalize(translation_vector_cp) ;
    vec3d u,v ;
    orthogonal_coords(translation_vector_cp,u,v) ;
    orthogonal_projection(u,v,BC1_nodes_pos,BC1_proj_pos) ;
    orthogonal_projection(u,v,BC2_nodes_pos,BC2_proj_pos) ;    
  }
  gettimeofday(&time_3D_to_2D_projection_end,NULL) ;
  gettimeofday(&time_essential_end,NULL) ;
  prog_essential_time +=
    time_essential_end - time_essential_start ;
  
  if(verbose)
    cout << "Done" << endl ;

  if(twoD_aux_viz) {
    string BC1_out_name = problem_name + "BC" +
      BC1name + ".2dgv" ;
    string BC2_out_name = problem_name + "BC" +
      BC2name + ".2dgv" ;
    if(verbose) {
      cout << "--------" << endl ;  
      cout << "Writing out "<<BC1name<<" projection into: "
           << BC1_out_name << "... " ;
    }
    TwodgvOutput(BC1_out_name.c_str(),BC1_proj_pos,BC1_N1,BC1_N2,
                 BC1_El,BC1_Er,BC1_edge_num, BC1_faces) ;
    if(verbose)
      cout << "Done" << endl ;

    if(verbose) {
      cout << "--------" << endl ;  
      cout << "Writing out "<<BC2name<<" projection into: "
           << BC2_out_name << "... " ;
    }
    TwodgvOutput(BC2_out_name.c_str(),BC2_proj_pos,BC2_N1,BC2_N2,
                 BC2_El,BC2_Er,BC2_edge_num, BC2_faces) ;
    if(verbose)
      cout << "Done" << endl ;
    
    string two_out = problem_name + "BC" +
      BC1name + BC2name + ".2dgv" ;

    if(verbose) {
      cout << "--------" << endl ;  
      cout << "Writing out both projections into: "
           << two_out << "... " ;
    }
    TwodgvOutputBOTH(two_out.c_str(),
                     BC1_proj_pos,BC2_proj_pos,
                     BC1_N1,BC1_N2,BC1_El,BC1_Er,
                     BC2_N1,BC2_N2,BC2_El,BC2_Er,
                     BC1_edge_num,BC2_edge_num,
                     BC1_faces,BC2_faces) ;
    if(verbose)
      cout << "Done" << endl ;
  } // end of if(twoD_aux_viz)
  
  gettimeofday(&time_essential_start,NULL) ;
  gettimeofday(&time_data_structure_start,NULL) ;
  // we get the node2edge relations
  multiMap BC1_node2edge ;
  get_node2edge(BC1_N1,BC1_N2,BC1_node2edge) ;
  multiMap BC2_node2edge ;
  get_node2edge(BC2_N1,BC2_N2,BC2_node2edge) ;
//   // we also need the node2node relations
//   multiMap BC1_node2node ;
//   get_node2node(BC1_nodes,BC1_N1,BC1_N2,BC1_node2edge,BC1_node2node) ;
//   multiMap BC2_node2node ;
//   get_node2node(BC2_nodes,BC2_N1,BC2_N2,BC2_node2edge,BC2_node2node) ;
  gettimeofday(&time_data_structure_end,NULL) ;
  data_structure_time +=
    time_data_structure_end - time_data_structure_start ;
  gettimeofday(&time_essential_end,NULL) ;
  prog_essential_time +=
    time_essential_end - time_essential_start ;

  // shifting points first
  if(verbose) {
    cout << "--------" << endl ;
    cout << "Shifting points... " ;
  }
  gettimeofday(&time_essential_start,NULL) ;
  gettimeofday(&time_data_structure_start,NULL) ;
  entitySet BC1_edge_nodes = get_edge_nodes(BC1_N1,BC1_N2,BC1_Er) ;
  entitySet BC2_edge_nodes = get_edge_nodes(BC2_N1,BC2_N2,BC2_Er) ;

#ifdef POINT_SHIFT_CONVEXITY_CHECK
  masterbc_edge_nodes = BC1_edge_nodes ;
  slavebc_edge_nodes = BC2_edge_nodes ;
#endif

  gettimeofday(&time_data_structure_end,NULL) ;
  data_structure_time +=
    time_data_structure_end - time_data_structure_start ;
  
  dMap BC1p_2_BC2p, BC2p_2_BC1p ;
  
  gettimeofday(&time_shift_point_total_start,NULL) ;
  entitySet shifted = 
    shift_points2(BC1_faces,BC2_faces,
                  BC1_nodes,BC1_proj_pos,
                  BC2_nodes,BC2_proj_pos,
                  face2node,node2face,
                  BC1p_2_BC2p,BC2p_2_BC1p) ;
  gettimeofday(&time_shift_point_total_end,NULL) ;
  gettimeofday(&time_essential_end,NULL) ;
  prog_essential_time +=
    time_essential_end - time_essential_start ;
  
  if(verbose) {
    cout << "[" << shifted.size() << " points shifted on "
         << BC2name << "]" << endl ;
  }

  if(twoD_aux_viz) {
    string shift_bc2_out = problem_name + "BC" +
      BC2name + "-shift.2dgv" ;
    if(verbose) {
      cout << "--------" << endl ;  
      cout << "Writing out shifted "<<BC2name<<" projection into: "
           << shift_bc2_out << "... " ;
    }
    TwodgvOutput(shift_bc2_out.c_str(),BC2_proj_pos,BC2_N1,BC2_N2,
                 BC2_El,BC2_Er,BC2_edge_num, BC2_faces) ;
    if(verbose)
      cout << "Done" << endl ;
    
    string shift_two_out = problem_name + "BC"+
      BC1name+ BC2name + "-shift.2dgv" ;
    if(verbose) {
      cout << "--------" << endl ;  
      cout << "Writing out shifted projection (overlapped) into: "
           << shift_two_out << "... " ;
    }
    TwodgvOutputBOTH(shift_two_out.c_str(),
                     BC1_proj_pos,BC2_proj_pos,
                     BC1_N1,BC1_N2,BC1_El,BC1_Er,
                     BC2_N1,BC2_N2,BC2_El,BC2_Er,
                     BC1_edge_num,BC2_edge_num,
                     BC1_faces,BC2_faces) ;
    if(verbose)
      cout << "Done" << endl ;
  }
  
  if(verbose) {
    cout << "-------" << endl ;
    cout << "Computing face splits..." << endl ;
  }

  gettimeofday(&time_essential_start,NULL) ;
  store<vec2d> BC1_new_nodes ;
  multiMap BC1_new_face2node ;
  multiMap BC1_new_interior_face2node ;
  dstore<alphaBeta> BC1_inverse_proj_info ;
  dMap BC1_new_cl, BC1_new_cr ;
  // we will need to first get the edge->interior faces map
  multiMap BC1_edge2iface ;

  gettimeofday(&time_data_structure_start,NULL) ;
  
  get_edge2interiorFace(BC1_N1,BC1_N2,BC1_El,BC1_Er,
                        node2face,BC1_edge2iface) ;

  gettimeofday(&time_data_structure_end,NULL) ;
  data_structure_time +=
    time_data_structure_end - time_data_structure_start ;
  
  // we also need to compute the maximum node index number
  Entity new_node_index = node2face.domain().Max() + 1 ;
  Entity new_face_index = face2node.domain().Max() + 1 ;
  // first we set the global mesh_new_node_index
  mesh_new_node_index = new_node_index ;
  mesh_new_face_index = new_face_index ;

  gettimeofday(&time_face_split_total_start,NULL) ;
  int total_split_pairs =
    face_split(BC1_faces, BC2_faces, BC1_nodes, BC1_proj_pos,
               BC2_nodes, BC2_proj_pos, face2node, node2face,
               BC1_edge_nodes, BC2_edge_nodes, BC1p_2_BC2p) ;
  gettimeofday(&time_face_split_total_end,NULL) ;
  gettimeofday(&time_essential_end,NULL) ;
  prog_essential_time +=
    time_essential_end - time_essential_start ;
  
  if(verbose) {
    cout << "  total pairs of face split: "
         << total_split_pairs << endl ;
  }

  if(verbose) {
    cout << "--------" << endl ;
    cout << "Reconstructing " << BC1name << " topology..." ;
    cout.flush() ;
  }

  gettimeofday(&time_essential_start,NULL) ;
  gettimeofday(&time_BC1_reconstruct_start,NULL) ;
  entitySet BC1_faces_split =
    gen_boundary_topo(global_bc1_split_records,
                      BC1_nodes, BC1_proj_pos,
                      BC2_nodes, BC2_proj_pos,
                      face2node,
                      BC1_node2edge, BC1_edge2iface,
                      cl, cr,
                      // these are the outputs stores
                      BC1_new_nodes, BC1_new_face2node,
                      BC1_new_interior_face2node,
                      BC1_new_cl, BC1_new_cr,
                      BC1_inverse_proj_info) ;
  gettimeofday(&time_BC1_reconstruct_end,NULL) ;
  gettimeofday(&time_essential_end,NULL) ;
  prog_essential_time +=
    time_essential_end - time_essential_start ;
  
  if(verbose) {
    cout << " Done" << endl ;
  }

  gettimeofday(&time_essential_start,NULL) ;
  gettimeofday(&time_data_structure_start,NULL) ;
  // update face2node to include the change of interior faces
  face2node = merge_multiMap(face2node,BC1_new_interior_face2node) ;
  
  gettimeofday(&time_data_structure_end,NULL) ;
  data_structure_time +=
    time_data_structure_end - time_data_structure_start ;
  // compute new BC1 faces set
  entitySet BC1_new_faces =
    (BC1_faces - BC1_faces_split) + BC1_new_face2node.domain() ;
  
  gettimeofday(&time_essential_end,NULL) ;
  prog_essential_time +=
    time_essential_end - time_essential_start ;
    
  if(verbose) {
    cout << "  new added nodes number: "
         << BC1_new_nodes.domain().size() << endl ;
    cout << "  new total faces number: " << BC1_new_faces.size()
         << endl ;
    cout << endl ;
  }
  
  if(twoD_aux_viz) {
    string split_bc1_poly = problem_name + "BC" + BC1name +
      "-new.poly" ;
    string split_bc1_2dgv = problem_name + "BC" + BC1name +
      "-new.2dgv" ;
    if(verbose) {
      cout << "  Writing out face split results into: " << endl
           << '\t' << split_bc1_poly << " (for \"showme\" use) and"
           << endl
           << '\t' << split_bc1_2dgv << " (for \"2dgv\" use)"
           << endl ;
      cout << endl ;
    }
    viz_split_showme_2dgv(BC1_nodes,BC1_proj_pos,
                          BC1_new_nodes,BC1_new_face2node,
                          split_bc1_poly.c_str(),
                          split_bc1_2dgv.c_str()) ;
  }

  if(after_split_check) {
    if(verbose)
      cout << "  Checking new topology validity..." << endl ;

    if(verbose)
      cout << "  --Checking split face definition..." ;
    if(!face_split_polygon_check(BC1_new_face2node)) {
      cout << endl ;
      cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++"
           << endl ;
      cout << "polygon check FAILED!! Please report when this happens"
           << endl ;
      cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++"
           << endl ;
      Loci::Abort() ;
    }else {
      if(verbose)
        cout << "\tPASSED!" << endl ;
    }

    if(verbose)
      cout << "  --Checking edge map..." ;

    gettimeofday(&time_essential_start,NULL) ;
    // get inverse map new_node2face
    multiMap BC1_new_node2face ;
    gettimeofday(&time_data_structure_start,NULL) ;

    get_node2face(BC1_new_face2node,BC1_new_node2face) ;
  
    gettimeofday(&time_data_structure_end,NULL) ;
    data_structure_time +=
      time_data_structure_end - time_data_structure_start ;
    gettimeofday(&time_essential_end,NULL) ;
    prog_essential_time +=
      time_essential_end - time_essential_start ;
  
    if(!face_split_edgeMap_check(BC1_new_face2node,BC1_new_node2face)) {
      cout << endl ;
      cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++"
           << endl ;
      cout << "edge map check FAILED!! Please report when this happens"
           << endl ;
      cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++"
           << endl ;
      Loci::Abort() ;
    }else {
      if(verbose)
        cout << "\t\tPASSED!" << endl ;
    }

    if(verbose)
      cout << "  --Checking area sum..." ;
    if(!face_split_areaSum_check(BC1_faces_split,BC1_proj_pos,face2node,
                                 BC1_new_nodes,BC1_new_face2node)) {
      cout << endl ;
      cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++"
           << endl ;
      cout << "area sum check FAILED!! Please report when this happens"
           << endl ;
      cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++"
           << endl ;
      Loci::Abort() ;
    }else {
      if(verbose)
        cout << "\t\tPASSED!" << endl ;
    }

#ifdef FACE_SPLIT_PROPERTY_CHECK
    cout << "  --Checking face convexity..." ;
    {
      double smallest_area ;
      int small_faces ;
      bool all_convex =
        face_split_convexity_area_scan(global_bc1_split_records,
                                       smallest_area,small_faces) ;
      cout << "\t\tDone" << endl ;
      cout << "    --smallest face area:         " << smallest_area << endl ;
      cout << "    --small faces (area < 1e-12): " << small_faces << endl ;
      cout << "    --all faces convex?:          "
           << (all_convex?"true":"false") << endl ;
    }
#endif
  } // end of if(after_split_check)
  
  gettimeofday(&time_essential_start,NULL) ;
  store<vec2d> BC2_new_nodes ;
  multiMap BC2_new_face2node ;
  multiMap BC2_new_interior_face2node ;
  dstore<alphaBeta> BC2_inverse_proj_info ;
  dMap BC2_new_cl, BC2_new_cr ;
  // we will need to first get the edge -> interior face map
  multiMap BC2_edge2iface ;
  gettimeofday(&time_data_structure_start,NULL) ;

  get_edge2interiorFace(BC2_N1,BC2_N2,BC2_El,BC2_Er,
                        node2face,BC2_edge2iface) ;

  gettimeofday(&time_data_structure_end,NULL) ;
  data_structure_time +=
    time_data_structure_end - time_data_structure_start ;
  gettimeofday(&time_essential_end,NULL) ;
  prog_essential_time +=
    time_essential_end - time_essential_start ;
  
  if(verbose) {
    cout << "--------" << endl ;
    cout << "Reconstructing " << BC2name << " topology..." ;
    cout.flush() ;
  }

  gettimeofday(&time_essential_start,NULL) ;
  gettimeofday(&time_BC2_reconstruct_start,NULL) ;
  entitySet BC2_faces_split =
    gen_boundary_topo(global_bc2_split_records,
                      BC2_nodes, BC2_proj_pos,
                      BC1_nodes, BC1_proj_pos,
                      face2node,
                      BC2_node2edge, BC2_edge2iface,
                      cl, cr,
                      // these are the outputs stores
                      BC2_new_nodes, BC2_new_face2node,
                      BC2_new_interior_face2node,
                      BC2_new_cl, BC2_new_cr,
                      BC2_inverse_proj_info) ;
  gettimeofday(&time_BC2_reconstruct_end,NULL) ;
  gettimeofday(&time_essential_end,NULL) ;
  prog_essential_time +=
    time_essential_end - time_essential_start ;
  
  if(verbose) {
    cout << " Done" << endl ;
  }
  // update face2node to include the change of interior faces
  gettimeofday(&time_essential_start,NULL) ;
  gettimeofday(&time_data_structure_start,NULL) ;
  
  face2node = merge_multiMap(face2node,BC2_new_interior_face2node) ;
  
  gettimeofday(&time_data_structure_end,NULL) ;
  data_structure_time +=
    time_data_structure_end - time_data_structure_start ;
  // compute new BC2 faces set
  entitySet BC2_new_faces =
    (BC2_faces - BC2_faces_split) + BC2_new_face2node.domain() ;
  gettimeofday(&time_essential_end,NULL) ;
  prog_essential_time +=
    time_essential_end - time_essential_start ;

  if(verbose) {
    cout << "  new added nodes number: "
         << BC2_new_nodes.domain().size() << endl ;
    cout << "  new total faces number: " << BC2_new_faces.size()
         << endl ;
    cout << endl ;
  }
  
  if(twoD_aux_viz) {
    string split_bc2_poly = problem_name + "BC" + BC2name +
      "-new.poly" ;
    string split_bc2_2dgv = problem_name + "BC" + BC2name +
      "-new.2dgv" ;
    if(verbose) {
      cout << "  Writing out face split results into: " << endl
           << '\t' << split_bc2_poly << " (for \"showme\" use) and"
           << endl
           << '\t' << split_bc2_2dgv << " (for \"2dgv\" use)"
           << endl ;
      cout << endl ;
    }
    viz_split_showme_2dgv(BC2_nodes,BC2_proj_pos,
                          BC2_new_nodes,BC2_new_face2node,
                          split_bc2_poly.c_str(),
                          split_bc2_2dgv.c_str()) ;
  }

  if(after_split_check) {
    if(verbose)
      cout << "  Checking new topology validity..." << endl ;
    
    if(verbose)
      cout << "  --Checking split face definition..." ;
    if(!face_split_polygon_check(BC2_new_face2node)) {
      cout << endl ;
      cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++"
           << endl ;
      cout << "polygon check FAILED!! Please report when this happens"
           << endl ;
      cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++"
           << endl ;
      Loci::Abort() ;
    }else {
      if(verbose)
        cout << "\tPASSED!" << endl ;
    }
    
    if(verbose)
      cout << "  --Checking edge map..." ;
    // get inverse map new_node2face
    gettimeofday(&time_essential_start,NULL) ;
    multiMap BC2_new_node2face ;
    gettimeofday(&time_data_structure_start,NULL) ;
    
    get_node2face(BC2_new_face2node,BC2_new_node2face) ;
    
    gettimeofday(&time_data_structure_end,NULL) ;
    data_structure_time +=
      time_data_structure_end - time_data_structure_start ;
    gettimeofday(&time_essential_end,NULL) ;
    prog_essential_time +=
      time_essential_end - time_essential_start ;
    
    if(!face_split_edgeMap_check(BC2_new_face2node,BC2_new_node2face)) {
      cout << endl ;
      cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++"
           << endl ;
      cout << "edge map check FAILED!! Please report when this happens"
           << endl ;
      cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++"
           << endl ;
      Loci::Abort() ;
    }else {
      if(verbose)
        cout << "\t\tPASSED!" << endl ;
    }
    
    if(verbose)
      cout << "  --Checking area sum..." ;
    if(!face_split_areaSum_check(BC2_faces_split,BC2_proj_pos,face2node,
                                 BC2_new_nodes,BC2_new_face2node)) {
      cout << endl ;
      cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++"
           << endl ;
      cout << "area sum check FAILED!! Please report when this happens"
           << endl ;
      cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++"
           << endl ;
      Loci::Abort() ;
    }else {
      if(verbose)
        cout << "\t\tPASSED!" << endl ;
    }
    
#ifdef FACE_SPLIT_PROPERTY_CHECK
    cout << "  --Checking face convexity..." ;
    {
      double smallest_area ;
      int small_faces ;
      bool all_convex =
        face_split_convexity_area_scan(global_bc2_split_records,
                                       smallest_area,small_faces) ;
      cout << "\t\tDone" << endl ;
      cout << "    --smallest face area:         " << smallest_area << endl ;
      cout << "    --small faces (area < 1e-12): " << small_faces << endl ;
      cout << "    --all faces convex?:          "
           << (all_convex?"true":"false") << endl ;
    }
#endif
  } // end of if(after_split_check)
  
  if(verbose) {
    cout << "--------" << endl ;
    cout << "Matching points on boundaries... " ;
    cout.flush() ;
  }
  
  gettimeofday(&time_essential_start,NULL) ;
  entitySet bc1_not_matched, bc2_not_matched ;
  dMap BC2_2_BC1_node_map ;

  gettimeofday(&time_boundary_nodes_matching_start,NULL) ;
  twoD_point_match(BC1_proj_pos, BC2_proj_pos,
                   BC1_new_nodes, BC2_new_nodes,
                   BC1_new_faces,face2node,
                   BC1_new_face2node, BC2_2_BC1_node_map,
                   bc1_not_matched, bc2_not_matched) ;
  gettimeofday(&time_boundary_nodes_matching_end, NULL) ;
  gettimeofday(&time_essential_end,NULL) ;
  prog_essential_time +=
    time_essential_end - time_essential_start ;
  
  if( (bc1_not_matched.size() != 0) ||
      (bc2_not_matched.size() != 0)) {
    cout << "Error: boundary points not matched!" << endl ;
    cout << "       please report this problem!" << endl ;
    Loci::Abort() ;
  }
  if(verbose) {
    cout << "Done (PASSED!)" << endl ;
  }
  
  if(verbose) {
    cout << "--------" << endl ;
    cout << "Projecting boundaries to 3D space... " ;
    cout.flush() ;
  }

  gettimeofday(&time_essential_start,NULL) ;
  store<vec3d> BC1_new_nodes_3D, BC2_new_nodes_3D ;
  rigid_transform rt(rotation_center,rotation_axis,
                     deg2rad(rotation_angle_in_degree),
                     translation_vector) ;

  gettimeofday(&time_2D_to_3D_projection_start,NULL) ;
  get_3D_nodes_pos(BC1_inverse_proj_info,BC2_2_BC1_node_map,rt,
                   BC1_nodes,BC1_inverse_proj_info.domain(),
                   BC2_nodes,BC2_inverse_proj_info.domain(),
                   pos,BC1_new_nodes_3D,BC2_new_nodes_3D) ;
  gettimeofday(&time_2D_to_3D_projection_end, NULL) ;
  gettimeofday(&time_essential_end,NULL) ;
  gettimeofday(&time_essential_end,NULL) ;
  prog_essential_time +=
    time_essential_end - time_essential_start ;
  
  if(verbose) {
    cout << "Done" << endl ;
  }

  if(output_name == "") {
    output_name = "pb_"+gridfile ;
  }
  if(output_name.size() < 5 ||
     output_name.substr(output_name.size()-4,4) != ".vog") {
    output_name+= ".vog" ;
  }
  if(output_name == gridfile)
    output_name = "pb_"+gridfile ;
    
  // we are now ready to write out the new  grid
  if(verbose) {
    cout << "--------" << endl ;
    cout << "Generating new grid: " << output_name << endl ;
  }

  gettimeofday(&time_new_grid_write_start,NULL) ;
  
  //outfile : output_name
  
  int npoints, nfaces, ncells ;
  
  npoints = pos.domain().size() + BC1_new_nodes.domain().size()
    + BC2_new_nodes.domain().size() ;

  nfaces = face2node.domain().size() -
    BC1_faces_split.size() +
    BC1_new_face2node.domain().size() -
    BC2_faces_split.size() +
    BC2_new_face2node.domain().size() ;

  entitySet range ;
  entitySet domain = cl.domain() ;
  range += cl.image(domain) ;
  domain = cr.domain() ;
  range += cr.image(domain) ;
  
  range &= interval(0,Loci::UNIVERSE_MAX) ;
  ncells = range.size() ;
  
  
  if(verbose) {
    cout << "  total number of points:               " << npoints << endl ;
    cout << "  total number of faces:                " << nfaces << endl ;
    cout << "  total number of cells:                " << ncells << endl ;
  }

  // we will need to re-number the cells
  // cells are numbered from 1
  store<int> cell_index ;
  cell_index.allocate(range) ;
  int count = 1 ;
  for(entitySet::const_iterator ei=range.begin();
      ei!=range.end();++ei) {
    cell_index[*ei] = count ;
    ++count ;
  }
  
  //define the new pos
  store<vector3d<double> > newPos;
  entitySet newPosDom = interval(0, npoints-1);
  newPos.allocate(newPosDom);

  // we will need to re-number the nodes
  // nodes are numbered from 0
  store<int> node_index ;
  entitySet all_nodes_dom = pos.domain() |
    BC1_new_nodes_3D.domain() | BC2_new_nodes_3D.domain() ;
  node_index.allocate(all_nodes_dom) ;
  count = 0 ;
  // we write out the nodes
  domain = pos.domain() ;
  entitySet::const_iterator ni = newPosDom.begin();
  for(entitySet::const_iterator ei=domain.begin();
      ei!=domain.end();++ei,++count, ni++) {
    // set index first
    node_index[*ei] = count ;
    newPos[*ni] = pos[*ei];
    // then write out
  }
  domain = BC1_new_nodes_3D.domain() ;
  for(entitySet::const_iterator ei=domain.begin();
      ei!=domain.end();++ei,++count, ni++) {
    // set index first
    node_index[*ei] = count ;
    newPos[*ni] =BC1_new_nodes_3D[*ei];
    // then write out
  }
  domain = BC2_new_nodes_3D.domain() ;
  for(entitySet::const_iterator ei=domain.begin();
      ei!=domain.end();++ei,++count, ni++) {
    // set index first
    node_index[*ei] = count ;
    newPos[*ni] =BC2_new_nodes_3D[*ei];
  }

  // we then write out face offset and cl, cr
  // first write out faces other than BC1 and BC2 new faces
  // and some of the potential changed boundary faces
  entitySet newFaceDom = interval(0, nfaces-1);
  Map newCl;
  Map newCr;
  store<int> faceCount;
  multiMap newFace2node;
  
  newCl.allocate(newFaceDom);
  newCr.allocate(newFaceDom);
  faceCount.allocate(newFaceDom);
  
  entitySet::const_iterator fi = newFaceDom.begin();
  
  domain = face2node.domain() - BC1_faces_split -
    BC2_faces_split ;
  
   count = 0 ;
  for(entitySet::const_iterator ei=domain.begin();
      ei!=domain.end();++ei, fi++) {
    int offset = face2node.num_elems(*ei) ;
    faceCount[*fi] = offset;
    Entity c1 = cl[*ei] ; if(c1>=0) c1=cell_index[c1] ;
    Entity c2 = cr[*ei] ; if(c2>=0) c2=cell_index[c2] ;
    newCl[*fi] = c1;
    newCr[*fi] = c2;
    count += offset ; 
  }
  // then BC1 new faces
  domain = BC1_new_face2node.domain() ;
  for(entitySet::const_iterator ei=domain.begin();
      ei!=domain.end();++ei, fi++) {
    int offset = BC1_new_face2node.num_elems(*ei) ;
    faceCount[*fi] = offset;
    Entity c1 = BC1_new_cl[*ei] ; if(c1>=0) c1=cell_index[c1] ;
    Entity c2 = BC1_new_cr[*ei] ; if(c2>=0) c2=cell_index[c2] ;
    newCl[*fi] = c1;
    newCr[*fi] = c2;
    count += offset ; 
  }
  // then BC2 new faces
  domain = BC2_new_face2node.domain() ;
  for(entitySet::const_iterator ei=domain.begin();
      ei!=domain.end();++ei, fi++) {
    int offset = BC2_new_face2node.num_elems(*ei) ;
     faceCount[*fi] = offset;
    Entity c1 = BC2_new_cl[*ei] ; if(c1>=0) c1=cell_index[c1] ;
    Entity c2 = BC2_new_cr[*ei] ; if(c2>=0) c2=cell_index[c2] ;
    newCl[*fi] = c1;
    newCr[*fi] = c2;
    count += offset ;
  }
  
  
  if(verbose) {
    cout << "  total number of face defining points: "
         << count << endl ;
  }
  newFace2node.allocate(faceCount);
  

  fi = newFaceDom.begin();
  // then we write out face2node array
  // first write out faces other than BC1 and BC2 new faces
  // and some of the potential changed boundary faces
  domain = face2node.domain() - BC1_faces_split -
    BC2_faces_split ;
  for(entitySet::const_iterator ei=domain.begin();
      ei!=domain.end();++ei, fi++) {
    vector<Entity> nodes = get_face_nodes(face2node,*ei) ;
    for(unsigned int i = 0; i < nodes.size(); i++)
      newFace2node[*fi][i] = node_index[nodes[i]] ;
  }
  // then BC1 new faces
  domain = BC1_new_face2node.domain() ;
  for(entitySet::const_iterator ei=domain.begin();
      ei!=domain.end();++ei, fi++) {
    vector<Entity> nodes = get_face_nodes(BC1_new_face2node,*ei) ;
    for(unsigned int i = 0; i < nodes.size(); i++)
      newFace2node[*fi][i] = node_index[nodes[i]] ;
  }
  // then BC2 new faces
  domain = BC2_new_face2node.domain() ;
  for(entitySet::const_iterator ei=domain.begin();
      ei!=domain.end();++ei, fi++) {
    vector<Entity> nodes = get_face_nodes(BC2_new_face2node,*ei) ;
    for(unsigned int i = 0; i < nodes.size(); i++)
      newFace2node[*fi][i] = node_index[nodes[i]] ;
  }
  
  // AND WE ARE FINALLY DONE
  // establish face left-right orientation
  if(verbose)
    cerr << "  --orienting faces" << endl ;
  VOG::orientFaces(newPos,newCl,newCr,newFace2node) ;
  
  if(verbose)
    cerr << "  --coloring matrix" << endl ;
  VOG::colorMatrix(newPos,newCl,newCr,newFace2node) ;
 
  if(optimize) {
    if(verbose)
      cerr << "  --optimizing mesh layout" << endl ;
    VOG::optimizeMesh(newPos,newCl,newCr,newFace2node) ;
  }
  
  if(verbose)
    cerr << "  --writing VOG file: " << output_name << endl ;
  
  //get boundary names
  vector<pair<int,string> > surf_ids = boundary_ids ;
  Loci::writeVOG(output_name, newPos, newCl, newCr, newFace2node,surf_ids) ;

  gettimeofday(&time_new_grid_write_end,NULL) ;
  gettimeofday(&time_prog_end,NULL) ;

  // generate a timing report
  if(report_timing) {
    cout << "--------" << endl ;
    cout << "timing report (all in seconds): " << endl ;
    cout << "total time used:                "
         << time_prog_end - time_prog_start << endl ;
    cout << "pure time (no IO and check):    "
         << prog_essential_time << endl ;

    cout << "  grid reading:                 "
         << time_grid_read_end - time_grid_read_start << endl ;

    cout << "  topology data building:       "
         << data_structure_time << endl ;

    cout << "  point shifting (total):       "
         << time_shift_point_total_end - time_shift_point_total_start
         << endl ;
    cout << "    quad tree building:         " << "  "
         << time_shift_point_qt_end - time_shift_point_qt_start << endl ;
    cout << "    point matching:             " << "  "
         << time_shift_point_end - time_shift_point_start << endl ;
    
    cout << "  face splitting (total):       "
         << time_face_split_total_end - time_face_split_total_start
         << endl ;
    cout << "    overlap face computing:     " << "  "
         << time_face_split_frm_end - time_face_split_frm_start << endl ;
    cout << "    quad tree building:         " << "  "
         << time_face_split_qt_end - time_face_split_qt_start << endl ;
    cout << "    face splitting:             " << "  "
         << time_face_split_end - time_face_split_start << endl ;

    cout << "  boundary 1 reconstruction:    "
         << time_BC1_reconstruct_end - time_BC1_reconstruct_start << endl ;

    cout << "  boundary 2 reconstruction:    "
         << time_BC2_reconstruct_end - time_BC2_reconstruct_start << endl ;

    cout << "  boundary nodes matching:      "
         << time_boundary_nodes_matching_end -
      time_boundary_nodes_matching_start << endl ;
    
    cout << "  3D -> 2D projection:          "
         << time_3D_to_2D_projection_end - time_3D_to_2D_projection_start
         << endl ;
    cout << "  2D -> 3D projection:          "
         << time_2D_to_3D_projection_end - time_2D_to_3D_projection_start
         << endl ;

    cout << "  new grid generation:          "
         << time_new_grid_write_end - time_new_grid_write_start << endl ;
  }
  
  if(verbose)
    cout << "..ALL DONE.." << endl ;

#endif // matching #ifdef NORMAL

  Loci::Finalize() ;
}
// The End //
