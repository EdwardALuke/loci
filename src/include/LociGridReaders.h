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
#ifndef LOCI_GRID_READERS_H
#define LOCI_GRID_READERS_H
#include <fact_db.h>
#include <Tools/basic_types.h>
#include <string>

namespace Loci {

  /// Define struct Area which data members of normal vector and area of the area
  template <class REAL> struct Area_t {
    vector3d<REAL> n ;  //normal vector of the face
    REAL sada ; //area of the face
  } ;

  /// Overload ostream (output) operator for struct Area
  template<class REAL>
  inline std::ostream & operator<<(std::ostream &s, const Area_t<REAL> &v)
  {
    s << v.n << ' ' << v.sada << ' ' ;
    return s ;
  }

  /// Overload istream (input) operator for struct Area
  template <class REAL>
  inline std::istream &operator>>(std::istream &s, Area_t<REAL> &v)
  {
    s >> v.n >> v.sada  ;
    return s ;
  }

  typedef Area_t<real_t> Area ;

  template<class T> struct data_schema_traits<Loci::Area_t<T> > {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(Loci::Area_t<T>()) ;
      LOCI_INSERT_TYPE(ct,Loci::Area_t<T>,n) ;
      LOCI_INSERT_TYPE(ct,Loci::Area_t<T>,sada) ;
      return DatatypeP(ct) ;
    }
  } ;


  /// Defines a rigid body transformation consisting of a rotation about an
  /// axis and a translation: v' = R * (v + t1) + t2
  struct rigid_transform {
    vector3d<double> t1,t2 ;
    tensor3d<double> R,Rinv ;
    rigid_transform() {
      t1 = vector3d<double>(0,0,0) ;
      t2 = t1 ;
      R.x = vector3d<double>(1,0,0) ;
      R.y = vector3d<double>(0,1,0) ;
      R.z = vector3d<double>(0,0,1) ;
      Rinv = R ;
    }
    rigid_transform(vector3d<double> center, vector3d<double> v, double angle, vector3d<double> translate) {
      t1 = -1.*center ;
      t2 = center + translate ;
      R.x = vector3d<double>(1,0,0) ;
      R.y = vector3d<double>(0,1,0) ;
      R.z = vector3d<double>(0,0,1) ;
      Rinv = R ;
      if(angle == 0)
        return ;
      double s = sin(angle) ;
      double c = cos(angle) ;
      double C = 1.-c ;
      double x = v.x ;
      double y = v.y ;
      double z = v.z ;
      R.x += vector3d<double>(C*(x*x-1.),C*x*y-z*s,C*x*z+y*s) ;
      R.y += vector3d<double>(C*x*y+z*s,C*(y*y-1.),C*y*z-x*s) ;
      R.z += vector3d<double>(C*x*z-y*s,C*y*z+x*s,C*(z*z-1.)) ;
      s = -s ;
      Rinv.x += vector3d<double>(C*(x*x-1.),C*x*y-z*s,C*x*z+y*s) ;
      Rinv.y += vector3d<double>(C*x*y+z*s,C*(y*y-1.),C*y*z-x*s) ;
      Rinv.z += vector3d<double>(C*x*z-y*s,C*y*z+x*s,C*(z*z-1.)) ;
    }
    template<class realT> vector3d<realT> rotate_vec(vector3d<realT> v) const {
      return vector3d<realT>(R.x.x*v.x+R.x.y*v.y+R.x.z*v.z,
			     R.y.x*v.x+R.y.y*v.y+R.y.z*v.z,
			     R.z.x*v.x+R.z.y*v.y+R.z.z*v.z) ;
      // return dot(R,v) ;
    }
    template<class realT> tensor3d<realT> rotate_tensor(const tensor3d<realT> &t) const {
      tensor3d<realT> tRinv ;
      tRinv.x.x = t.x.x*Rinv.x.x+t.x.y*Rinv.y.x+t.x.z*Rinv.z.x ;
      tRinv.y.x = t.y.x*Rinv.x.x+t.y.y*Rinv.y.x+t.y.z*Rinv.z.x ;
      tRinv.z.x = t.z.x*Rinv.x.x+t.z.y*Rinv.y.x+t.z.z*Rinv.z.x ;

      tRinv.x.y = t.x.x*Rinv.x.y+t.x.y*Rinv.y.y+t.x.z*Rinv.z.y ;
      tRinv.y.y = t.y.x*Rinv.x.y+t.y.y*Rinv.y.y+t.y.z*Rinv.z.y ;
      tRinv.z.y = t.z.x*Rinv.x.y+t.z.y*Rinv.y.y+t.z.z*Rinv.z.y ;

      tRinv.x.z = t.x.x*Rinv.x.z+t.x.y*Rinv.y.z+t.x.z*Rinv.z.z ;
      tRinv.y.z = t.y.x*Rinv.x.z+t.y.y*Rinv.y.z+t.y.z*Rinv.z.z ;
      tRinv.z.z = t.z.x*Rinv.x.z+t.z.y*Rinv.y.z+t.z.z*Rinv.z.z ;

      tensor3d<realT> RtRinv ;
      RtRinv.x.x = R.x.x*tRinv.x.x+R.x.y*tRinv.y.x+R.x.z*tRinv.z.x ;
      RtRinv.y.x = R.y.x*tRinv.x.x+R.y.y*tRinv.y.x+R.y.z*tRinv.z.x ;
      RtRinv.z.x = R.z.x*tRinv.x.x+R.z.y*tRinv.y.x+R.z.z*tRinv.z.x ;

      RtRinv.x.y = R.x.x*tRinv.x.y+R.x.y*tRinv.y.y+R.x.z*tRinv.z.y ;
      RtRinv.y.y = R.y.x*tRinv.x.y+R.y.y*tRinv.y.y+R.y.z*tRinv.z.y ;
      RtRinv.z.y = R.z.x*tRinv.x.y+R.z.y*tRinv.y.y+R.z.z*tRinv.z.y ;

      RtRinv.x.z = R.x.x*tRinv.x.z+R.x.y*tRinv.y.z+R.x.z*tRinv.z.z ;
      RtRinv.y.z = R.y.x*tRinv.x.z+R.y.y*tRinv.y.z+R.y.z*tRinv.z.z ;
      RtRinv.z.z = R.z.x*tRinv.x.z+R.z.y*tRinv.y.z+R.z.z*tRinv.z.z ;
      return RtRinv ;
      //      return product(R,product(t,Rinv)) ;
    }
    template<class realT> vector3d<realT> transform(vector3d<realT> v) const {
      const vector3d<realT> vpt1(v.x+t1.x,v.y+t1.y,v.z+t1.z) ;
      const vector3d<realT> rvpt1=rotate_vec(vpt1) ;
      return vector3d<realT>(rvpt1.x+t2.x,rvpt1.y+t2.y,rvpt1.z+t2.z) ;
      //      return vector3d<realT>(rotate_vec(vpt1)+t2 ;
    }
  } ;

  struct periodic_info {
    std::string name ;
    bool master, processed ;
    vector3d<double> center, v, translate ;
    double angle ;
    entitySet bset ;
    Entity bc_num ;
    periodic_info() {
      name = "PERIODIC" ;
      master = false ;
      processed = false ;
      center = vector3d<double>(0.,0.,0.) ;
      v = vector3d<double>(1.,0.,0.) ;
      translate = center ;
      angle = 0 ;
      bset = EMPTY ;
      bc_num = 0 ;
    }
  } ;

  template<> struct data_schema_traits<Loci::rigid_transform> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(Loci::rigid_transform()) ;
      LOCI_INSERT_TYPE(ct,Loci::rigid_transform,t1) ;
      LOCI_INSERT_TYPE(ct,Loci::rigid_transform,t2) ;
      LOCI_INSERT_TYPE(ct,Loci::rigid_transform,R) ;
      LOCI_INSERT_TYPE(ct,Loci::rigid_transform,Rinv) ;
      return DatatypeP(ct) ;

    }
  } ;

  /// Reads an FVM grid and generates facts from the grid information.
  ///
  /// Creates facts for the `cl`, `cr`, `pos`, `face2node`, `boundary_names`,
  /// and `boundary_tags`. The `volumeTag(<tagName>)` fact is created for all
  /// volume tags found in the grid.
  /// @param[in,out] facts Fact database to add grid-related facts into.
  /// @param[in] filename Input grid filename.
  /// @param[in] cellwts Optional cell weights for partitioning.
  /// @return true if success
  bool readFVMGrid(fact_db &facts, std::string filename,storeRepP cellwts=0) ;

  /// Sets up FVM grid in the fact database.
  /// @param[in,out] facts Fact database.
  /// @param[in] filename Input grid filename.
  /// @param[in] cellwts Optional cell weights for partitioning.
  /// @return true if success
  bool setupFVMGrid(fact_db &facts, std::string filename,storeRepP cellwts=0) ;

  /// Sets up FVM grid in the fact database for the case where cell weights are
  /// stored in a storeRep.
  /// @param[in,out] facts Fact database to set up grid into.
  /// @param[in] filename Input grid filename.
  /// @param[in] cellwts  cell weights for partitioning.
  /// @return true if success
  bool setupFVMGridWithWeightInStore(fact_db &facts, std::string filename, storeRepP cellwt );

  /// Sets up FVM grid in the fact database for the case where cell weights are
  /// stored in a file.
  /// @param[in,out] facts Fact database to set up grid into.
  /// @param[in] filename Input grid filename.
  /// @param[in] weightfile  filename for file containing cell weights.
  /// @return true if success
  bool setupFVMGridWithWeightInFile(fact_db &facts, std::string filename, std::string weightfile);

  /// Gets the BC information from a VOG grid.
  /// @param[in] filename Name of the VOG grid file.
  /// @param[out] boundary_ids Vector of pairs of BC ID and name.
  /// @return true if success
  bool readBCfromVOG(std::string filename,
                     std::vector<std::pair<int,std::string> > &boundary_ids) ;

  /// Sets up facts related to grid boundary conditions.
  /// @param[in,out] facts Fact database to set up grid into.
  void setupBoundaryConditions(fact_db &facts) ;

  /// Creates the `lower`, `upper`, and `boundary_map` facts.
  /// @param[in,out] facts Fact database to set up grid into.
  void createLowerUpper(fact_db &facts) ;

  /// Creates `edges`, `edge2node` and `node2edge` facts.
  /// @param[in,out] facts Fact database to set up grid into.
  void createEdgesPar(fact_db& facts) ;


  extern void writeVOG(std::string filename,store<vector3d<double> > &pos,
                       Map &cl, Map &cr, multiMap &face2node,
                       std::vector<std::pair<int,std::string> > surfaceids) ;
  extern void writeVOG(std::string filename,store<vector3d<double> > &pos,
                       Map &cl, Map &cr, multiMap &face2node,
                       std::vector<std::pair<int,std::string> >& surfaceids,
                       std::vector<std::pair<std::string,entitySet> >& volTags) ;
  bool readGridVOG(std::vector<entitySet> &local_nodes,
                   std::vector<entitySet> &local_faces,
                   std::vector<entitySet> &local_cells,
                   store<vector3d<double> > &pos, Map &cl, Map &cr,
                   multiMap &face2node,
                   store<std::string> &boundary_names,
                   store<std::string> &boundary_tags,
                   std::vector<std::pair<std::string,entitySet> > &volTags,
                   int max_alloc, std::string filename) ;

  void setupOverset(fact_db &facts) ;

  template < class T> inline void setupPosAutoDiff(fact_db &facts, T val) {
#if defined(USE_AUTODIFF) || defined(MULTIFAD)
    store<vector3d<T> > pout ;
    {
      store<vector3d<double> > pin ;
      pin.setRep(facts.get_fact("pos")) ;

      entitySet dom = pin.domain() ;
      pout.allocate(dom) ;
      FORALL(dom,ii) {
	pout[ii] = vector3d<T>(pin[ii].x,
			       pin[ii].y,
			       pin[ii].z) ;
      } ENDFORALL ;
    }
    facts.replace_fact("pos",pout.Rep()) ;
#endif
  }

  template<class T> inline void setupPosAutoDiff(fact_db &facts,
						 std::string filename,
						 T val) {
#if defined(USE_AUTODIFF) || defined(MULTIFAD)
    setupPosAutoDiff(facts,val) ;
#endif
  }

  typedef enum {SIMPLE=0,SFC=1,ORB=2,GRAPH=3,RANDOM=4} partitionerSelector ;
  extern partitionerSelector partitionerMethod ;

}



#endif
