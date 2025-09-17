//#############################################################################
//#
//# Copyright 2015-2025, Mississippi State University
//#
//# This file is part of the flowPsi computational fluid dynamics solver.
//#
//# The flowPsi solver is free software: you can redistribute it and/or modify
//# it under the terms of the GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The flowPsi solver is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License
//# along with the flowPsi solver.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#ifndef FLOWPSIIO_H
#define FLOWPSIIO_H
#include <Loci>
#include <fstream>

namespace flowPsi {

  class integratedFileDBManager {
  public:
    // register routine to clean up errors on completion
    integratedFileDBManager() ;
    // clean up opened files on destruction
    ~integratedFileDBManager() ;
  } ;
  
  extern std::ofstream *getStreamFluxFile(const string filename, bool truncate) ;
  
  class variableOperatorList {
  public:
    struct listset {
      const char *name ;
      const char *nodalValue ;
      int val ;
      listset *next ;
    } ;
    static listset *varlist ;
    variableOperatorList(const char *name, int val) {
      listset *p = new listset ;
      p->name = name ;
      p->nodalValue = 0 ;
      p->val = val ;
      p->next = varlist ;
      varlist = p ;
    }
    variableOperatorList(const char *name, const char *nodalval, int val) {
      listset *p = new listset ;
      p->name = name ;
      p->nodalValue = nodalval ;
      p->val = val ;
      p->next = varlist ;
      varlist = p ;
    }
  } ;

  class scalar_node_output : public Loci::pointwise_rule {
    std::string var_name ;
    std::string value_name ;
    Loci::const_store<double> c2n ;
    Loci::const_param<std::string> plot_postfix ;
    Loci::const_param<std::string> modelName ;
    Loci::param<bool> OUTPUT ;
  public:
    scalar_node_output(const char *vname, const char *valname) ;
    virtual void compute(const Loci::sequence &seq) ;
  } ;

  class vector_node_output : public Loci::pointwise_rule {
    std::string var_name ;
    std::string value_name ;
    Loci::const_store<vector3d<double> > c2n ;
    Loci::const_param<std::string> plot_postfix ;
    Loci::const_param<std::string> modelName ;
    Loci::param<bool> OUTPUT ;
  public:
    vector_node_output(const char *vname, const char *valname) ;
    virtual void compute(const Loci::sequence &seq) ;
  } ;

  class dump_boundary_scalar : public Loci::pointwise_rule {
    std::string var_name ;
    std::string value_name ;
    Loci::const_store<real> var ;
    Loci::const_param<string> plot_postfix ;
    Loci::const_param<string> modelName ;
    Loci::param<bool> OUTPUT ;
  public:
    dump_boundary_scalar(const char *vname, const char *valname) ;
    virtual void compute(const Loci::sequence &seq) ;
  } ;

  class dump_boundary_vector : public Loci::pointwise_rule {
    std::string var_name ;
    std::string value_name ;
    Loci::const_store<vect3d> var ;
    Loci::const_param<string> plot_postfix ;
    Loci::const_param<string> modelName ;
    Loci::param<bool> OUTPUT ;
  public:
    dump_boundary_vector(const char *vname, const char *valname) ;
    virtual void compute(const Loci::sequence &seq) ;
  } ;



  // note, these are bits
  enum boundaryConditionErrorCode { OUTFLOW_BECAME_INFLOW=1,
                                    YPLUS_TOO_LARGE_IN_WALL_LAW=2,
                                    FIXED_MASS_FAILED_CONVERGENCE=4,
                                    ISENTROPIC_INFLOW_FAILED_CONVERGENCE=8} ;


  extern unsigned long BCErrorCode ;
}

#define OUTPUT_SCALAR(X,Y) class OUT_##Y : public scalar_node_output {\
  public:                                                             \
  OUT_##Y() : scalar_node_output(X,#Y){}                              \
  }; register_rule<OUT_##Y> register_OUT_##Y ;                        \
  variableOperatorList OUT_INFO_##Y(#Y,X,0)

#define OUTPUT_SCALAR_ALWAYS(X,Y) class OUT_##Y : public scalar_node_output {\
  public:                                                               \
  OUT_##Y() : scalar_node_output(X,#Y){}                                \
  }; register_rule<OUT_##Y> register_OUT_##Y ;                          \
  variableOperatorList OUT_INFO_##Y(#Y,X,1) 

#define OUTPUT_VECTOR(X,Y) class VECOUT_##Y : public vector_node_output {  \
  public:                                                                  \
  VECOUT_##Y() : vector_node_output(X,#Y){}                                \
  }; register_rule<VECOUT_##Y> register_VECOUT_##Y ;                       \
  variableOperatorList OUT_INFO_##Y(#Y,0)
#define OUTPUT_VECTOR_ALWAYS(X,Y) class VECOUT_##Y : public vector_node_output {  \
  public:                                                                  \
  VECOUT_##Y() : vector_node_output(X,#Y){}                                \
  }; register_rule<VECOUT_##Y> register_VECOUT_##Y ;                       \
  variableOperatorList OUT_INFO_##Y(#Y,1)

#define OUTPUT_BNDRY_SCALAR(X,Y,Z) class OUTB_##Y : public dump_boundary_scalar { \
  public:								\
  OUTB_##Y() : dump_boundary_scalar(X,#Y){				\
      constraint(Z);constraint(X);}					\
  }; register_rule<OUTB_##Y> register_OUTB_##Y ;			\
  variableOperatorList OUT_INFO_##Y(#Y,1)

#define OUTPUT_BNDRY_VECTOR(X,Y,Z) class OUTB_##Y : public dump_boundary_vector { \
  public:								\
  OUTB_##Y() : dump_boundary_vector(X,#Y){				\
      constraint(Z);}							\
  }; register_rule<OUTB_##Y> register_OUTB_##Y ;			\
  variableOperatorList OUT_INFO_##Y(#Y,1)




#endif
