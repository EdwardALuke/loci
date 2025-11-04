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
#include <Tools/options_list.h>
#include <Tools/parse.h>
#include <Tools/debug.h>
#include <Tools/except.h>
#include <Tools/tools.h>

#include <sstream>

using std::string ;
using std::ostream ;
using std::istream ;
using std::ostringstream ;
using std::cerr ;
using std::endl ;


using std::make_pair ;

namespace Loci {

  options_list::options_list(const string &s) {
    int sz = s.length() ;
    restrict_set = true ;
    if(sz == 0) {
      restrict_set = false ;
      return ;
    }
    string option ;
    for(int i=0;i<sz;++i) {
      if(s[i] == ':') {
        if(option.size() != 0)
          set_of_options.insert(option) ;
        option = "" ;
      } else 
        option += s[i] ;
    }
    if(option.size() != 0)
      set_of_options.insert(option) ;
  }

  option_value_type
    options_list::getOptionValueType(const string &option) const
    {
      option_map::const_iterator tmp ;
      if((tmp = options_db.find(option)) == options_db.end())
        return NOT_ASSIGNED ;
      return (*tmp).second.value_type ;
    }

  option_values
    options_list::getOption(const string &option) const
    {
      option_map::const_iterator tmp ;
      if((tmp = options_db.find(option)) == options_db.end())
        return option_values() ;
      return (*tmp).second ;
    }  
  
  void options_list::getOption(const string &option, bool &value) const {
    option_map::const_iterator tmp ;
    if((tmp = options_db.find(option)) == options_db.end()) {

        string s = "WARNING:attempt to retrieve BOOLEAN type option " ;
        s+= option ;
        s+= " failed." ;
        
        throw StringError(s) ;
        return ;
    }
    warn((*tmp).second.value_type != BOOLEAN) ;
    if((*tmp).second.value_type == BOOLEAN)
      value = (*tmp).second.boolean_value ;
  }
  
  void options_list::getOptionUnits(const string &option,
                                    const string &units,
                                    double &value) const {
    if(optionExists(option)) {
      if(getOptionValueType(option) == Loci::REAL) {
        getOption(option,value) ;
      } else if(getOptionValueType(option) == Loci::UNIT_VALUE) {
        Loci::UNIT_type Tu ;
        getOption(option,Tu) ;
        if(!Tu.is_compatible(units)) {
            ostringstream oss ;
          oss << "wrong type of unit for " << option <<": " << Tu << endl ;
          oss << "should have units compatible with " << units << endl ;
          throw StringError(oss.str()) ;
        } else {
          value = Tu.get_value_in(units) ;
        }
      } else {
          ostringstream oss ;
          oss << "incorrect type for "<< option << endl ;
          throw StringError(oss.str()) ;
      }
    } else {
      ostringstream oss ;
      oss << "options list cannot find option " << option << endl ;
      throw StringError(oss.str()) ;
    }

  }

  void options_list::getOptionUnits(const string &option,
                                    const string &units,
                                    FAD2d &value) const {
    if(optionExists(option)) {
      if(getOptionValueType(option) == Loci::REAL) {
        getOption(option,value) ;
      } else if(getOptionValueType(option) == Loci::UNIT_VALUE) {
        Loci::UNIT_type Tu ;
        getOption(option,Tu) ;
        if(!Tu.is_compatible(units)) {
            ostringstream oss ;
          oss << "wrong type of unit for " << option <<": " << Tu << endl ;
          oss << "should have units compatible with " << units << endl ;
          throw StringError(oss.str()) ;
        } else {
          value = Tu.get_value_inD(units) ;
        }
      } else {
          ostringstream oss ;
          oss << "incorrect type for "<< option << endl ;
          throw StringError(oss.str()) ;
      }
    } else {
      ostringstream oss ;
      oss << "options list cannot find option " << option << endl ;
      throw StringError(oss.str()) ;
    }

  }

  void options_list::getOptionUnits(const string &option,
                                    const string &units,
                                    VFAD &value, int batch) const {
    if(optionExists(option)) {
      if(getOptionValueType(option) == Loci::REAL) {
        getOption(option,value,batch) ;
      } else if(getOptionValueType(option) == Loci::UNIT_VALUE) {
        Loci::UNIT_type Tu ;
        getOption(option,Tu) ;
        if(!Tu.is_compatible(units)) {
            ostringstream oss ;
          oss << "wrong type of unit for " << option <<": " << Tu << endl ;
          oss << "should have units compatible with " << units << endl ;
          throw StringError(oss.str()) ;
        } else {
          value = Tu.get_value_inVF(units,batch) ;
        }
      } else {
          ostringstream oss ;
          oss << "incorrect type for "<< option << endl ;
          throw StringError(oss.str()) ;
      }
    } else {
      ostringstream oss ;
      oss << "options list cannot find option " << option << endl ;
      throw StringError(oss.str()) ;
    }

  }
  
  void options_list::getOptionUnits(const std::string &vname, 
				    const std::string &units,
				    vector3d<double> &vec,
				    double scale) const {
    Loci::option_value_type ovt= getOptionValueType(vname) ;
    if(ovt == Loci::REAL) {
      double v ;
      getOption(vname,v) ;
      vec = vector3d<double>(v*scale,0,0) ;
    } else if(getOptionValueType(vname) == Loci::UNIT_VALUE) {
      Loci::UNIT_type vu ;
      getOption(vname,vu) ;
      if(!vu.is_compatible(units)) {
	ostringstream oss ;
        oss << "wrong type of units for vector " << vname
	    << ": " << vu << std::endl ;
	throw StringError(oss.str()) ;
      } else {
        double v ;
        v = vu.get_value_in(units) ;
        vec = vector3d<double>(v,0,0) ;
      }
    } else if(ovt == Loci::LIST) {
      Loci::options_list::arg_list value_list ;
      getOption(vname,value_list) ;
      if(value_list.size() != 3) {
	ostringstream oss ;
        oss << "error on reading '" << vname
	    <<"': vector input must contain 3 terms"
	    << std::endl ;
	throw StringError(oss.str()) ;
      }
      for(int i=0;i<3;++i)
        if(value_list[i].type_of() != Loci::REAL &&
           value_list[i].type_of() != Loci::UNIT_VALUE) {
	  ostringstream oss ;
          oss << "improper vector specification for '"
	      << vname << std::endl ;
	  throw StringError(oss.str()) ;
        }
      double vecval[3] ;
      for(int i=0;i<3;++i) {
        if(value_list[i].type_of() == Loci::UNIT_VALUE) {
          Loci::UNIT_type vu ;
          value_list[i].get_value(vu) ;
          if(!vu.is_compatible(units)) {
	    ostringstream oss ;
            oss << "wrong type of units for vector " << vname
		<< ": " << vu << std::endl ;
	    throw StringError(oss.str()) ;
          }
          vecval[i] = vu.get_value_in(units) ;
        } else {
          value_list[i].get_value(vecval[i]) ;
          vecval[i] *= scale ;
        }
      }
      vec.x = vecval[0] ;
      vec.y = vecval[1] ;
      vec.z = vecval[2] ;
    } else if(ovt == Loci::FUNCTION) {
      string name ;
      Loci::options_list::arg_list value_list ;
      getOption(vname,name,value_list) ;
      if(name != "polar") {
	ostringstream oss ;
        oss << "don't know coordinate function '" << name
	    <<"', defaulting to polar" << std::endl ;
	throw StringError(oss.str()) ;
      }
      if(value_list.size() != 3) {
	ostringstream oss ;
        oss << "error on reading '"
	    << vname << "': vector input must contain 3 terms"
	    << std::endl ;
	throw StringError(oss.str()) ;
      }
      for(int i=0;i<3;++i)
        if(value_list[i].type_of() != Loci::REAL &&
           value_list[i].type_of() != Loci::UNIT_VALUE) {
	  ostringstream oss ;
          oss << "improper vector specification for '"
	      << vname << std::endl ;
	  throw StringError(oss.str()) ;
        }
      double r=1 ,theta=0 ,eta=0 ;
      double conv = M_PI/180.0 ;
      if(value_list[0].type_of() == Loci::UNIT_VALUE) {
        Loci::UNIT_type vu ;
        value_list[0].get_value(vu) ;
        if(!vu.is_compatible(units)) {
	  ostringstream oss ;
          oss << "wrong type of units for vector " << vname
	      << ": " << vu << std::endl ;
	  throw StringError(oss.str()) ;
        }
        r = vu.get_value_in(units) ;
      } else {
        value_list[0].get_value(r) ;
        r *= scale ;
      }
      if(value_list[1].type_of() == Loci::UNIT_VALUE) {
        Loci::UNIT_type vu ;
        value_list[1].get_value(vu) ;
        if(!vu.is_compatible("radians")) {
	  ostringstream oss ;
          oss << "wrong type of units for vector " << vname
                    << ": " << vu << std::endl ;
	  throw StringError(oss.str()) ;
        }
        theta = vu.get_value_in("radians") ;
      } else {
        value_list[1].get_value(theta) ;
        theta *= conv  ;
      }
      if(value_list[2].type_of() == Loci::UNIT_VALUE) {
        Loci::UNIT_type vu ;
        value_list[2].get_value(vu) ;
        if(!vu.is_compatible("radians")) {
	  ostringstream oss ;
          oss << "wrong type of units for vector " << vname
                    << ": " << vu << std::endl ;
	  throw StringError(oss.str()) ;
        }
        eta = vu.get_value_in("radians") ;
      } else {
        value_list[2].get_value(eta) ;
        eta *= conv  ;
      }
      
      vec.x = r*cos(theta)*cos(eta) ;
      vec.y = r*sin(theta)*cos(eta) ;
      vec.z = r*sin(eta) ;
    } else {
      ostringstream oss ;
      oss << "unable to get vector type!" << std::endl ;
      throw StringError(oss.str()) ;
    }
  }

    void options_list::getOptionUnits(const std::string &vname, 
				    const std::string &units,
				    vector3d<FAD2d> &vec, 
				    FAD2d scale) const {
    Loci::option_value_type ovt= getOptionValueType(vname) ;
    if(ovt == Loci::REAL) {
      FAD2d v ;
      getOption(vname,v) ;
      vec = vector3d<FAD2d>(v*scale,FAD2d(0.0,0.0,0.0),FAD2d(0.0,0.0,0.0)) ;
    } else if(getOptionValueType(vname) == Loci::UNIT_VALUE) {
      Loci::UNIT_type vu ;
      getOption(vname,vu) ;
      if(!vu.is_compatible(units)) {
	ostringstream oss ;
        oss << "wrong type of units for vector " << vname
	    << ": " << vu << std::endl ;
	throw StringError(oss.str()) ;
      } else {
        FAD2d v ;
        v = vu.get_value_inD(units) ;
        vec = vector3d<FAD2d>(v,0.0,0.0) ;
      }
    } else if(ovt == Loci::LIST) {
      Loci::options_list::arg_list value_list ;
      getOption(vname,value_list) ;
      if(value_list.size() != 3) {
	ostringstream oss ;
        oss << "error on reading '" << vname
	    <<"': vector input must contain 3 terms"
	    << std::endl ;
	throw StringError(oss.str()) ;
      }
      for(int i=0;i<3;++i)
        if(value_list[i].type_of() != Loci::REAL &&
           value_list[i].type_of() != Loci::UNIT_VALUE) {
	  ostringstream oss ;
          oss << "improper vector specification for '"
	      << vname << std::endl ;
	  throw StringError(oss.str()) ;
        }
      FAD2d vecval[3] ;
      for(int i=0;i<3;++i) {
        if(value_list[i].type_of() == Loci::UNIT_VALUE) {
          Loci::UNIT_type vu ;
          value_list[i].get_value(vu) ;
          if(!vu.is_compatible(units)) {
	    ostringstream oss ;
            oss << "wrong type of units for vector " << vname
		<< ": " << vu << std::endl ;
	    throw StringError(oss.str()) ;
          }
          vecval[i] = vu.get_value_inD(units) ;
        } else {
          value_list[i].get_value(vecval[i]) ;
          vecval[i] *= scale ;
        }
      }
      vec.x = vecval[0] ;
      vec.y = vecval[1] ;
      vec.z = vecval[2] ;
    } else if(ovt == Loci::FUNCTION) {
      string name ;
      Loci::options_list::arg_list value_list ;
      getOption(vname,name,value_list) ;
      if(name != "polar") {
	ostringstream oss ;
        oss << "don't know coordinate function '" << name
	    <<"', defaulting to polar" << std::endl ;
	throw StringError(oss.str()) ;
      }
      if(value_list.size() != 3) {
	ostringstream oss ;
        oss << "error on reading '"
	    << vname << "': vector input must contain 3 terms"
	    << std::endl ;
	throw StringError(oss.str()) ;
      }
      for(int i=0;i<3;++i)
        if(value_list[i].type_of() != Loci::REAL &&
           value_list[i].type_of() != Loci::UNIT_VALUE) {
	  ostringstream oss ;
          oss << "improper vector specification for '"
	      << vname << std::endl ;
	  throw StringError(oss.str()) ;
        }
      FAD2d r=1 ,theta=0 ,eta=0 ;
      double conv = M_PI/180.0 ;
      if(value_list[0].type_of() == Loci::UNIT_VALUE) {
        Loci::UNIT_type vu ;
        value_list[0].get_value(vu) ;
        if(!vu.is_compatible(units)) {
	  ostringstream oss ;
          oss << "wrong type of units for vector " << vname
	      << ": " << vu << std::endl ;
	  throw StringError(oss.str()) ;
        }
        r = vu.get_value_inD(units) ;
      } else {
        value_list[0].get_value(r) ;
        r *= scale ;
      }
      if(value_list[1].type_of() == Loci::UNIT_VALUE) {
        Loci::UNIT_type vu ;
        value_list[1].get_value(vu) ;
        if(!vu.is_compatible("radians")) {
	  ostringstream oss ;
          oss << "wrong type of units for vector " << vname
                    << ": " << vu << std::endl ;
	  throw StringError(oss.str()) ;
        }
        theta = vu.get_value_inD("radians") ;
      } else {
        value_list[1].get_value(theta) ;
        theta *= conv  ;
      }
      if(value_list[2].type_of() == Loci::UNIT_VALUE) {
        Loci::UNIT_type vu ;
        value_list[2].get_value(vu) ;
        if(!vu.is_compatible("radians")) {
	  ostringstream oss ;
          oss << "wrong type of units for vector " << vname
                    << ": " << vu << std::endl ;
	  throw StringError(oss.str()) ;
        }
        eta = vu.get_value_inD("radians") ;
      } else {
        value_list[2].get_value(eta) ;
        eta *= conv  ;
      }
      
      vec.x = r*cos(theta)*cos(eta) ;
      vec.y = r*sin(theta)*cos(eta) ;
      vec.z = r*sin(eta) ;
    } else {
      ostringstream oss ;
      oss << "unable to get vector type!" << std::endl ;
      throw StringError(oss.str()) ;
    }
  }
  
  void options_list::getOptionUnits(const std::string &vname, 
				    const std::string &units,
				    vector3d<VFAD> &vec, 
				    VFAD scale,
                                    int batch) const {
    Loci::option_value_type ovt= getOptionValueType(vname) ;
    if(ovt == Loci::REAL) {
      VFAD v ;
      getOption(vname,v,batch) ;
      vec = vector3d<VFAD>(v*scale,VFAD(0.0),VFAD(0.0)) ;
    } else if(getOptionValueType(vname) == Loci::UNIT_VALUE) {
      Loci::UNIT_type vu ;
      getOption(vname,vu) ;
      if(!vu.is_compatible(units)) {
	ostringstream oss ;
        oss << "wrong type of units for vector " << vname
	    << ": " << vu << std::endl ;
	throw StringError(oss.str()) ;
      } else {
        VFAD v ;
        v = vu.get_value_inVF(units,batch) ;
        vec = vector3d<VFAD>(v,VFAD(0.0),VFAD(0.0)) ;
      }
    } else if(ovt == Loci::LIST) {
      Loci::options_list::arg_list value_list ;
      getOption(vname,value_list) ;
      if(value_list.size() != 3) {
	ostringstream oss ;
        oss << "error on reading '" << vname
	    <<"': vector input must contain 3 terms"
	    << std::endl ;
	throw StringError(oss.str()) ;
      }
      for(int i=0;i<3;++i)
        if(value_list[i].type_of() != Loci::REAL &&
           value_list[i].type_of() != Loci::UNIT_VALUE) {
	  ostringstream oss ;
          oss << "improper vector specification for '"
	      << vname << std::endl ;
	  throw StringError(oss.str()) ;
        }
      VFAD vecval[3] ;
      for(int i=0;i<3;++i) {
        if(value_list[i].type_of() == Loci::UNIT_VALUE) {
          Loci::UNIT_type vu ;
          value_list[i].get_value(vu) ;
          if(!vu.is_compatible(units)) {
	    ostringstream oss ;
            oss << "wrong type of units for vector " << vname
		<< ": " << vu << std::endl ;
	    throw StringError(oss.str()) ;
          }
          vecval[i] = vu.get_value_inVF(units,batch) ;
        } else {
          value_list[i].get_value(vecval[i],batch) ;
          vecval[i] *= scale ;
        }
      }
      vec.x = vecval[0] ;
      vec.y = vecval[1] ;
      vec.z = vecval[2] ;
    } else if(ovt == Loci::FUNCTION) {
      string name ;
      Loci::options_list::arg_list value_list ;
      getOption(vname,name,value_list) ;
      if(name != "polar") {
	ostringstream oss ;
        oss << "don't know coordinate function '" << name
	    <<"', defaulting to polar" << std::endl ;
	throw StringError(oss.str()) ;
      }
      if(value_list.size() != 3) {
	ostringstream oss ;
        oss << "error on reading '"
	    << vname << "': vector input must contain 3 terms"
	    << std::endl ;
	throw StringError(oss.str()) ;
      }
      for(int i=0;i<3;++i)
        if(value_list[i].type_of() != Loci::REAL &&
           value_list[i].type_of() != Loci::UNIT_VALUE) {
	  ostringstream oss ;
          oss << "improper vector specification for '"
	      << vname << std::endl ;
	  throw StringError(oss.str()) ;
        }
      VFAD r(1.0) ,theta(0.0) ,eta(0.0) ;
      double conv = M_PI/180.0 ;
      if(value_list[0].type_of() == Loci::UNIT_VALUE) {
        Loci::UNIT_type vu ;
        value_list[0].get_value(vu) ;
        if(!vu.is_compatible(units)) {
	  ostringstream oss ;
          oss << "wrong type of units for vector " << vname
	      << ": " << vu << std::endl ;
	  throw StringError(oss.str()) ;
        }
        r = vu.get_value_inVF(units,batch) ;
      } else {
        value_list[0].get_value(r,batch) ;
        r *= scale ;
      }
      if(value_list[1].type_of() == Loci::UNIT_VALUE) {
        Loci::UNIT_type vu ;
        value_list[1].get_value(vu) ;
        if(!vu.is_compatible("radians")) {
	  ostringstream oss ;
          oss << "wrong type of units for vector " << vname
                    << ": " << vu << std::endl ;
	  throw StringError(oss.str()) ;
        }
        theta = vu.get_value_inVF("radians",batch) ;
      } else {
        value_list[1].get_value(theta,batch) ;
        theta *= conv  ;
      }
      if(value_list[2].type_of() == Loci::UNIT_VALUE) {
        Loci::UNIT_type vu ;
        value_list[2].get_value(vu) ;
        if(!vu.is_compatible("radians")) {
	  ostringstream oss ;
          oss << "wrong type of units for vector " << vname
                    << ": " << vu << std::endl ;
	  throw StringError(oss.str()) ;
        }
        eta = vu.get_value_inVF("radians",batch) ;
      } else {
        value_list[2].get_value(eta,batch) ;
        eta *= conv  ;
      }
      
      vec.x = r*cos(theta)*cos(eta) ;
      vec.y = r*sin(theta)*cos(eta) ;
      vec.z = r*sin(eta) ;
    } else {
      ostringstream oss ;
      oss << "unable to get vector type!" << std::endl ;
      throw StringError(oss.str()) ;
    }
  }
  
  void options_list::getOptionUnits(const std::string &vname, 
				    const std::string &units,
				    vector3d<FAD2d> &vec) const {
    return getOptionUnits(vname,units,vec,FAD2d(1.0,0.0,0.0)) ;
  }

  void options_list::getOptionUnits(const std::string &vname, 
				    const std::string &units,
				    vector3d<VFAD> &vec, int batch) const {
    return getOptionUnits(vname,units,vec,VFAD(1.0),batch) ;
  }

  void options_list::getOption(const string &option, double &value) const {
    option_map::const_iterator tmp ;
    if((tmp = options_db.find(option)) == options_db.end()) {
        ostringstream oss ;
        oss << "attempt to retrieve REAL type option " << option
            << " failed." << endl ;
        throw StringError(oss.str()) ;
      return ;
    }
    warn((*tmp).second.value_type != REAL) ;
    if((*tmp).second.value_type == REAL)
      value = (*tmp).second.real_value ;
  }

  void options_list::getOption(const string &option, FAD2d &value) const {
    option_map::const_iterator tmp ;
    if((tmp = options_db.find(option)) == options_db.end()) {
        ostringstream oss ;
        oss << "attempt to retrieve REAL type option " << option
            << " failed." << endl ;
        throw StringError(oss.str()) ;
      return ;
    }
    warn((*tmp).second.value_type != REAL) ;
    if((*tmp).second.value_type == REAL)
      value = FAD2d((*tmp).second.real_value,
		    (*tmp).second.real_grad, 
		    (*tmp).second.real_grad2) ;
  }

  void options_list::getOption(const string &option, VFAD &value,int batch) const {
    option_map::const_iterator tmp ;
    if((tmp = options_db.find(option)) == options_db.end()) {
        ostringstream oss ;
        oss << "attempt to retrieve REAL type option " << option
            << " failed." << endl ;
        throw StringError(oss.str()) ;
      return ;
    }
    warn((*tmp).second.value_type != REAL) ;
    if((*tmp).second.value_type == REAL) {
      value.data.value = (*tmp).second.real_value ;
      for(size_t i=0;i<VFAD::maxN;++i)
        value.data.grad[i] = 0 ;

      if(batch == 0) {
        value.data.grad[0] = (*tmp).second.real_grad  ;
        for(size_t i=1;i<min((*tmp).second.gradN.size()+1,VFAD::maxN);++i)
          value.data.grad[i] = (*tmp).second.gradN[i-1] ;
      } else {
        size_t offset = batch*VFAD::maxN ;
        size_t num = min(VFAD::maxN, offset-((*tmp).second.gradN.size()+1)) ;
        for(size_t i=0;i<num;++i)
          value.data.grad[i] = (*tmp).second.gradN[i+offset-1] ;
      }
    } else { 
      cerr << "(*tmp).second.value_type = " <<
        (*tmp).second.value_type <<endl  ;
    }

  }
  
  void options_list::getOption(const string &option, UNIT_type &uvalue) const {
    option_map::const_iterator tmp ;
    if((tmp = options_db.find(option)) == options_db.end()) {
        ostringstream oss ;
        oss << "attempt to retrieve UNIT_VALUE type option " << option
            << " failed." << endl ;
        throw StringError(oss.str()) ;
        return ;
    }
    warn((*tmp).second.value_type != UNIT_VALUE) ;
    if((*tmp).second.value_type == UNIT_VALUE)
      uvalue = (*tmp).second.units_value ;
  }

  void options_list::getOption(const string &option, string &name) const {
    option_map::const_iterator tmp ;
    if((tmp = options_db.find(option)) == options_db.end()) {
        ostringstream oss ;
        oss << "WARNING:attempt to retrieve NAME type option " << option
            << " failed." << endl ;
        throw StringError(oss.str()) ;
        return ;
    }
    warn((*tmp).second.value_type != NAME &&
         (*tmp).second.value_type != STRING) ;
    name = (*tmp).second.name ;
  }

  void options_list::getOption(const string &option, string &name,
                               arg_list &value_list) const {
    option_map::const_iterator tmp ;
    if((tmp = options_db.find(option)) == options_db.end()) {
        ostringstream oss ;
        oss << "attempt to retrieve FUNCTION type option " << option
            << " failed." << endl ;
        throw StringError(oss.str()) ;
        return ;
    }
    warn((*tmp).second.value_type != FUNCTION) ;
    
    name = (*tmp).second.name ;
    value_list = (*tmp).second.value_list ;
  }

  void options_list::getOption(const string &option, arg_list &value_list) const {
    option_map::const_iterator tmp ;
    if((tmp = options_db.find(option)) == options_db.end()) {
        ostringstream oss ;
        oss << "attempt to retrieve LIST type option " << option
            << " failed." << endl ;
        throw StringError(oss.str()) ;
        return ;
    }
    warn((*tmp).second.value_type != LIST) ;
    
    value_list = (*tmp).second.value_list ;
  }

  void options_list::setOption(const string &option, bool value) {
    option_map::iterator tmp ;

    if((tmp = options_db.find(option)) == options_db.end()) {
      option_values v ;
      tmp = options_db.insert(make_pair(option,v)).first ;
    }
    (*tmp).second.value_type = BOOLEAN ;
    (*tmp).second.boolean_value = value ;
  }
  
  void options_list::setOption(const string &option, double value) {
    option_map::iterator tmp ;
    
    if((tmp = options_db.find(option)) == options_db.end()) {
      option_values v ;
      tmp = options_db.insert(make_pair(option,v)).first ;
    }
    (*tmp).second.value_type = REAL ;
    (*tmp).second.real_value = value ;
  }

  void options_list::setOption(const string &option, UNIT_type uvalue) {
    option_map::iterator tmp ;
    
    if((tmp = options_db.find(option)) == options_db.end()) {
      option_values v ;
      tmp = options_db.insert(make_pair(option,v)).first ;
    }
    (*tmp).second.value_type = UNIT_VALUE ;
    (*tmp).second.units_value = uvalue ;
  }

  void options_list::setOption(const string &option, const string &name) {
    option_map::iterator tmp ;
    
    if((tmp = options_db.find(option)) == options_db.end()) {
      option_values v ;
      tmp = options_db.insert(make_pair(option,v)).first ;
    }
    (*tmp).second.value_type = NAME ;
    (*tmp).second.name = name ;
    arg_list l ;
    (*tmp).second.value_list = l ;
  }

  void options_list::setOption(const string &option, const arg_list &value_list) {
    option_map::iterator tmp ;
    
    if((tmp = options_db.find(option)) == options_db.end()) {
      option_values v ;
      tmp = options_db.insert(make_pair(option,v)).first ;
    }
    (*tmp).second.value_type = LIST ;
    (*tmp).second.value_list = value_list ;
  }

  void options_list::setOption(const string &option, const string &name,
                               const arg_list &value_list) {
    option_map::iterator tmp ;
    
    if((tmp = options_db.find(option)) == options_db.end()) {
      option_values v ;
      tmp = options_db.insert(make_pair(option,v)).first ;
    }
    (*tmp).second.value_type = FUNCTION ;
    (*tmp).second.name = name ;
    (*tmp).second.value_list = value_list ;
  }

  bool
    options_list::checkOption(const string &option, const string &name) const
    {
      option_map::const_iterator tmp ;
      if((tmp = options_db.find(option)) == options_db.end()) {
        return false ;
      }
      warn((*tmp).second.value_type != NAME) ;
      return name == (*tmp).second.name ;
    }

  bool options_list::optionExists(const string &option) const {
    option_map::const_iterator tmp ;
    if((tmp = options_db.find(option)) == options_db.end()) {
      return false ;
    }
    return (*tmp).second.value_type != NOT_ASSIGNED ;
  }

  ostream &option_values::Print(ostream &s) const {
    switch(value_type) {
    case BOOLEAN:
      s << (boolean_value?"$true":"$false") ;
      break ;
    case REAL:
      s << real_value ;
      if(gradN.size() > 0) {
        s << "^[" ;
        s << real_grad ;
        for(size_t i=0;i<gradN.size();++i) 
          s << ' ' << gradN[i] ;
        s << "]" ;
      } else {
        if(real_grad != 0 || real_grad2 != 0)
          s << "^" << real_grad ;
        if(real_grad2 != 0) {
          s << "^^" << real_grad2 ;
        }
      }
        
      break ;
    case UNIT_VALUE:
      s << units_value ;
      break ;
    case NAME:
      s << name ;
      break ;
    case NAME_ASSIGN:
      s << name << '=' ;
      {
        value_list_type::const_iterator i=value_list.begin() ;
        i->Print(s) ;
        ++i ;
        if(i!=value_list.end()) 
          throw StringError("confused setup in NAME_ASSIGN") ;
      }
      break ;
    case STRING:
      s << "\"" << name << "\"" ;
      break ;
    case FUNCTION:
      s << name << "(" ;
      for(value_list_type::const_iterator i=value_list.begin();;) {
        if(i==value_list.end())
          break ;
        s << *i ;
        ++i ;
        if(i==value_list.end())
          break ;
        s << "," ;
      }
      s << ")" ;
      break ;
    case LIST:
      s << "[" ;
      for(value_list_type::const_iterator i=value_list.begin();;) {
        if(i==value_list.end())
          break ;
        i->Print(s) ;
        ++i ;
        if(i==value_list.end())
          break ;
        s << "," ;
      }
      s << "]" ;
      break ;
    default:
      break ;
    }
          
    return s ;
  }

  istream &option_values::Input(istream &s) {
    value_list_type l ;
    value_list = l ;
    parse::kill_white_space(s) ;
    if(parse::is_real(s)) {
      real_value = parse::get_real(s) ;
      real_grad = 0 ;
      gradN = std::vector<double>() ;
      if(s.peek()=='^') {
	s.get() ;
        if(s.peek() == '[') { // multiple gradients
          s.get() ;
          size_t N = 0 ;
          while(s.peek() != ']') {
            if(parse::is_real(s)) {
              N++ ;
              double val = parse::get_real(s) ;
              if(N == 1)
                real_grad = val ;
              else
                gradN.push_back(val) ;
              
              parse::kill_white_space(s) ;
              if(s.peek() == ',') {
                s.get() ;
                parse::kill_white_space(s) ;
              }
            } else
              throw StringError("Expecting real in derivative list!") ;
          }
          s.get() ;
        } else {
          real_grad = parse::get_real(s) ;
        }
      }
      real_grad2 = 0 ;
      if(s.peek()=='^') {
	while(s.peek()=='^')
	  s.get() ;
	real_grad2 = parse::get_real(s) ;
      }

      parse::kill_white_space(s) ;
      if(parse::is_name(s)) {
        string units ;
        int ch = EOF ;
        int opens = 0 ;
        do {
	  int rc = 0 ;
          while(!s.eof() &&((ch=s.peek()) != EOF) &&
                (isalnum(ch) || ch=='*' || ch == '/' )) {
            units += s.get() ;
	    rc++ ;
	  }
          if(ch == '(') {
            units += s.get() ;
	    rc++ ;
            opens++ ;
          }
          if(opens != 0 && ch == ')') {
            units += s.get() ;
	    rc++ ;
            opens-- ;
          }
	  if(rc == 0) 
	    throw StringError("having trouble reading units in input") ;
        } while(opens!=0) ;

        units_value = UNIT_type(UNIT_type::MKS,"general",units,
                                real_value, real_grad, real_grad2,gradN) ;
        
        value_type = UNIT_VALUE ;
      } else
        value_type = REAL ;
    } else if(parse::is_name(s)) {
      name = parse::get_name(s) ;
      parse::kill_white_space(s) ;
      value_type = NAME ;
      if(s.peek() == '(') {
        value_type = FUNCTION ;
        s.get() ;
        parse::kill_white_space(s) ;
        while(!s.eof() && s.peek()!=')') {
          option_values ov ;
          s >> ov ;
          value_list.push_back(ov) ;
          parse::kill_white_space(s) ;
          if(s.peek() == ',') {
            s.get() ;
            parse::kill_white_space(s) ;
          }
        }
        if(s.peek() == ')')
          s.get() ;
      }
      if(s.peek() == '=') {
        s.get() ;
        parse::kill_white_space(s) ;
        value_type = NAME_ASSIGN ;
        option_values ov ;
        s >> ov ;
        value_list.push_back(ov) ;
        parse::kill_white_space(s) ;
      }
    } else if(s.peek() == '[') {
      s.get() ;
      value_type = LIST ;
      parse::kill_white_space(s) ;
      while(!s.eof() && s.peek()!=']') {
        option_values ov ;
        s >> ov ;
        parse::kill_white_space(s) ;
        if(s.peek() == '=') {
          s.get() ;
          parse::kill_white_space(s) ;
          
          if(ov.value_type != NAME && ov.value_type != STRING) {
              throw StringError("improper assignement in list") ;
              ov.value_type = NOT_ASSIGNED ;
          } else
            ov.value_type = NAME_ASSIGN ;
          option_values ov2 ;
          s >> ov2 ;
          if(ov2.value_type == LIST) {
            ov.value_list = ov2.value_list ;
          } else
            ov.value_list.push_back(ov2) ;
        }
        value_list.push_back(ov) ;
        if(s.peek() == ',') {
          s.get() ;
          parse::kill_white_space(s) ;
        }
      }
      if(s.peek() == ']')
        s.get() ;
    } else if(s.peek() == '"') {
      s.get() ;
      value_type = STRING ;
      char ch = s.get() ;
      while(ch != '"') {
        name += ch ;
        ch = s.get() ;
      }
    } else if(s.peek() == '$') {
      bool v = true ;
      if(s.peek() == '$') {
        s.get() ;
        string bvs = parse::get_name(s) ;
        if(bvs == "false")
          v = false ;
        else if(bvs != "true") {
            ostringstream oss ;
            oss << "option_list warning:" << endl ;
            oss << "boolean value can only be \"$true\"  or \"$false\""
                << endl;
            throw StringError(oss.str()) ;
        }
      }
      value_type = BOOLEAN ;
      boolean_value = v ;
    } else {
      throw StringError("error reading option_values") ;
    }
    return s ;
  }

  ostream &options_list::Print(ostream &s) const {
    s << "<" ;
    option_map::const_iterator tmp = options_db.begin() ;
    if(tmp != options_db.end()) {
      s<<(*tmp).first<<"="<<(*tmp).second ;
      for(++tmp;tmp!=options_db.end();++tmp)
        s<<","<<(*tmp).first<<"="<<(*tmp).second ;
    }
    s << ">" ;
    return s ;
  }

  istream &options_list::Input(istream &s) {
    try {
      std::set<string> option_parsed ;
      parse::kill_white_space(s) ;
      if(s.peek() != '<') {
          ostringstream oss ;
          oss << "format error in options_list::Input" << endl ;
          oss << "expected '<', got '" << s.peek() << "'" << endl ;
          throw StringError(oss.str()) ;
        return s ;
      }

      s.get() ;
    
      for(;;) {
        parse::kill_white_space(s) ;
        if(s.peek()=='>') {
          s.get() ;
          break ;
        }
        if(!parse::is_name(s)) {
            throw StringError("format error while reading option in option_list::Input") ;
          return s ;
        }
        string option = parse::get_name(s) ;

        if(option_parsed.find(option) == option_parsed.end()) {
          option_parsed.insert(option) ;
        } else {
          ostringstream oss ;
          oss << "input option '" << option << "' is reassigned value in option list" << endl ;
          throw StringError(oss.str()) ;
        }          
          
        if(set_of_options.find(option) == set_of_options.end()) {
          if(restrict_set) {
            ostringstream oss ;
            oss << "Invalid option name " << option << " in read options"
                << endl ;
            throw StringError(oss.str()) ;
          } else {
            set_of_options.insert(option) ;
          }
        } 

        try {
            parse::kill_white_space(s) ;
            option_values v ;

            if(s.peek() == '=') {
                s.get() ;
                parse::kill_white_space(s) ;
                
                s >> v ;
            } else {
                v.value_type = BOOLEAN ;
                v.boolean_value = true ;
            }

            option_map::iterator tmp ; 
            if((tmp = options_db.find(option)) == options_db.end()) {
                tmp = options_db.insert(make_pair(option,v)).first ;
            } else {
              (*tmp).second = v ;
            }
            parse::kill_white_space(s) ;
            if(s.peek()=='>') {
                s.get() ;
                break ;
            }
            if(s.peek()!=',') {
                ostringstream oss ;
                oss << "error reading option " << option << endl ;
                throw StringError(oss.str()) ;
            } else
              s.get() ;
        }
        catch(const BasicException &err) {
            err.Print(cerr) ;
            ostringstream oss ;
            oss << "error occurred while parsing option " << option
                << " in options list" ;
            throw StringError(oss.str()) ;
        }
      }
    }
    catch (const BasicException &err) {
        err.Print(cerr) ;
        throw StringError("options list parse error") ;
    }
    return s ;
  }

  void options_list::Input(const arg_list &l) {
    options_db.clear() ;

    std::set<string> option_parsed ;
    for(arg_list::const_iterator li=l.begin();li!=l.end();++li) {
      option_values v  ;
      if(li->type_of() == NAME_ASSIGN) {
        arg_list lv ;
        li->get_value(lv) ;
        warn(lv.size() != 1) ;
        v = *lv.begin() ;
      } else if(li->type_of() == NAME) {
        v.value_type = BOOLEAN ;
        v.boolean_value = true ;
      } else {
          ostringstream oss ;
          oss << "unable to parse input value " << *li << endl ;
          throw StringError(oss.str()) ;
          continue ;
      }
      string option  ;
      li->get_value(option) ;

      if(option_parsed.find(option) == option_parsed.end()) {
        option_parsed.insert(option) ;
      } else {
        ostringstream oss ;
        oss << "option '" << option << "' is reassigned value in option list" << endl ;
        throw StringError(oss.str()) ;
      }

      if(set_of_options.find(option) == set_of_options.end()) {
        if(restrict_set) {
          continue ;
        } else {
          set_of_options.insert(option) ;
        }
      } 

      option_map::iterator tmp ; 
      if((tmp = options_db.find(option)) == options_db.end()) {
        tmp = options_db.insert(make_pair(option,v)).first ;
      } else {
        (*tmp).second = v ;
      }
    }
  }
}
