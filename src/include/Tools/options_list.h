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
#ifndef OPTIONS_LIST_H
#define OPTIONS_LIST_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <Tools/unit_type.h>

#include <istream>
#include <ostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <Tools/basic_types.h>

namespace Loci {
  class options_list ;

  enum option_value_type {NOT_ASSIGNED,REAL,NAME,FUNCTION,LIST,STRING,BOOLEAN,
                          UNIT_VALUE,NAME_ASSIGN} ;

  class option_values {
  public:
    typedef std::vector<option_values> value_list_type ;
  private:
    option_value_type value_type ;
    value_list_type value_list ;
    std::string name ;
    double real_value ;
    double real_grad ;
    double real_grad2 ;
    int grad_size ;
    std::vector<double> gradN ;
    bool boolean_value ;
    UNIT_type units_value ;
    friend class options_list ;
  public:
    option_values() { value_type = NOT_ASSIGNED ; real_value = 0 ; real_grad = 0; boolean_value = false ; grad_size = MFAD_SIZE ; gradN.resize(grad_size) ; for(int i=0;i<grad_size;i++) gradN[i] = 0; }

    option_value_type type_of() const { return value_type ; }

    void get_value(bool &b) const { b = boolean_value; }
    void get_value(double &r) const { r = real_value ; }
    void get_value(MFADd &r) const { r = MFADd(real_value,&gradN[0],grad_size) ; }
    void get_value(FADd &r) const { r = FADd(real_value,real_grad) ; }
    void get_value(FAD2d &r) const { r = FAD2d(real_value,real_grad,real_grad2) ; }
    void get_value(value_list_type &l) const { l = value_list ; }
    void get_value(std::string &n) const { n = name ; }
    void get_value(UNIT_type &ut) const { ut = units_value ; }
    
    std::ostream & Print(std::ostream &s) const ;
    std::istream & Input(std::istream &s) ;
  } ;

  inline std::ostream & operator<<(std::ostream &s, const option_values &ov)
    { return ov.Print(s) ; }
  inline std::istream & operator>>(std::istream &s, option_values &ov)
    { return ov.Input(s) ; }    

  class options_list {
  public:
    typedef option_values::value_list_type arg_list ;
    typedef std::map<std::string,option_values> option_map ;
    typedef std::set<std::string> option_set ;
    typedef std::list<std::string> option_namelist ;

  private:
    option_set set_of_options ;
    bool restrict_set ;
    option_map options_db ;

  protected:

  public:
    options_list(const std::string &s) ;
    options_list() {restrict_set=false; } ;


    option_namelist getOptionNameList() const {
      option_namelist l ;
      for(option_map::const_iterator mi=options_db.begin();
          mi!=options_db.end();
          mi++) {
        l.push_back(mi->first) ;
      }
      return l ;
    }
    option_value_type getOptionValueType(const std::string &option) const ;
    option_values getOption(const std::string &option) const ;
    void getOption(const std::string &option, bool &value) const ;
    void getOption(const std::string &option, double &value) const ;
    void getOption(const std::string &option, MFADd & value) const ;
    void getOption(const std::string &option, FAD2d & value) const ;
    void getOption(const std::string &option, FADd & value) const {
      FAD2d v ;
      getOption(option,v) ;
      value = FADd(v.value,v.grad) ;
    }
    void getOption(const std::string &option, UNIT_type &uvalue) const ;
    void getOption(const std::string &option, std::string &name) const ;
    void getOption(const std::string &option, arg_list &value_list) const ;
    void getOption(const std::string &option, std::string &name,
                   arg_list &value_list) const ;
    void getOptionUnits(const std::string &option, const std::string &units,
                        double &value) const ;
    void getOptionUnits(const std::string &option, const std::string &units,
                        MFADd &value) const ;
    void getOptionUnits(const std::string &option, const std::string &units,
                        FAD2d &value) const ;
    void getOptionUnits(const std::string &option, const std::string &units,
                        FADd &value) const {
      FAD2d v ;
      getOptionUnits(option,units,v) ;
      value = v ;
    }
    void getOptionUnits(const std::string &option, const std::string &units,
			vector3d<double> &value, double scale=1.0) const ;
    void getOptionUnits(const std::string &option, const std::string &units,
			vector3d<MFADd> &value) const ;
    void getOptionUnits(const std::string &option, const std::string &units,
			vector3d<MFADd> &value,
			MFADd scale) const ;
    void getOptionUnits(const std::string &option, const std::string &units,
			vector3d<FAD2d> &value) const ;
    void getOptionUnits(const std::string &option, const std::string &units,
			vector3d<FAD2d> &value,
			FAD2d scale) const ;

    void getOptionUnits(const std::string &option, const std::string &units,
			vector3d<FADd> &value) const {
      vector3d<FAD2d> v ;
      getOptionUnits(option,units,v) ;
      value = vector3d<FADd>(FADd(v.x),FADd(v.y),FADd(v.z)) ;
    }

    void getOptionUnits(const std::string &option, const std::string &units,
			vector3d<FADd> &value,
			FADd scale) const {
      vector3d<FAD2d> v ;
      FAD2d scalec = FAD2d(scale.value,scale.grad,0.0)  ;
      getOptionUnits(option,units,v,scalec) ;
      value = vector3d<FADd>(v.x,v.y,v.z) ;
    }

    void setOption(const std::string &option, bool value) ;
    void setOption(const std::string &option, double value) ;
    void setOption(const std::string &option, UNIT_type uvalue) ;
    void setOption(const std::string &option, const std::string &name) ;
    void setOption(const std::string &option, const arg_list &value_list) ;
    void setOption(const std::string &option, const std::string &name,
                   const arg_list &value_list) ;


    bool checkOption(const std::string &option, const std::string &name) const ;
    bool optionExists(const std::string &option) const ;

    std::ostream & Print(std::ostream &s) const ;
    std::istream & Input(std::istream &s) ;
    void Input(const arg_list &l) ;
  } ;

  inline std::ostream & operator<<(std::ostream &s, const options_list &ol)
    { return ol.Print(s) ; }
  inline std::istream & operator>>(std::istream &s, options_list &ol)
    { return ol.Input(s) ; }
}
#endif
