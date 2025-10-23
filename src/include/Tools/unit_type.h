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
#ifndef UNIT_TYPE_H
#define UNIT_TYPE_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <Tools/except.h>

#include <Tools/parse.h>
#include <Tools/expr.h>
#include <Tools/autodiff.h>

#include <stack>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>

namespace Loci {

  class UNIT_type{

  public:
    std::string unit_kind;//the kind of unit ,such as pressure, time and so on
    std::string input_unit;// the unit which you input
    double input_value ;
    double gradient ;
    double secondGradient ;
    std::vector<double> gradientList ;

    enum unit_mode{MKS, CGS, check_available};
    // MKS: use MKS system
    // CGS: use CGS system
    //check_available: check the system if the unit available
    unit_mode mode;

    std::map<std::string,int> unit_num_map,unit_den_map;//containers of numerator and dnominator
    enum basic_unit_type {Length,Mass,Time,Temperature,Electric_current,
			  Amount_of_substance,Luminous_intensity, Angle, NoDim};
    double conversion_factor;

    //three tables of unit type - basic, composite, reference types----//
    struct basic_units{
      const char* name;
      basic_unit_type unit_type;
      double convert_factor;};  
    static basic_units basic_unit_table[] ;//MKS
    static basic_units cgs_basic_unit_table[];//CGS
   
    struct composite_units{
      const char* name;
      const char* derived_unit;
      double convert_factor;};  
    static composite_units composite_unit_table[] ;//MKS
    static composite_units cgs_composite_unit_table[] ;//CGS
  
    struct reference_units{
      const char* name;
      const char* refer_unit;
      double convert_factor;};  
    static reference_units reference_unit_table[] ;

    struct default_units{//If no unit input, check the default table 
      const char* default_type;//according to the unit_kind
      const char* default_name;
    };
    static default_units default_unit_table[] ;
  
  private:
    //check there is single temperature or temperature internal
    int is_single_temperature(const exprP input_expr);
    // if single temperature, do the conversion
    template<class T>
    void calculate_temperature(exprP &input_expr, T &val) {
    seperate_unit(unit_num_map,unit_den_map,input_expr); 
    get_conversion(unit_num_map,unit_den_map,conversion_factor);
    if(conversion_factor>0){
      conversion_factor=1;
      switch(is_single_temperature(input_expr)){
      case 1:
	val=(val+459.67)/1.8;
	break;
      case 2:
	val=val+273.15;
	break;
      case 3:
	val=val/1.8;
	break;
      }
    }else
      std::cerr<<"Not a unit"<<std::endl;
    }
    // when you need to convert to other temperature, reserve do it.
    template<class T>
    void reverse_calculate_temperature(exprP &input_expr,T &val) {
    seperate_unit(unit_num_map,unit_den_map,input_expr); 
    get_conversion(unit_num_map,unit_den_map,conversion_factor);
    if(conversion_factor>0){
      conversion_factor=1;
      switch(is_single_temperature(input_expr)){
      case 1:
	val=val*1.8-459.67;
	break;
      case 2:
	val=val-273.15;
	break;
      case 3:
	val=val*1.8;
	break;
      }
    }else
      std::cerr<<"Not a unit"<<std::endl;
    }

  public:
    exprP input(std::istream &in);// get the input unit
    void output(exprP &in_exp);// ouput the converted basic unit

    bool is_compatible(const std::string unit_str);//check two units compatible

    //get the value in converted unit
    double get_value_in(const std::string unit_str);
    FAD2d get_value_inD(const std::string unit_str);
    MFADd get_value_inM(const std::string unit_str);
    VFAD get_value_inVF(const std::string unit_str);
    void get_values_in(const std::string unit_str,double &val, double &grad, double &secondGrad, std::vector<double> &gradlist) ;
    UNIT_type(unit_mode in_mode, std::string in_kind, double in_value, std::string in_unit) {
      mode=in_mode,unit_kind=in_kind,input_value=in_value,input_unit=in_unit;
      input_value=in_value;
      gradient = 0 ;
      secondGradient = 0 ;
      exprP exp;
      exp=expression::create(input_unit);
      output(exp);
    }
    UNIT_type(unit_mode in_mode, std::string in_kind, std::string in_unit,
              double in_value, double in_gradient, double in_gradient2,
              const std::vector<double> &gradList) {
      mode=in_mode; unit_kind=in_kind; 
      input_value=in_value ;
      gradient = in_gradient ;
      secondGradient = in_gradient2 ;
      if(gradList.size() >0)
        gradient = gradList[0] ;
      for(size_t i=1;i<gradList.size();++i)
        gradientList.push_back(gradList[i]) ;
      input_unit=in_unit;
      exprP exp;
      exp=expression::create(input_unit);
      output(exp);
    }
    //#define UNIT_CONSTRUCTORS
#ifdef UNIT_CONSTRUCTORS
    UNIT_type(unit_mode in_mode, std::string in_kind, FADd in_value, std::string in_unit) {
      mode=in_mode; unit_kind=in_kind; 
      input_value=in_value.value ;
      gradient = in_value.grad ;
      secondGradient = 0 ;
      input_unit=in_unit;
      exprP exp;
      exp=expression::create(input_unit);
      output(exp);
    }
    UNIT_type(unit_mode in_mode, std::string in_kind, FAD2d in_value, std::string in_unit) {
      mode=in_mode,unit_kind=in_kind,input_unit=in_unit;
      input_value=in_value.value;
      gradient=in_value.grad ;
      secondGradient = in_value.grad2 ;
      exprP exp;
      exp=expression::create(input_unit);
      output(exp);
    }
    UNIT_type(unit_mode in_mode, std::string in_kind, MFADd in_value, std::string in_unit) {
      mode=in_mode,unit_kind=in_kind,input_unit=in_unit;
      gradient = in_value.grad[0] ;
      secondGradient = 0 ;
      for(int i=1;i<MFAD_SIZE;++i)
        gradientList.push_back(in_value.grad[i]) ;
      exprP exp;
      exp=expression::create(input_unit);
      output(exp);
    }
    UNIT_type(unit_mode in_mode, std::string in_kind, VFAD in_value, std::string in_unit) {
      mode=in_mode ;
      unit_kind=in_kind ;
      input_value = in_value.data.value ;
      gradient = in_value.data.grad[0] ;
      secondGradient = 0 ;
      for(int i=1;i<VFAD_SIZE;++i)
        gradientList.push_back(in_value.data.grad[i]) ;
      
      input_unit=in_unit ;
      exprP exp;
      exp=expression::create(input_unit);
      output(exp);
    }
#endif
    
    UNIT_type() {
      mode=MKS; unit_kind=""; 
      conversion_factor=1; input_value = 0;
      gradient = 0 ; secondGradient = 0 ;}

  private:
    bool is_reference_unit(std::string str);
    int where_reference_unit(std::string str);
    bool is_composite_unit(std::string str);
    int where_composite_unit(std::string str);
    bool is_basic_unit(std::string str);
    int where_basic_unit(std::string str);

    void build_lists(exprList &numerator, exprList &denominator, exprP input, bool isnum=true);
    void count_dim(std::map<std::string,int> &c_map, const exprList c_list);
    void rem_dup(std::map<std::string,int> &num_map,std::map<std::string,int> &den_map);
    std::map<std::string,int> combine_units(std::map<std::string,int> num_map,std::map<std::string,int> den_map);
    void seperate_unit(std::map<std::string,int> &num_map, std::map<std::string,int> &den_map,exprP p);

    void change_to_basic_unit(std::map<std::string,int>initial_map,std::map<std::string,int>&num_map,std::map<std::string,int>&den_map,double &conversion_factor);
    void get_conversion(std::map<std::string,int> &num_map, std::map<std::string,int> &den_map,double &conversion_factor);

    exprP set_default_unit(exprP &in_exp);
    int in_unit_kind();// check the unit_kind is available in db

  };

  inline std::ostream &operator<<(std::ostream &s, const UNIT_type &o_unit){
    s << o_unit.input_value  ;
    if(o_unit.gradientList.size() > 0) {
      s << "^[" << o_unit.gradient ;
      for(size_t i=0;i<o_unit.gradientList.size();++i)
        s << " " << o_unit.gradientList[i] ;
      s << "]" ;
    } else {
      if(o_unit.gradient != 0 || o_unit.secondGradient!=0) {
        s << "^" << o_unit.gradient ;
        if(o_unit.secondGradient != 0)
          s << "^^" <<o_unit.secondGradient ;
      }
    }
    if(o_unit.input_unit != "")
      s << ' ' << o_unit.input_unit ;
    return s ;
  }

  inline std::istream &operator>>(std::istream &s, UNIT_type &unit){
    exprP exp1;
    exp1=unit.input(s);//get the expression for unit
    unit.output(exp1); // get the unit exressed by basic units
    return s;
  }

  inline void unit_error(const int e_no, const std::string &err)
    {
        throw StringError(err) ;
    }
}
#endif
