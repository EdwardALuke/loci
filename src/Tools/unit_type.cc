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
#include <Tools/unit_type.h>
#include <Tools/tools.h>

#include <list>

using namespace std ;

namespace Loci {

  //three tables of unit type - basic, composite, reference types----//
  UNIT_type::basic_units UNIT_type::basic_unit_table[]={
    {"meter",Length,1},
    {"kilogram",Mass,1},
    {"second",Time,1},
    {"kelvin",Temperature,1},
    {"ampere",Electric_current,1},
    {"mole",Amount_of_substance,1},
    {"candela",Luminous_intensity,1},
    {"radians", Angle,1},
    {"none", NoDim,1},
    {0,NoDim,0}
  };
  
  UNIT_type::basic_units UNIT_type::cgs_basic_unit_table[]={
    {"centimeter",Length,1},
    {"gram",Mass,1},
    {"second",Time,1},
    {"kelvin",Temperature,1},
    {"ampere",Electric_current,1},
    {"mole",Amount_of_substance,1},
    {"candela",Luminous_intensity,1},
    {"radians", Angle,1},
    {"none", NoDim,1},
    {0,NoDim,0}
  };

  UNIT_type::composite_units UNIT_type::composite_unit_table[]={
    {"rad","radians",1},
    {"deg","radians",0.01745329251994329576923690768488612713443},
    {"degrees","radians",0.01745329251994329576923690768488612713443},
    {"rotations","radians",2*M_PI},
    {"one","meter/meter",1},
      
    //abbreviation of SI
    {"m","meter",1},
    {"kg","kilogram",1},
    {"g","kilogram",0.001},
    {"gram","kilogram",0.001},
    {"s","second",1},
    {"sec","second",1},
    {"K","kelvin",1},
    {"A","ampere",1},
    {"mol","mole",1},
    {"kmol","mole",1000},
    {"kmole","mole",1000},
    {"cd","candela",1},

    //metric system
    {"centimeter","meter",0.01},//length
    {"cm","meter",0.01},
    {"micrometer","meter",1e-6},
    {"micron","meter",1e-6},
    {"angstrom","meter",1e-10},
    {"kilometer","meter",1000},
    {"km","meter",1000},
    {"millimeter","meter",0.001},
    {"mm","meter",0.001},

    {"minute","second",60},//time
    {"min","second",60},
    {"hour","second",3600},
    {"h","second",3600},
    {"day","second",86400},
    {"d","second",86400},
    {"year","second",31536000},
    {"y","second",31536000},
    {"fortnight","second", 1209600}, // 14 days
    {"shake","second",1e-8},
    {"millisecond","second",1e-3},
    {"ms","second",1e-3},
    {"microsecond","second",1e-6},
    {"us","second",1e-6},
    {"usec","second",1e-6},
    {"nanosecond","second",1e-9},
    {"ns","second",1e-9},
    {"picosecond","second",1e-12},

    {"Fahrenheit","kelvin",0.5555556},//temperature interval
    {"F","kelvin",0.5555556},
    {"Rankine","kelvin",0.5555556},
    {"R","kelvin",0.5555556},
    {"Celsius","kelvin",1},
    {"C","kelvin",1},

    {"centiare","m*m",1},//area

    //composite unit
    {"hertz","second/second/second",1},
    {"newton","kilogram*meter/second/second",1},
    {"N","kilogram*meter/second/second",1},
    {"joule","kilogram*meter/second/second*meter",1},
    {"J","kilogram*meter/second/second*meter",1},
    {"watt","kilogram*meter/second/second*meter/second",1},
    {"W","kilogram*meter/second/second*meter/second",1},
    {"pascal","kilogram*meter/second/second/meter/meter",1},
    {"Pa","kilogram*meter/second/second/meter/meter",1},

    {"volt","kilogram*meter*meter/second/second/second/ampere",1},
    {"V","kilogram*meter*meter/second/second/second/ampere",1},
    {"tesla","kilogram/second/second/ampere",1},
    {"T","kilogram/second/second/ampere",1},

    {"siemens","second*second*second*ampere*ampere/kilogram/meter/meter",1},
    {"S","second*second*second*ampere*ampere/kilogram/meter/meter",1},
    {"ohm","kilogram*meter*meter/second/second/second/ampere/ampere",1},

    {"Hz",  "second/second/second",1},
    {"kHz", "second/second/second",1e3},
    {"MHz", "second/second/second",1e6},
    {"GHz", "second/second/second",1e9},

    {"rph","radians/second",1.7453292519943295769236907684887e-3},
    {"rpm","radians/second",0.10471975511965977461542144610932},
    {"rps","radians/second",6.2831853071795864769252867665592},

    {0,0,0}
  };

  UNIT_type::reference_units UNIT_type::reference_unit_table[]={
    {"acre","m*m",404.6873},//area
    {"are","m*m",100},
    {"a","m*m",100},
    {"barn","m*m",1e-28},
    {"b","m*m",1e-28},
    {"hectare","m*m",1.0e4},
    {"ha","m*m",1.0e4},

    {"calorie","joule",4.184},//thermochemical calorie
    {"cal","joule",4.184},
    {"electronvolt","joule",1.602177e-19},
    {"eV","joule",1.602177e-19},
    {"erg","joule",1.0e-7},
    {"kilocalorie","joule",4.184e3},
    {"kcal","joule",4.184e3},
    {"kW","watt",1e3},
    {"kWh","joule",3.6e6},
    {"kW*h","joule",3.6e6},
    {"quad","joule",1.055056e18},
    {"kiloton","joule",4.184e12},
     
    {"dyne","N",1.0e-5},//force
    {"dyn","N",1.0e-5},
    {"kgf","N",9.80665},
    {"kp","N",9.80665},
    {"kip","N",4.448222e3},
    {"ozf","N",2.780139e-1},
    {"poundal","N",1.38255e-1},
    {"lbf","N",4.448222},
    {"tonforce","N",8.896443e3},

    {"chain","m",2.011684e1},//length
    {"ch","m",2.011684e1},
    {"fathom","m",1.828804},
    {"fermi","m",1e-15},
    {"feet","m",0.3048},
    {"foot","m",0.3048},
    {"ft","m",0.3048},
    {"inch","m",2.54e-2},
    {"in","m",2.54e-2},
    {"microinch","m",2.54e-8},
    {"micro","m",1e-6},
    {"mil","m",2.54e-5},
    {"miles","m",1609.344},
    {"mile","m",1609.344},
    {"mi","m",1609.344},
    {"rod","m",5.029210},
    {"rd","m",5.029210},
    {"yard","m",0.9144},
    {"yd","m",0.9144},

    {"carat","kg",2.0e-4},//mass
    {"grain","kg",6.479891e-5},
    {"gr","kg",6.479891e-5},
    {"ounce","kg",2.834952e-2},
    {"oz","kg",2.834952e-2},
    {"pennyweight","kg",1.555174e-3},
    {"dwt","kg",1.555174e-3},
    {"pound","kg",0.4535924},
    {"lb","kg",0.4535924},
    {"lbm","kg",0.4535924},
    {"slug","kg",1.459390e1},
    {"tonne","kg",1e3},

    {"denier","kg/m",1.111111e-7},//mass divided by length
    {"tex","kg/m",1e-6},

    {"darcy","m*m",9.869233e-13},//permeability
    {"perm","kg/Pa/s/m/m",5.72135e-11},//(0C)

    {"horsepower","W",7.354988e2},//power(metric)
    {"hp","W",7.354988e2},

    {"atmosphere","Pa",1.01325e5},//pressure or stress
    {"bar","Pa",1e5},
    {"MPa","Pa",1e6}, // MegaPascals
    {"GPa","Pa",1e9}, // GigaPascals
    {"cmHg","Pa",1.333224e3},
    {"cmH2O","Pa",9.80665e1},
    {"ftHg","Pa",4.063666e4},
    {"ftH2O","Pa",2.989067e3},
    {"inHg","Pa",3.686389e3},
    {"inH2O","Pa",2.490889e2},
    {"ksi","Pa",6.894757e6},
    {"millibar","Pa",1e2},
    {"psi","Pa",6.894757e3},
    {"psia","Pa",6.894757e3},
    {"torr","Pa",1.333224e2},
    {"Torr","Pa",1.333224e2},

    {"centipoise","Pa*s",1.0e-3},//viscosity,dynamic
    {"cP","Pa*s",1e-3},
    {"poise","Pa*s",1e-1},
    {"P","Pa*s",1e-1},
    {"rhe","one/(Pa*s)",1e1},
    
    {"centistokes","m*m/s",1e-6},//viscosity,kinematic
    {"stokes","m*m/s",1e-4},
    {"St","m*m/s",1e-4},
    
    {"acrefoot","m*m*m",1.233489e3},//volume
    {"barrel","m*m*m",1.589873e-1},
    {"bll","m*m*m",1.589873e-1},
    {"bushel","m*m*m",3.523907e-2},
    {"cord","m*m*m",3.624556},
    {"cup","m*m*m",2.365882e-4},
    {"gallon","m*m*m",4.54609e-3},
    {"gal","m*m*m",4.54609e-3},
    {"gill","m*m*m",1.420653e-4},
    {"gi","m*m*m",1.420653e-4},
    {"liter","m*m*m",1e-3},
    {"peck","m*m*m",8.809768e-3},
    {"pk","m*m*m",8.809768e-3},
    {"pint","m*m*m",4.731765e-4},
    {"quart","m*m*m",9.463529e-4},
    {"stere","m*m*m",1},
    {"st","m*m*m",1},
    {"tablespoon","m*m*m",1.478676e-5},
    {"teaspoon","m*m*m",4.928922e-6},
    {"cc","m*m*m",1e-6}, // cubic centimeter

    {"clo","m*m*K/W",1.55e-1},//thermal insulance

    {"atm","pascal",101325},
    {"bar","pascal",100000},
    {"kPa","pascal",1000},
    {"Btu","joule",1055.87},

    {0,0,0}
  };

  UNIT_type::composite_units UNIT_type::cgs_composite_unit_table[]={
    {"one","centimeter/centimeter",1},
    //abbreviation of SI
    {"m","centimeter",0.01},
    {"kilogram","gram",1000},
    {"g","gram",1},
    {"kg","gram",1000},
    {"s","second",1},
    {"K","kelvin",1},
    {"A","ampere",1},
    {"mol","mole",1},
    {"cd","candela",1},

    //metric system
    {"kilometer","centimeter",100000},
    {"km","centimeter",100000},
    {"meter","centimeter",100},
    {"millimeter","centimeter",0.1},
    {"mm","centimeter",0.1},
    {"cm","centimeter",1},

    //time
    {"minute","second",60},
    {"min","second",60},
    {"hour","second",3600},
    {"h","second",3600},
    {"day","second",86400},
    {"d","second",86400},
    {"year","second",31536000},
    {"y","second",31536000},
    {"shake","second",1e-8},
    {"millisecond","second",1e-3},
    {"ms","second",1e-3},
    {"microsecond","second",1e-6},
    {"us","second",1e-6},
    {"usec","second",1e-6},
    {"nanosecond","second",1e-9},
    {"ns","second",1e-9},
    {"picosecond","second",1e-12},

    //temperature interval
    {"Fahrenheit","kelvin",0.5555556},
    {"F","kelvin",0.5555556},
    {"Rankine","kelvin",0.5555556},
    {"R","kelvin",0.5555556},
    {"Celsius","kelvin",1},
    {"C","kelvin",1},

    //composite unit
    {"newton","gram*centimeter/second/second",1e5},
    {"N","gram*centimeter/second/second",1e5},
    {"joule","gram*centimeter/second/second*centimeter",1e7},
    {"J","gram*centimeter/second/second*centimeter",1e7},
    {"watt","gram*centimeter/second/second*centimeter/second",1e7},
    {"W","gram*centimeter/second/second*centimeter/second",1e7},
    {"pascal","gram*centimeter/second/second/centimeter/centimeter",10},
    {"Pa","gram*centimeter/second/second/centimeter/centimeter",10},

    {"rph","radians/second",1.7453292519943295769236907684887e-3},
    {"rpm","radians/second",0.10471975511965977461542144610932},
    {"rps","radians/second",6.2831853071795864769252867665592},

    {0,0,0}

  };

  UNIT_type::default_units UNIT_type::default_unit_table[]={
    {"general"," "},
    {"length","meter"},
    {"Length","meter"},
    {"mass","kilogram"},
    {"Mass","kilogram"},
    {"time","second"},
    {"Time","second"},
    {"Temperature","K"},
    {"temperature","K"},
    {"electric_current","A"},
    {"Electric_current","A"},
    {"Amount_of_substance","mole"},
    {"amount_of_substance","mole"},
    {"Luminous_intensity","candela"},
    {"luminous_intensity","candela"},
    {"area","m*m"},
    {"density","kg/m/m/m"},
    {"energy","J"},
    {"power","W"},
    {"stress","Pa"},
    {"capacity","m*m*m"},
    {"density","kg/m/m/m"},
    {"flow","kg/s"},
    {"force","N"},
    {"Power","J"},
    {"power","J"},
    {"Velocity","m/s"},
    {"velocity","m/s"},
    {"speed","m/s"},
    {"Speed","m/s"},
    {"viscosityD","Pa*s"},
    {"viscosityK","m*m/s"},
    {"thermalConductivity","W/m/K"},
    {"volume","m*m*m"},
    {"Volume","m*m*m"},
    {"pressure","Pa"},
    {"Pressure","Pa"},
    {"heat","J"},
    {"Heat","J"},
    {0,0}
  };

  bool UNIT_type::is_reference_unit(string str){
    for(int i=0;reference_unit_table[i].name!=0;++i){
      if(reference_unit_table[i].name==str)
	return true;
    }
    return false;
  }

  int UNIT_type::where_reference_unit(string str){
    for(int i=0;reference_unit_table[i].name != 0;++i){
      if(reference_unit_table[i].name==str)
	return i;
    }
    return -1;
  }

  bool UNIT_type::is_composite_unit(string str){
    if(mode==MKS){
      for(int i=0;composite_unit_table[i].name!=0;++i){
	if(composite_unit_table[i].name==str)
	  return true;
      }
    }
    else if(mode==CGS){
      for(int i=0;cgs_composite_unit_table[i].name!=0;++i){
	if(cgs_composite_unit_table[i].name==str)
	  return true;
      }
    }
    else
      cout<<"No this mode "<<mode<<endl;
    return false;
  }

  int UNIT_type::where_composite_unit(string str){
    if(mode==MKS){
      for(int i=0;composite_unit_table[i].name!=0;++i){
	if(composite_unit_table[i].name==str)
	  return i;
      }
    }
    else if(mode==CGS){
      for(int i=0;cgs_composite_unit_table[i].name != 0;++i){
	if(cgs_composite_unit_table[i].name==str)
	  return i;
      }
    }
    else
      cout<<"No this mode "<<mode<<endl;
    return -1;
  }

  bool UNIT_type::is_basic_unit(string str){
    if(mode==MKS){
      for(int i=0;basic_unit_table[i].name != 0;++i){
	if(basic_unit_table[i].name==str)
	  return true;
      }
    }
    else if(mode==CGS){
      for(int i=0;cgs_basic_unit_table[i].name != 0;++i){
	if(cgs_basic_unit_table[i].name==str)
	  return true;
      }
    }
    else
      cout<<"No this mode "<<mode<<endl;
    return false;
  }

  int UNIT_type::where_basic_unit(string str){
    if(mode==MKS){
      for(int i=0;basic_unit_table[i].name!=0;++i){
	if(basic_unit_table[i].name==str)
	  return i;
      }
    }
    else if(mode==CGS){
      for(int i=0;cgs_basic_unit_table[i].name!=0;++i){
	if(cgs_basic_unit_table[i].name==str)
	  return i;
      }
    }
    else
      cout<<"No this mode "<<mode<<endl;
    return -1;
  }

  //------build lists of numerator and denominator-----------//
  void UNIT_type::build_lists(exprList &numerator, exprList &denominator, exprP input, bool isnum)
    {
      exprList::const_iterator li,lj;
      switch(input->op){ 
      case OP_NAME:
	if(isnum)
	  numerator.push_back(input);
	else
	  denominator.push_back(input);
	break;
      case OP_TIMES:
	for(li= input->expr_list.begin();li!=input->expr_list.end();++li)
	  build_lists(numerator,denominator,(*li),isnum);
	break;
      case OP_EXOR:
	for(li= input->expr_list.begin();li!=input->expr_list.end();++li){
	  if((*li)->op==OP_INT){
	    //	    cout<<endl<<(*li)->int_val<<endl;
	    lj= input->expr_list.begin();
	    for(int i=1; i!=(*li)->int_val;i++){
	      //cout<<endl<<"i  "<<i<<"  "<<(*lj)->name<<endl;
	      build_lists(numerator,denominator,(*lj),isnum);
	    }
	  }else if((*li)->op==OP_NAME){
	    //cout<<endl<<(*li)->name<<endl;
	    build_lists(numerator,denominator,(*li),isnum);
	  }else
	    cout<<"Wrong exponent!"<<endl;
	}
	break;
      case OP_DIVIDE:
        {
          bool new_isnum=!isnum;
          for(li= input->expr_list.begin();li!=input->expr_list.end();++li){
            build_lists(numerator,denominator,(*li),isnum);
            isnum=new_isnum;
          }
        }
	break;
      default:
          string s = "unknown operator type in unit_type expression parsing!" ;
          cerr << "op = " << input->op << endl ;
          unit_error(4,s) ;
      }
      //cout<<"size of numerator "<<numerator.size()<<endl;
    }

  //count the dimension of numerator or denominator
  void UNIT_type::count_dim(map<string,int> &c_map, const exprList c_list)
    {
      exprList::const_iterator li;
      /*cout<<"content of c_map"<<endl;
      map<string,int>::const_iterator mj=c_map.begin();
      for(;mj!=c_map.end();++mj)
      cout<<mj->first<<"  "<<mj->second<<endl;*/

      for(li= c_list.begin();li!=c_list.end();++li)
	{
	  int dim=1;
	  //cout<<"content of c_list"<<(*li)->name<<endl;
	  if(c_map.find((*li)->name)!=c_map.end())
	    c_map[(*li)->name]+=dim;
	  else
	    c_map[(*li)->name]=1;
	}
      /*map<string,int>::const_iterator mi=c_map.begin();
      for(;mi!=c_map.end();++mi)
      cout<<mi->first<<"  "<<mi->second<<endl;*/
    }

  //remove the duplicated units 
  void UNIT_type::rem_dup(map<string,int> &num_map,map<string,int> &den_map)
    {
      map<string,int>::iterator mi,mj;
      list<map<string,int>::iterator> delnum,delden ;
      for(mi= num_map.begin();mi!=num_map.end();++mi)
	{
	  for(mj= den_map.begin();mj!=den_map.end();++mj)
	    {
	      if((*mj).first==(*mi).first)
		{
		  int tmp=(*mi).second-(*mj).second;
		  if(tmp>0){
		    num_map[(*mi).first]=tmp;
                    delden.push_back(mj) ;
		  }
		  else if(tmp<0){
		    den_map[(*mi).first]=-tmp;
                    delnum.push_back(mi) ;
		  }
		  else{
                    delnum.push_back(mi) ;
                    delden.push_back(mj) ;
		  }
		}
	    }
	}
      while(!delnum.empty()) {
        num_map.erase(delnum.back()) ;
        delnum.pop_back() ;
      }
      while(!delden.empty()) {
        den_map.erase(delden.back()) ;
        delden.pop_back() ;
      }
      /*cout<<"Numerator: "<<endl;
      for(mi= num_map.begin();mi!=num_map.end();++mi)
	cout<<mi->first<<"  "<<mi->second<<endl;
      cout<<"Denominator: "<<endl;
      for(mj= den_map.begin();mj!=den_map.end();++mj)
      cout<<mj->first<<"  "<<mj->second<<endl;*/
    }

  //combine the elements of units//
  map<string,int> UNIT_type::combine_units(map<string,int> num_map,map<string,int> den_map)
    {
      map<string,int>::iterator mi,mj;
      if(num_map.size()!=0){
	for(mi= num_map.begin();mi!=num_map.end();++mi)
	  {
	    for(mj= den_map.begin();mj!=den_map.end();++mj)
	      {
		if((*mj).first==(*mi).first){
		  int tmp=(*mi).second+(*mj).second;
		  if(tmp>0)
		    num_map[(*mi).first]=tmp;
		} else
		  num_map.insert(*mj);
	      }
	  }
      }else{//if there is no element in num_map, insert the elements of den_map
	for(mj= den_map.begin();mj!=den_map.end();++mj){
	  num_map.insert(*mj);
	}
      }
      return num_map;
    }

  //seperate the input units into maps of numerator and denominator//
  void UNIT_type::seperate_unit(map<string,int> &num_map, map<string,int> &den_map,exprP p)    
    {
      exprList numerator,denominator;
      build_lists(numerator,denominator,p);
      //cout<<"Numerator:  ";
      //printing(cout,numerator);
      count_dim(num_map,numerator);
      // cout<<"Denominator:  ";
      //printing(cout,denominator);
      count_dim(den_map,denominator);
      //cout<<endl;
      //cout<<"Remove the duplicated units"<<endl;
      rem_dup(num_map,den_map);
    }

  //change map of numerator or denominator to basic type
  void UNIT_type::change_to_basic_unit(map<string,int>initial_map,map<string,int>&num_map,map<string,int>&den_map,double &conv_factor)
    {
      map<string,int>::iterator mi;
      exprP exp2 = 0 ;
      exprList numerator, denominator;
      conv_factor=1;
      //for MKS system
      if(mode==MKS){
	for(mi= initial_map.begin();mi!=initial_map.end();++mi) {
	  if(is_reference_unit((*mi).first)){ 
	    for(int i=0;i!=(*mi).second;i++){//for several same units
	      //cout<<"IS reference unit"<<endl;       
	      int where=where_reference_unit((*mi).first);
	      conv_factor=conv_factor*reference_unit_table[where].convert_factor;
	      exp2 = expression::create(reference_unit_table[where].refer_unit) ;
	      map<string,int> ref_num_map, ref_den_map ;
	      seperate_unit(ref_num_map,ref_den_map,exp2) ;
	      double ref_conversion = 1.0 ;
	      get_conversion(ref_num_map,ref_den_map,ref_conversion) ;
	      conv_factor *= ref_conversion ;
	      for(map<string,int>::const_iterator ri=ref_num_map.begin();
		  ri!=ref_num_map.end();++ri) {
		if(num_map.find((*ri).first)!=num_map.end())
		  num_map[(*ri).first] += (*ri).second ;
		else
		  num_map[(*ri).first] = (*ri).second ;
	      }
	      for(map<string,int>::const_iterator ri=ref_den_map.begin();
		  ri!=ref_den_map.end();++ri) {
		if(den_map.find((*ri).first)!=den_map.end())
		  den_map[(*ri).first] += (*ri).second ;
		else
		  den_map[(*ri).first] = (*ri).second ;
	      }
	      
	      //cout<<"conversion factor:"<<conv_factor<<endl;
	    }
	  } else if(is_composite_unit((*mi).first)){
	    
	    for(int i=0;i!=(*mi).second;i++){
	      //cout<<"IS composite unit"<<endl;
	      //cout<<(*mi).first<<(*mi).second<<endl;
	      int where=where_composite_unit((*mi).first);
	      conv_factor=conv_factor*composite_unit_table[where].convert_factor;
	      exp2=expression::create(composite_unit_table[where].derived_unit);
	      seperate_unit(num_map,den_map,exp2);
	      //cout<<"conversion factor:"<<conv_factor<<endl;
	    }
            
	  } else if(is_basic_unit((*mi).first)){
	    //cout<<"IS basic unit"<<endl;
	    if(num_map.find((*mi).first)!=num_map.end())
	      num_map[(*mi).first]=num_map[(*mi).first]+(*mi).second;
	    else
	      num_map[(*mi).first]=(*mi).second;
	    //cout<<(*mi).first<<(*mi).second<<endl;
	    //conv_factor=conv_factor*1;
	  } else{
              string tmp = (*mi).first + ": Not in MKS database" ;
              unit_error(4,tmp);
	  }
      }
      }
      // for CGS system
      else if(mode==CGS){
	for(mi= initial_map.begin();mi!=initial_map.end();++mi)
	  {
	  if(is_reference_unit((*mi).first)){ 
	    for(int i=0;i!=(*mi).second;i++){//for several same units
	      //cout<<"IS reference unit"<<endl;       
	      int where=where_reference_unit((*mi).first);
	      conv_factor=conv_factor*reference_unit_table[where].convert_factor;
	      exp2 = expression::create(reference_unit_table[where].refer_unit) ;
	      map<string,int> ref_num_map, ref_den_map ;
	      seperate_unit(ref_num_map,ref_den_map,exp2) ;
	      double ref_conversion = 1.0 ;
	      get_conversion(ref_num_map,ref_den_map,ref_conversion) ;
	      conv_factor *= ref_conversion ;
	      for(map<string,int>::const_iterator ri=ref_num_map.begin();
		  ri!=ref_num_map.end();++ri) {
		if(num_map.find((*ri).first)!=num_map.end())
		  num_map[(*ri).first] += (*ri).second ;
		else
		  num_map[(*ri).first] = (*ri).second ;
	      }
	      for(map<string,int>::const_iterator ri=ref_den_map.begin();
		  ri!=ref_den_map.end();++ri) {
		if(den_map.find((*ri).first)!=den_map.end())
		  den_map[(*ri).first] += (*ri).second ;
		else
		  den_map[(*ri).first] = (*ri).second ;
	      }
	      
	      //cout<<"conversion factor:"<<conv_factor<<endl;
	    }
	  }
	  else if(is_composite_unit((*mi).first)){
	    
	    for(int i=0;i!=(*mi).second;i++){
	      //cout<<"IS composite unit"<<endl;
	      //cout<<(*mi).first<<(*mi).second<<endl;
	      int where=where_composite_unit((*mi).first);
	      conv_factor=conv_factor*cgs_composite_unit_table[where].convert_factor;
	      exp2=expression::create(cgs_composite_unit_table[where].derived_unit);
	      seperate_unit(num_map,den_map,exp2);
	      //cout<<"conversion factor:"<<conv_factor<<endl;
	    }

	  }
	  else if(is_basic_unit((*mi).first)){
	    //cout<<"IS basic unit"<<endl;
	    if(num_map.find((*mi).first)!=num_map.end())
	      num_map[(*mi).first]=num_map[(*mi).first]+(*mi).second;
	    else
	      num_map[(*mi).first]=(*mi).second;
	    //cout<<(*mi).first<<(*mi).second<<endl;
	    conv_factor=conv_factor*1;
	  }
	  else{
              string tmp = (*mi).first + ": Not in CGS database!" ;
	    unit_error(5,tmp);
	    conv_factor=0;
	  }
	}
	}
    }

  //get the units convert to baisc types//
  void UNIT_type::get_conversion(map<string,int> &num_map, map<string,int> &den_map,double &conv_factor)
    {
      map<string,int> num_map2,den_map2,num_map3,den_map3;
      double num_conversion_factor,den_conversion_factor;
      int num_size,den_size;
  
      num_size=num_map.size();
      den_size=den_map.size();
      change_to_basic_unit(num_map,num_map2,den_map2,num_conversion_factor);
      //cout<<endl;
      change_to_basic_unit(den_map,num_map3,den_map3,den_conversion_factor);
  
      if(num_size>0&&den_size>0)
	conv_factor=num_conversion_factor/den_conversion_factor;
      else if(num_size>0&&den_size<=0)
	conv_factor=num_conversion_factor;
      else if(den_size>0&&num_size<=0)
	conv_factor=1/den_conversion_factor;
  
      //cout<<endl;
      num_map=combine_units(num_map2,den_map3); 
      den_map=combine_units(den_map2,num_map3); 
    }

  //------unit check---------//

  int UNIT_type::in_unit_kind(){
    for(int i=0;default_unit_table[i].default_type!=0;i++){
      if(unit_kind==default_unit_table[i].default_type)
	return i;
    }
    return -1;
  }

  //---if no unit input, set the default unit-------//
  exprP UNIT_type::set_default_unit(exprP &in_exp){
    int tmp=in_unit_kind();
    if(tmp!=-1)
      in_exp=expression::create(default_unit_table[tmp].default_name);
    else{
        unit_error(6,"Unit is not compatible with the kind of unit! ");
    }
    if(unit_kind=="general"){
      unit_error(7,"You must input a unit!");
    }
    return in_exp;
  }

  //------input the expression-------//
  exprP UNIT_type::input(istream &in){
    unit_kind="" ;
    //    input_unit="" ; existing input_unit will be default if no unit is
    // specified.
    mode = MKS ;
    input_value = 0 ;
    gradient = 0 ;
    secondGradient = 0 ;
    gradientList = std::vector<double>() ;

    while(in.peek() == ' ' || in.peek() == '\t')
      in.get() ; ;
    if(parse::is_real(in)) {
      input_value = parse::get_real(in) ;

      while(in.peek() == ' ' || in.peek() == '\t')
        in.get() ; ;
      if(in.peek() == '^') {
        // get gradient information
        in.get() ;
        while(in.peek() == ' ' || in.peek() == '\t')
          in.get() ;
        if(in.peek() == '[') {
          // Reading in vector of values
          in.get();
          while(in.peek() == ' ' || in.peek() == '\t')
            in.get() ;
          if(parse::is_real(in)) {
            gradient = parse::get_real(in) ;
            while(in.peek() == ' ' || in.peek() == '\t')
              in.get() ;
          }

          while(parse::is_real(in)) {
            gradientList.push_back(parse::get_real(in)) ;
            while(in.peek() == ' ' || in.peek() == '\t')
              in.get() ;
          }
          if(in.peek() != ']') {
            throw exprError("Syntax",
                            "Unit Type: Error parsing derivative block",
                            ERR_SYNTAX) ;
          }  else
            in.get() ;
        } else {
          if(parse::is_real(in)) {
            gradient = parse::get_real(in) ;
            while(in.peek() == ' ' || in.peek() == '\t')
              in.get() ;
            if(in.peek() == '^') {
              in.get() ;
              if(in.peek() == '^')
                in.get() ;
              if(parse::is_real(in))
                secondGradient = parse::get_real(in) ;
              else
                throw exprError("Syntax",
                                "Unit Type: Error parsing second derivative",
                                ERR_SYNTAX) ;
            }
          } else {
            throw exprError("Syntax",
                            "Unit Type: Error parsing gradient information",
                            ERR_SYNTAX) ;
          }
        }
      }
    } else {
      throw Loci::exprError("Syntax",
                            "Unit Value Type: value not parsed!",
                            ERR_SYNTAX) ;
    }
    while(in.peek() == ' ' || in.peek() == '\t')
      in.get() ;

    // parse unit information
    if(parse::is_name(in) || parse::is_token(in,"(")) {
      string parsed ;
      while(true) {
        parse::kill_white_space(in) ;
        if(parse::is_name(in)) {
          parsed += parse::get_name(in) ;
        } else if(parse::is_token(in,"(")) {
          parse::get_token(in,"(") ;
          parsed += '(' ;
          int cnt = 1 ;
          while(cnt > 0) {
            int ch = in.get() ;
            if(ch == EOF || in.eof()) {
              throw Loci::exprError("Syntax",
                                    "Unit Type: Error parsing unit expression",
                                    ERR_SYNTAX) ;
            }
            if(ch == '(')
              cnt++ ;
            if(ch == ')')
              cnt-- ;
            parsed += ch ;
          }
        } else {
          break ;
        }

        parse::kill_white_space(in) ;
        if(parse::is_token(in,"*")) {
          parse::get_token(in,"*") ;
          parsed += "*" ;
          continue ;
        }
        if(parse::is_token(in,"/")) {
          parse::get_token(in,"/") ;
          parsed += "/" ;
          continue ;
        }
        break ;
      }

      if(parsed.empty() || parsed[parsed.size()-1] == '*' || parsed[parsed.size()-1] == '/') {
        throw Loci::exprError("Syntax",
                              "Unit Type: Error parsing unit expression",
                              ERR_SYNTAX) ;
      }
      input_unit = parsed ;
    }
    exprP exp=expression::create(input_unit);
    return exp;
  }

  // check if the unit is single temperature //
  int UNIT_type::is_single_temperature(const exprP input_expr){
    if(input_expr->name=="Fahrenheit"||input_expr->name=="F")
      return 1;
    else if(input_expr->name=="Celsius"||input_expr->name=="C")
      return 2;
    else if(input_expr->name=="Rankine"||input_expr->name=="R")
      return 3;
    return 0;
  }


  //--- input unit expression and output basic units ---//
  void UNIT_type::output(exprP &in_exp){
    unit_den_map = std::map<std::string,int>() ;
    unit_num_map = std::map<std::string,int>() ;
    if(is_single_temperature(in_exp)!=0)
      calculate_temperature(in_exp, input_value);
    else{
      //cout<<endl;
      //cout<<"Loci expression:  ";
      //in_exp->Print(cout);
	
      /*cout<<endl;
        cout<<endl;
        cout<<"print the unit:  "<<endl;*/
      seperate_unit(unit_num_map,unit_den_map,in_exp); 
      
      /*cout<<endl;
        cout<<endl;
        cout<<"change the unit to basic type:"<<endl;*/
      
      //change the reference type and composite type to the basic type//
      get_conversion(unit_num_map,unit_den_map,conversion_factor);
      if(conversion_factor>0){
        //cout<<endl;
        //cout<<"Final unit conversion is: "<<endl;
        //cout<<"Convert  ";
        //in_exp->Print(cout);
        //cout<< "  to basic units: "<<endl;
        rem_dup(unit_num_map,unit_den_map);
        /*cout<<"conversion factor is: "<<conversion_factor<<endl;
          cout<<"Input value is: "<<get_value()<<endl;
          cout<<"final value is: "<<conversion_factor*val<<endl;*/
      }else
        cout<<"Not an unit!"<<endl;  
    }
  }

  //compare the two units, check if they are comparable
  bool UNIT_type::is_compatible(const std::string unit_str){
    UNIT_type sec_unit;
    exprP sec_exp = 0;
    std::map<std::string,int>::iterator fst_mi;
    std::map<std::string,int> fst_n_map=(*this).unit_num_map,fst_d_map=(*this).unit_den_map;
    sec_unit.mode=(*this).mode;
    sec_unit.unit_kind=(*this).unit_kind;
    sec_unit.input_value=1;
    sec_unit.gradient=0 ;
    sec_unit.secondGradient = 0 ;

    sec_exp=expression::create(unit_str);
    sec_unit.output(sec_exp);
    //cout<<sec_unit;

    std::map<std::string,int>::iterator sec_mi;
    std::map<std::string,int> sec_n_map=sec_unit.unit_num_map,sec_d_map=sec_unit.unit_den_map;
    bool num_flag=false,den_flag=false;

    //    std::cerr << "unit_num_map ="  ;
    //    for(auto i=unit_num_map.begin();i!=unit_num_map.end();++i)
    //      std::cerr << "(" << i->first << "," << i->second << ")" ;
    //    std::cerr << endl ;
    //    std::cerr << "unit_den_map ="  ;
    //    for(auto i=unit_den_map.begin();i!=unit_den_map.end();++i)
    //      std::cerr << "(" << i->first << "," << i->second << ")" ;
    //    std::cerr << endl ;
    for(fst_mi= fst_n_map.begin();fst_mi!=fst_n_map.end();++fst_mi) {
      for(sec_mi= sec_n_map.begin();sec_mi!=sec_n_map.end();++sec_mi) {
        if(((*sec_mi).first==(*fst_mi).first)&&((*sec_mi).second==(*fst_mi).second))
          num_flag=true;
      }
    }
    for(fst_mi= fst_d_map.begin();fst_mi!=fst_d_map.end();++fst_mi) {
      for(sec_mi= sec_d_map.begin();sec_mi!=sec_d_map.end();++sec_mi) {
        if(((*sec_mi).first==(*fst_mi).first)&&((*sec_mi).second==(*fst_mi).second))
          den_flag=true;
      }
    }
    if(sec_d_map.size()==0&&fst_d_map.size()==0)
      den_flag=true;
    if(sec_n_map.size()==0&&fst_n_map.size()==0)
      num_flag=true ;

    //    std::cerr << "num_flag=" << num_flag <<",den_flag=" << den_flag
    //              << ",sec_d_map.size()=" << sec_d_map.size()
    //              << ",fst_d_map.size()=" << fst_d_map.size()
    //              << ",sec_n_map.size()=" << sec_n_map.size()
    //              << ",fst_n_map.size()=" << fst_n_map.size() << std::endl ;
    if((num_flag==true)&&(den_flag==true)&&(sec_d_map.size()==fst_d_map.size())&&(sec_n_map.size()==fst_n_map.size()))
      return true;
    return false;
  }

  double UNIT_type::get_value_in(const std::string unit_str){
    UNIT_type sec_unit;
    exprP sec_exp = 0 ;

    sec_unit.mode=(*this).mode;
    sec_unit.unit_kind=(*this).unit_kind;

    sec_unit.input_value=1;
    sec_unit.gradient=0 ;
    sec_unit.secondGradient = 0 ;

    sec_exp=expression::create(unit_str);
    if(is_single_temperature(sec_exp)!=0){
      double val = input_value ;
      reverse_calculate_temperature(sec_exp,val);
      return val ;
    }
    else
      sec_unit.output(sec_exp);
    //cout<<sec_unit;

    return (input_value)*(*this).conversion_factor/sec_unit.conversion_factor;
  }

  FAD2d UNIT_type::get_value_inD(const std::string unit_str){
    UNIT_type sec_unit;
    exprP sec_exp = 0 ;

    sec_unit.mode=(*this).mode;
    sec_unit.unit_kind=(*this).unit_kind;

    sec_unit.input_value=1;
    sec_unit.gradient=0 ;
    sec_unit.secondGradient = 0 ;

    sec_exp=expression::create(unit_str);
    if(is_single_temperature(sec_exp)!=0){
      FAD2d val(input_value,gradient,secondGradient) ;
      reverse_calculate_temperature(sec_exp,val);
      return val;
    }
    else
      sec_unit.output(sec_exp);
    //cout<<sec_unit;
    double factor = (*this).conversion_factor/sec_unit.conversion_factor;
    return FAD2d(input_value,gradient,secondGradient)*factor;
  }

  VFAD UNIT_type::get_value_inVF(const std::string unit_str, int batch){
    UNIT_type sec_unit;

    sec_unit.mode=(*this).mode;
    sec_unit.unit_kind=(*this).unit_kind;

    sec_unit.input_value=1;
    sec_unit.gradient=0 ;
    sec_unit.secondGradient = 0 ;

    exprP sec_exp=expression::create(unit_str);
    VFAD val ;
    val.data.value = input_value ;
    for(size_t i=0;i<VFAD::maxN;++i)
      val.data.grad[i] = 0 ;
    if(batch == 0) {
      val.data.grad[0] = gradient ;
      for(size_t i=1;i<min(gradientList.size()+1,VFAD::maxN);++i)
        val.data.grad[i] = gradientList[i-1] ;
    } else {
      size_t offset = batch*VFAD::maxN ;
      size_t num = min(VFAD::maxN, offset-(gradientList.size()+1)) ;
      for(size_t i=0;i<num;++i)
        val.data.grad[i] = gradientList[i+offset-1] ;
    }
    
    if(is_single_temperature(sec_exp)!=0){
      reverse_calculate_temperature(sec_exp,val);
      return val;

    }
    else
      sec_unit.output(sec_exp);
    //cout<<sec_unit;

    val *= (*this).conversion_factor/sec_unit.conversion_factor;
    return val;
    
  }

  void UNIT_type::get_values_in(const std::string unit_str,
                                double &val, double &grad, double &secondGrad,
                                std::vector<double> &gradlist){
    UNIT_type sec_unit;

    sec_unit.mode=(*this).mode;
    sec_unit.unit_kind=(*this).unit_kind;

    sec_unit.input_value=1;
    sec_unit.gradient=0 ;
    sec_unit.secondGradient = 0 ;

    exprP sec_exp=expression::create(unit_str);

    if(is_single_temperature(sec_exp)!=0){
      FADd v(input_value,1.0) ;
      reverse_calculate_temperature(sec_exp,v);
      val = realToDouble(v) ;
      double coef = v.grad ;
      grad = gradient*coef ;
      secondGrad = secondGradient*coef ;
      gradlist = gradientList ;
      for(size_t i=0;i<gradlist.size();++i)
        gradlist[i] *= coef ;
      return ;
    }  else
      sec_unit.output(sec_exp);
    //cout<<sec_unit;
    double coef = conversion_factor/sec_unit.conversion_factor;
    val = coef*input_value ;
    grad = gradient*coef ;
    secondGrad = secondGradient*coef ;
    gradlist = gradientList ;
    for(size_t i=0;i<gradlist.size();++i)
      gradlist[i] *= coef ;
  }
  

}
