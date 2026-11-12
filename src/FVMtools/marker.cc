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

#include "adapt_input.h"
#include <iostream>
#include <stdexcept>

int main(int argc,char **argv) {
  Loci::Init(&argc,&argv) ; int status=0 ;
  try {
    if(argc==1 || (argc==2 && std::string(argv[1])=="-help")) FVMTools::showAdaptHelp(false) ;
    else {
      auto input=FVMTools::readAdaptInput(argc,argv,false) ;
      Loci::rule_db rules ; rules.add_rules(Loci::global_rule_list) ; Loci::load_module("fvmadapt",rules) ;
      namespace F=Loci::FVMAdapt ; std::string error ;
      if(!F::markParametricPlan(rules,input.state,input.mesh,input.options,input.tagFile,error) ||
         !F::writeParametricPlan(input.outputPlan,input.state,input.options,error)) throw std::runtime_error(error) ;
    }
  } catch(const std::exception &error) { if(Loci::MPI_rank==0) std::cerr<<"marker: "<<error.what()<<'\n' ; status=1 ; }
  Loci::Finalize() ; return status ;
}
