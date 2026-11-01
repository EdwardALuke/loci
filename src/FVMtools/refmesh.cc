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
#include <FVMAdapt/parametric_plan.h>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>

int main(int argc,char **argv) {
  Loci::Init(&argc,&argv) ; int status=0 ;
  try {
    if(argc==1 || (argc==2 && std::string(argv[1])=="-help")) {
      if(Loci::MPI_rank==0) std::cout<<
        "refmesh -g original.vog -r saved.plan [-o out.vog] [-c2p parents.txt] [-pr previous.plan]\n"
        "Cell-parent pairs are one-based. Without -pr, parents are original cells.\n"
        "Plans retain their split and geometry controls.\n" ;
    } else {
      const std::set<std::string> known={"-g","-r","-o","-c2p","-pr"} ; std::map<std::string,std::string> args ;
      for(int i=1;i<argc;++i) {
        const std::string key=argv[i] ;
        if(!known.count(key) || i+1==argc || args.count(key)) throw std::invalid_argument("unknown, repeated or incomplete option: "+key) ;
        args[key]=argv[++i] ;
      }
      if(!args.count("-g") || !args.count("-r")) throw std::invalid_argument("-g original.vog and -r saved.plan are required") ;
      if(args.count("-pr") && !args.count("-c2p")) throw std::invalid_argument("-pr requires -c2p") ;
      namespace F=Loci::FVMAdapt ; Loci::RefinementState state,parent ; Loci::RefinementOptions options,parentOptions ; std::string error ;
      if(!F::readParametricPlan(args.at("-r"),args.at("-g"),state,options,error)) throw std::runtime_error(error) ;
      if(args.count("-pr")) {
        if(!F::readParametricPlan(args.at("-pr"),args.at("-g"),parent,parentOptions,error)) throw std::runtime_error(error) ;
        if(options.splitFamily!=parentOptions.splitFamily || options.extrusionAxis!=parentOptions.extrusionAxis)
          throw std::invalid_argument("parent and output plans use different split families or extrusion axes") ;
      }
      Loci::rule_db rules ; rules.add_rules(Loci::global_rule_list) ; Loci::load_module("fvmadapt",rules) ;
      if(!F::writeParametricMesh(rules,state,args.count("-pr") ? &parent:nullptr,options,
          args.count("-o") ? args.at("-o"):"out.vog",args.count("-c2p") ? args.at("-c2p"):"",error)) throw std::runtime_error(error) ;
    }
  } catch(const std::exception &error) { if(Loci::MPI_rank==0) std::cerr<<"refmesh: "<<error.what()<<'\n' ; status=1 ; }
  Loci::Finalize() ; return status ;
}
