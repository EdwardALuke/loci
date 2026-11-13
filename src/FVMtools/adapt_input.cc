// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#include "adapt_input.h"
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>

namespace FVMTools {
namespace {
  constexpr int isotropicMode=2 ;

  int integer(const std::string &value) {
    std::size_t used=0 ; const int result=std::stoi(value,&used) ;
    if(used!=value.size()) throw std::invalid_argument("expected an integer: "+value) ;
    return result ;
  }
  double number(const std::string &value) {
    std::size_t used=0 ; const double result=std::stod(value,&used) ;
    if(used!=value.size() || !std::isfinite(result)) throw std::invalid_argument("expected a finite number: "+value) ;
    return result ;
  }
  void applyMode(Loci::RefinementOptions &options,int value) {
    switch(value) {
    case 0:
      options.automaticSplit=true ;
      options.extrusionAxis=-1 ;
      options.splitFamily="native" ;
      options.centerMethod="reference" ;
      break ;
    case 1:
      options.automaticSplit=false ;
      options.extrusionAxis=2 ;
      options.splitFamily="extruded" ;
      options.centerMethod="weighted" ;
      break ;
    case isotropicMode:
      options.automaticSplit=false ;
      options.extrusionAxis=-1 ;
      options.splitFamily="corner" ;
      options.centerMethod="weighted" ;
      break ;
    default:
      throw std::invalid_argument("-mode must be 0, 1 or 2") ;
    }
  }
}

  void applyQualityOptions(const std::map<std::string,std::string> &args,Loci::RefinementOptions &options) {
    auto &q=options.quality ;
    if(args.count("-quality")) q.mode=args.at("-quality") ;
    if(args.count("-quality-centroid")) q.centroid=args.at("-quality-centroid") ;
    if(args.count("-max-nonorthogonality")) q.maxNonOrthogonality=number(args.at("-max-nonorthogonality")) ;
    if(args.count("-max-face-warpage")) q.maxFaceWarpage=number(args.at("-max-face-warpage")) ;
    if(args.count("-max-child-volume-ratio")) q.maxChildVolumeRatio=number(args.at("-max-child-volume-ratio")) ;
    q.validate() ;
  }
  void showQualityHelp() {
    if(Loci::MPI_rank==0) std::cout<<
      "Quality: -quality off|report|reject (default off) -quality-centroid wireframe|exact\n"
      "  -max-nonorthogonality degrees -max-face-warpage degrees -max-child-volume-ratio ratio\n" ;
  }

  AdaptInput readAdaptInput(int argc,char **argv,bool meshOutput) {
    std::set<std::string> known={"-g","-r","-o","-ctag","-tag","-xml","-par","-mode","-family","-center","-axis","-automatic",
      "-levels","-tol","-balance","-grading","-difference","-factor","-max-depth","-max-nodes","-max-rounds","-quality","-quality-centroid","-max-nonorthogonality","-max-face-warpage","-max-child-volume-ratio"} ;
    if(meshOutput) { known.insert("-oplan") ; known.insert("-c2p") ; }
    std::map<std::string,std::string> args ;
    for(int i=1;i<argc;++i) {
      const std::string key=argv[i] ;
      if(!known.count(key) || i+1==argc || args.count(key)) throw std::invalid_argument("unknown, repeated or incomplete option: "+key) ;
      args[key]=argv[++i] ;
    }
    if(!args.count("-g")) throw std::invalid_argument("-g original.vog is required") ;
    int inputs=0 ; for(const char *key:{"-ctag","-tag","-xml","-par"}) inputs+=args.count(key) ;
    if(inputs!=1) throw std::invalid_argument("choose one of -ctag, -tag, -xml or -par") ;
    namespace F=Loci::FVMAdapt ; AdaptInput input ;
    auto &state=input.state ; auto &options=input.options ;
    applyMode(options,isotropicMode) ;
    input.mesh=args.at("-g") ; std::string error ;
    if(args.count("-r") && !F::readParametricPlan(args.at("-r"),args.at("-g"),state,options,error)) throw std::runtime_error(error) ;
    if(args.count("-mode")) applyMode(options,integer(args.at("-mode"))) ;
    if(args.count("-family")) {
      options.splitFamily=args.at("-family") ; options.automaticSplit=false ;
      options.centerMethod=(options.splitFamily=="native" || options.splitFamily=="mixed") ? "reference":"weighted" ;
      options.extrusionAxis=options.splitFamily=="extruded" ? 2:-1 ;
    }
    if(args.count("-center")) options.centerMethod=args.at("-center") ;
    if(args.count("-axis")) options.extrusionAxis=integer(args.at("-axis")) ;
    if(args.count("-automatic")) {
      const int automatic=integer(args.at("-automatic")) ;
      if(automatic<0 || automatic>1) throw std::invalid_argument("-automatic must be 0 or 1") ;
      options.automaticSplit=automatic!=0 ;
    }
    if(args.count("-grading")) options.grading=args.at("-grading") ;
    if(args.count("-balance")) {
      const int balance=integer(args.at("-balance")) ;
      if(balance<0 || balance>2) throw std::invalid_argument("-balance must be 0, 1 or 2") ;
      options.faceBalance=balance==0 ? "none":balance==1 ? "majority":"opposite" ;
    }
    if(args.count("-levels")) options.refinementLevels=integer(args.at("-levels")) ;
    if(args.count("-difference")) options.levelDifference=integer(args.at("-difference")) ;
    if(args.count("-max-depth")) options.maxDepth=integer(args.at("-max-depth")) ;
    if(args.count("-max-nodes")) options.maxNodes=integer(args.at("-max-nodes")) ;
    if(args.count("-max-rounds")) options.maxRounds=integer(args.at("-max-rounds")) ;
    if(args.count("-factor")) options.splitAspectRatio=number(args.at("-factor")) ;
    if(args.count("-tol")) options.minEdgeLength=number(args.at("-tol")) ;
    if(args.count("-xml")) options.regionFile=args.at("-xml") ;
    if(args.count("-par")) options.sourceFile=args.at("-par") ;
    options.tagLocation=args.count("-tag") ? "nodes":"cells" ;
    input.tagFile=args.count("-tag") ? args.at("-tag"):args.count("-ctag") ? args.at("-ctag"):"" ;
    if(!input.tagFile.empty() && args.count("-tol") && !options.automaticSplit)
      throw std::invalid_argument("-tol limits geometric or automatic requests; use -levels for isotropic tags") ;
    if(meshOutput) {
      input.outputMesh=args.count("-o") ? args.at("-o"):"out.vog" ;
      if(args.count("-oplan")) input.outputPlan=args.at("-oplan") ;
      if(args.count("-c2p")) input.cellParents=args.at("-c2p") ;
    } else input.outputPlan=args.count("-o") ? args.at("-o"):"out.plan" ;
    applyQualityOptions(args,options) ;
    return input ;
  }
  void showAdaptHelp(bool meshOutput) {
    if(Loci::MPI_rank!=0) return ;
    std::cout<<(meshOutput ?
      "refine -g original.vog [-r previous.plan] [-o out.vog] [-oplan next.plan] [-c2p parents.txt] <marking> [controls]\n":
      "marker -g original.vog [-r previous.plan] [-o out.plan] <marking> [controls]\n") ;
    std::cout<<
      "Marking: -ctag cells.txt | -tag nodes.txt | -xml region.xml | -par sources.par\n"
      "Tags: positive refines, negative coarsens, zero keeps; current mesh file order.\n"
      "Controls: -mode 0|1|2 (automatic native | planar Z | isotropic; default 2)\n"
      "  -family native|corner|mixed|extruded  -center reference|weighted|integrated\n"
      "  -automatic 0|1 (choose native hex/prism split directions; also available with mixed)\n"
      "  -axis 0|1|2  -levels N  -tol length (XML requires positive length)\n"
      "  -balance 0|1|2 (level grading | also majority | also opposite faces)\n"
      "  -grading faces|edges|vertices  -difference N  -factor ratio\n"
      "  -max-depth N  -max-nodes N  -max-rounds N\n"
      "Saved plans retain controls; command-line controls override them.\n" ;
    showQualityHelp() ;
    if(meshOutput) std::cout<<"Cell-parent pairs are one-based and refer to the input plan, or original cells without -r.\n" ;
  }
}
