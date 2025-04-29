//#############################################################################
//#
//# Copyright 2008-2019, Mississippi State University
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
#include <map>
#include "defines.h"
#include <list>
using std::list ;
#include <string>
using std::string ;
#include <vector>
using std::vector ;
using std::pair ;
#include <algorithm>
using std::sort ;
using std::unique ;


using std::cout ;


using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::ofstream;
using Loci::MPI_processes;

class fileNumRule: public pointwise_rule {
  store<int> fileNumX ;
public:
  fileNumRule() {
    name_store("fileNumber(X)",fileNumX) ;
    output("fileNumber(X)") ;
    constraint("X") ;
    disable_threading() ;
  }
  void compute(const sequence &seq) {
   if(Loci::MPI_processes == 1) {
    
     for(sequence::const_iterator si=seq.begin();si!= seq.end();++si){
       fileNumX[*si] = *si ;
       // if(*si < 0) cout << "negative file Number " << *si << endl; 
     }
     return;
   }
   fact_db::distribute_infoP df = Loci::exec_current_fact_db->get_distribute_info() ;
    Map l2g ;
    l2g = df->l2g.Rep() ;
    store<unsigned char> key_domain ;
    key_domain = df->key_domain.Rep() ;
    //    dMap g2f ;
    //    g2f = df->g2f.Rep() ;

    for(sequence::const_iterator si=seq.begin();si!= seq.end();++si){
      int kd = key_domain[*si] ;
      fileNumX[*si] = df->g2fv[kd][l2g[*si]] ;
      //if(fileNumX[*si] < 0) cout << "negative file Number " << fileNumX[*si] << endl; 
    }
  }
  
} ;

register_rule<fileNumRule> register_fileNumRule ;






