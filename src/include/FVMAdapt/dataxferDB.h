//#############################################################################
//#
//# Copyright 2015-2025, Mississippi State University
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
#ifndef DATAXFERDB_H
#define DATAXFERDB_H

#include <Loci>
#include <map>
#include <vector>
#include <string>

namespace Loci {

  class DataXFERDatabase {
    std::map<std::string,Loci::storeRepP> dbmap ;
  public:
    Loci::storeRepP getItem(const string &s) {
      std::map<std::string,Loci::storeRepP>::iterator mi ;
      mi = dbmap.find(s) ;
      if(mi == dbmap.end()) 
	return 0 ;
      return mi->second ;
    }
    void insertItem(const string &s, Loci::storeRepP item) {
      dbmap[s] = item ;
    }
    void deleteItem(const string &s) {
      std::map<std::string,Loci::storeRepP>::iterator mi ;
      mi = dbmap.find(s) ;
      if(mi != dbmap.end()) 
	dbmap.erase(mi) ;
    }
    void clearDB() {
      dbmap = std::map<std::string,Loci::storeRepP>() ;
    }
    std::vector<std::string> nameList() const {
      std::vector<std::string> nlist ;
      std::map<std::string,Loci::storeRepP>::const_iterator mi ;
      for(mi=dbmap.begin();mi!=dbmap.end();++mi) {
	nlist.push_back(mi->first) ;
      }
      return nlist ;
    }
  } ;

  extern DataXFERDatabase DataXFER_DB ;
}
#endif
