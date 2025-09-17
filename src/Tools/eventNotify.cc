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
#include <Tools/eventNotify.h>

#include <vector>

namespace Loci {
  eventDispatcher::eventDispatcher() {}

  eventDispatcher::~eventDispatcher() {
    warn(notify_group.size()!=0) ;
  }

  void eventDispatcher::engage(eventNotify *p) {
    //bmutex l(mutex) ;
    mutex.lock();
    warn(!p) ;
    notify_group.push_back(p) ;
    mutex.unlock();
  }

  void eventDispatcher::disengage(eventNotify *p) {
    //bmutex l(mutex) ;
    mutex.lock();
    
    warn(!p) ;

    if(notify_group.begin() == notify_group.end()) {
      std::cerr << "disengage with empty list" << std::endl ;
      std::cerr << "p = " << p << std::endl ;
      return ;
    }
    notify_list::iterator nlp = notify_group.end() ;
    for(--nlp;nlp != notify_group.begin() && *nlp != p;--nlp)
      /* NULL STATEMENT */ ;

    warn(*nlp != p) ;
        
    if(*nlp == p)
      notify_group.erase(nlp) ;

    mutex.unlock();
  }    

  void eventDispatcher::dispatch_notify() {
    mutex.lock() ;
    // notify_list copy ;
    notify_list::iterator nlp ;
    for(nlp=notify_group.begin();nlp!=notify_group.end();++nlp)
      (*nlp)->notification();
      //copy.push_back(*nlp) ;
    mutex.unlock() ;
    // for(nlp=copy.begin();nlp!=copy.end();++nlp) 
    //   (*nlp)->notification() ;
  }

  eventNotify::~eventNotify() {}

}
