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
#include "hex_defines.h"


using std::cerr;
using std::endl;
Range2d  transfer_f2c(Range2d f, Point2d maxPc, char orientCode){
  Point2d p1 = transfer_f2c(f.minP, maxPc, orientCode);
  Point2d p2 = transfer_f2c(f.maxP, maxPc, orientCode);
  return Range2d(Point2d(min(p1.x, p2.x), min(p1.y, p2.y)),
                 Point2d(max(p1.x, p2.x), max(p1.y, p2.y)));
}
  

  
Point2d  transfer_f2c(Point2d p, Point2d maxPc, char orientCode){
  
  Point2d pc;
  
  switch(orientCode){
  case 0:
    pc.x = p.x;
    pc.y = p.y;
    break;
    
  case 1:
    pc.x = p.y;
    pc.y = maxPc.y - p.x;
    break;
    
  case 2:
    pc.x = maxPc.x - p.x;
    pc.y = maxPc.y - p.y;
    break;
    
  case 3:
    pc.x = maxPc.x - p.y;
    pc.y = p.x;
    break;
    
  case 4:
    pc.x = p.y;
    pc.y = p.x;
    break;
    
  case 5:
    pc.x = maxPc.x - p.x;
    pc.y = p.y;
    break;
    
  case 6:
    pc.x = maxPc.x - p.y;
    pc.y = maxPc.y - p.x;
    break;

  case 7:
    pc.x = p.x;
    pc.y = maxPc.y - p.y;
    break;

  default:
    cerr <<"WARNING: illegal orientCode in function transfer_f2c" << endl;
    break;
  }
  return pc; 
}
            
Point2d transfer_c2f(Point2d p, Point2d maxPf, char orientCode){
  
  Point2d pf;
  
  switch(orientCode){
  case 0:
    pf.x = p.x;
    pf.y = p.y;
    break;
    
  case 1:
    pf.x = maxPf.x - p.y;
    pf.y = p.x;
    break;
    
  case 2:
    pf.x = maxPf.x - p.x;
    pf.y = maxPf.y - p.y;
    break;
    
  case 3:
    pf.x = p.y;
    pf.y = maxPf.y -p.x;
    break;
    
  case 4:
    pf.x = p.y;
    pf.y = p.x;
    break;
    
  case 5:
    pf.x = maxPf.x - p.x;
    pf.y = p.y;
    break;
    
  case 6:
    pf.x = maxPf.x - p.y;
    pf.y = maxPf.y - p.x;
    break;

  case 7:
    pf.x = p.x;
    pf.y = maxPf.y - p.y;
    break;

  default:
    cerr << "WARNING: illegal orientCode in function transfer_c2f "<<endl;
    
    break;
  }
  return pf; 
}




