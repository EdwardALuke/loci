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
//##############################################################################

#ifndef FIELD_SORT_H
#define FIELD_SORT_H 1
#include <Tools/intervalSet.h>

namespace Loci{
 
  //used for gStore, gMultiStore etc.
  template<class T>  inline  bool fieldSort1(const std::pair<gEntity, T> &i1, 
                                             const std::pair<gEntity, T> &i2) {
    return i1.first < i2.first ;
  }

  template<class T>  inline  bool fieldSort1_unique(const std::pair<gEntity, T> &i1, 
                                                    const std::pair<gEntity, T> &i2) {
    if(i1.first < i2.first)return true ;
    else if(i1.first == i2.first)return i1.second < i2.second;
    return false;
  }
  //used for gMap, gMultiMap etc.
  inline bool fieldSort2(const std::pair<gEntity,gEntity> &i1,
                         const std::pair<gEntity,gEntity> &i2) {
    return i1.second < i2.second ;
  }

  //used  for general vector, not just for gContainers
  template<class T> inline  bool field_sort2(const std::pair<T,T> &p1,
                                             const std::pair<T,T> &p2) {
    return p1.second < p2.second ;
  }

  //sort the component according to different filed
  //gMultiStore component
  template<class T>  struct ms_comp{
    gEntity dom;
    T img;
    short ind;
    ms_comp(gEntity i1, T i2, short i3):dom(i1), img(i2), ind(i3){}
  };
  
  //sort according domain field and index
  template<class T> inline  bool field_sort_dom(const ms_comp<T> &i1, 
                                                const ms_comp<T> &i2) {
    if(i1.dom < i2.dom) return true ;
    if(i1.dom == i2.dom) return i1.ind < i2.ind;
    return false;
  }
  
  //sort according image field
   inline  bool field_sort_img(const ms_comp<gEntity> &i1, 
                               const ms_comp<gEntity> &i2) {
    return (i1.img < i2.img) ;
   }

  
  struct f2e_comp{
    gEntity f;
    gEntity n1;
    gEntity n2;
    short ind;
    f2e_comp(gEntity i1, gEntity i2, gEntity i3,  short i4):f(i1), n1(i2), n2(i3), ind(i4){}
  };
 //sort according domain field and index
  inline  bool f2e_field_sort_f(const f2e_comp &i1, 
                            const f2e_comp &i2) {
    if(i1.f < i2.f) return true;
    if(i1.f == i2.f) return i1.ind < i2.ind;
    return false;
  }
  inline  bool f2e_field_sort_e(const f2e_comp &i1, 
                            const f2e_comp &i2) {
    if(i1.n1 < i2.n1) return true;
    if(i1.n1 == i2.n1) return (i1.n2 < i2.n2);
    return false;
  }

  struct e2n_comp{
    gEntity e;
    gEntity n1;
    gEntity n2;
    e2n_comp(gEntity i1, gEntity i2, gEntity i3):e(i1), n1(i2), n2(i3){}
  };

  inline  bool e2n_field_sort_n(const e2n_comp &i1, 
                                const e2n_comp &i2) {
    if(i1.n1 < i2.n1) return true;
    if(i1.n1 == i2.n1) return (i1.n2 < i2.n2);
    return false;
  }


  
}


#endif
