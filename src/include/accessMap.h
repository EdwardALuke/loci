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
#ifndef ACCESSMAP_H
#define ACCESSMAP_H


// Note, this works for apply rules that are really doing the work of a blackbox
// It probably needs to be thought through more carefully, but it is needed
// when coupling interpolation routines between solvers
namespace Loci {

  class accessMapRepI : public MapRep {
  private:
    entitySet imageMap ;
  public:
    accessMapRepI() { }
    accessMapRepI(const entitySet &p) {imageMap=p; }
    virtual void allocate(const entitySet &ptn) {imageMap=ptn;}
    virtual ~accessMapRepI() {}
    virtual storeRep *new_store(const entitySet &p) const {
      return new accessMapRepI(p) ; }
    virtual storeRep *new_store(const entitySet &p, const int* cnt) const {
      return new accessMapRepI(p) ; }
    virtual storeRepP remap(const dMap &m) const
    { entitySet newImage = m.image(m.domain()&imageMap) ;
      return new accessMapRepI(newImage) ;
    } 
    virtual storeRepP MapRemap(const dMap &dm, const dMap &rm) const
    { entitySet newImage = dm.image(dm.domain()&imageMap) ;
      return new accessMapRepI(newImage) ;
    } 
    virtual void compose(const dMap &m, const entitySet &context) {}
    virtual void copy(storeRepP &st, const entitySet &context)
    { st = new_store(imageMap) ; }
    virtual void gather(const dMap &m, storeRepP &st,
                        const entitySet &context) 
    {st=new_store(imageMap) ; }
    virtual void scatter(const dMap &m, storeRepP &st,
                         const entitySet &context) 
    {st=new_store(imageMap) ; }
    
    virtual int pack_size(const entitySet& e, entitySet& packed) {return 0;}
    virtual int pack_size(const entitySet &e) {return 0 ; }
    virtual int estimated_pack_size(const entitySet &e) {return 0 ; }
    virtual void pack(void *ptr, int &loc, int &size, const entitySet &e) {}
    virtual void unpack(void *ptr, int &loc, int &size,  const sequence &seq)  {}
    
    virtual entitySet domain() const { return ~EMPTY ; } 

    virtual entitySet image(const entitySet &domain) const { return imageMap ; }
    virtual std::pair<entitySet,entitySet>
    preimage(const entitySet &codomain) const
    {entitySet d=~EMPTY; return std::make_pair(d,d) ; }
    virtual storeRepP get_map() {
      cerr << "getMap does not make sense for an accessMap" << endl ;
      multiMap result ;
      return result.Rep() ;
    }
    virtual std::ostream &Print(std::ostream &s) const {return s ;}
    virtual std::istream &Input(std::istream &s) {return s ;}
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en) {}
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en) const {}
#ifdef H5_HAVE_PARALLEL
    virtual void readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en, hid_t xfer_plist_id) {}
    virtual void writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en, hid_t xfer_plist_id) const {}
#endif
    virtual storeRepP expand(entitySet &out_of_dom, std::vector<entitySet> &init_ptn) {return new_store(EMPTY); }
    virtual storeRepP freeze() {return getRep() ; }
    virtual storeRepP thaw() {return getRep() ; }
    virtual DatatypeP getType() { return DatatypeP(new AtomicType(INT)) ; }


    virtual frame_info get_frame_info() {frame_info fi; return fi ; }
  } ;

  class const_accessMap ;
  class accessMap : public store_instance {
    friend  class const_accessMap ;
    typedef accessMapRepI accessMapType ;
    accessMap(const accessMap &var) { setRep(var.Rep()) ; }
    accessMap & operator=(const accessMap &str)
      { setRep(str.Rep()) ; return *this ;}
  public:
    accessMap() { setRep(new accessMapType) ; }
    accessMap(storeRepP rp) { setRep(rp) ; }
    virtual ~accessMap() {}
    virtual void notification() {}
    accessMap & operator=(storeRepP p) { setRep(p) ; return *this ;}
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    entitySet domain() const { return Rep()->domain() ; }
    entitySet image(const entitySet &dom) const {
      return MapRepP(Rep())->image(dom) ;
    }
    std::pair<entitySet,entitySet> preimage(const entitySet &codomain) const {
      return MapRepP(Rep())->preimage(codomain) ;
    }
    operator MapRepP() { MapRepP p(Rep()) ; fatal(p==0) ; return p ; }
    void elem(int indx) { } 
    void const_elem(int indx)  const {} 
    void operator[](int indx) { }
    void operator[](int indx) const { }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }
  } ;

  class const_accessMap : public store_instance {
    typedef accessMapRepI accessMapType ;
    const_accessMap(const const_accessMap &var) { setRep(var.Rep()) ; } 
    const_accessMap(const accessMap &var) { setRep(var.Rep()) ; }
    const_accessMap & operator=(const const_accessMap &str)
      { setRep(str.Rep()) ; return *this ;}
    const_accessMap & operator=(const accessMap &str)
    { setRep(str.Rep()) ; return *this ;}
  public:
    const_accessMap() { setRep(new accessMapType) ; }
    const_accessMap(storeRepP rp) { setRep(rp) ; }
    virtual ~const_accessMap() {}
    virtual void notification() {}
    virtual instance_type access() const { return READ_ONLY;}
    const_accessMap & operator=(storeRepP p) { setRep(p) ; return *this ;}
    const entitySet domain() const { return Rep()->domain() ; }
    entitySet image(const entitySet &dom) const {
      return MapRepP(Rep())->image(dom) ;
    }
    std::pair<entitySet,entitySet> preimage(const entitySet &codomain) const {
      return MapRepP(Rep())->preimage(codomain) ;
    }
    operator MapRepP() { MapRepP p(Rep()) ; fatal(p==0) ; return p ; }
    void const_elem(int indx)  const { }
    void operator[](int indx) const {  }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
  } ;  



  inline std::ostream & operator<<(std::ostream &s, const accessMap &m){
    return m.Print(s) ;
  }
    
  inline std::istream & operator>>(std::istream &s, accessMap &m) {
    return m.Input(s) ;
  }
  
  //*************************************************************************/

  inline std::ostream & operator<<(std::ostream &s,
				   const const_accessMap &m){
    return m.Print(s) ;
  }



}
    

#endif
