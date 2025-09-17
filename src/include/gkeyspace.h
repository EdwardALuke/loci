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
#ifndef GKEY_SPACE_H
#define GKEY_SPACE_H

#include <vector>
#include <map>
#include <set>
#include <string>
#include <mpi.h>
#include <Loci_types.h>
#include <Tools/intervalSet.h>
#include <Tools/cptr.h>
#include <Tools/expr.h>
#include <Tools/variable.h>
#include <Tools/simple_partition_long.h>
#include <gmap.h>
#include <store_rep.h>
#include <distribute_long.h>
#include <gkey_manager.h>
namespace Loci {

  
  enum KeyScope {FILESPACE, GLOBALSPACE, LOCALSPACE} ;
  class gKeySpace;
  typedef CPTR<gKeySpace> gKeySpaceP ;
  
  class gKeySpace: public CPTR_type {
  public:
    typedef CPTR<gKeySpace> gKeySpaceP ;
  protected:
    // the communicator that determines the group of
    // processess that share this set of keys
    MPI_Comm comm ;
    int rank ;                  // processes rank in comm
    int np ;             // communicator size

    // this is the key distribution across the processess
    // it is in file numbering
    //while key_ptn is set, keys and my_keys are also set
    std::vector<gEntitySet> key_ptn ;

    // the entire keyset over all processes
    gEntitySet keys ;
    
    // the keyset I own
    gEntitySet my_keys ;
    
    //the maps that map the keys between local, global and file numbering
    //everything related to local number is not coded yet, maybe only one map
    //is needed between local and global numbering since gMap can be sorted
    //according to either domain field or image field
    gMap g2f ;

    //the following two vectors are used in partition stage,
    //and since they are memory-consuming, should be cleaned up after partition is done 
    std::vector<gEntitySet> send_ptn;//the keys I send to each process
    std::vector<gEntitySet> recv_ptn;//the keys I receive from each process

    //is the current keyset in file numbering, global numbering or local numbering
    KeyScope scope;
    
    //the variables in out_var and invar need to be unique.
    //i.e., no two variable in out_var point to the same gStoreRepP
    variableSet out_vars; //the variable whose domain is this space
    variableSet in_vars; //the variable whose image is this space
  private:

    //singleton class, create and store gKeySpace, the name of gKeySpace are hard-coded,
    //it can be "NodeSpace", "FaceSpace", "CellSpace", "EdgeSpace", "BcSpace" and "UniverseSpace".
    //"BcSpace" is for boundary surfaces. And "UniverseSpace" is a special gKeySpace designed for
    //gParams, gConstraints and gBlackBoxs that are defined over all keys in all spaces. 
    class gKeySpaceManager{
      static std::map<std::string, gKeySpaceP> space_map;
    public:
      void add_space(const std::string& spacename, const std::string& casename, gKeySpaceP space){
        std::string name = spacename;
        if(casename !="") name += '.'+casename;
        std::map<std::string,gKeySpaceP> ::iterator mi = space_map.find(name);
        if(mi!= space_map.end()){
          debugout<<"WARNING: space " << name<< " exists, nothing done" << endl;
          // mi->second = space;
        }else{
          space_map[name] = space;
        }
      }
    
      static gKeySpaceP get_space(const std::string& spacename,const std::string& casename) {
        std::string name = spacename;
        if(casename !="") name += '.'+casename;
        std::map<std::string,gKeySpaceP> ::const_iterator mi = space_map.find(name);
        if(mi!= space_map.end())return mi->second;
        else{
          if(spacename=="UniverseSpace"){
            //std::cerr << "space " << name <<" does not exist, create one " << endl; 
            //return gKeySpaceP(0);
            gKeySpaceP new_space = new gKeySpace();
            new_space->register_space(spacename, casename);
            return new_space;
          }else{
            debugout << "space " << name <<" does not exist, return null space " << endl; 
            return gKeySpaceP(0);
          }
        }
      }

      static std::vector<gKeySpaceP> get_all_spaces() {
        std::vector<gKeySpaceP> result;
        for(std::map<std::string,gKeySpaceP> ::const_iterator mi = space_map.begin();
            mi != space_map.end(); mi++){
          std::string spacename = mi->first;
          std::size_t found = spacename.find("UniverseSpace");
          if(found == std::string::npos){
            result.push_back(mi->second);
          }
        }
        return result;
      }
    };
    static gKeySpaceManager *ksm ;
    void create_ksm() {
      if(0==ksm) { 
	ksm = new gKeySpace::gKeySpaceManager ;
      }
    }
      
  protected:    
    // set the keyspace's MPI_COMM_GROUP   
    void 
    MPI_comm_group(MPI_Comm c) {
      comm = c ;
      MPI_Comm_rank(comm, &rank) ;
      MPI_Comm_size(comm, &np) ;
    }
  private:
    // a copy function used in the copy constructor and
    // assignment operator
    void
    copy_from(const gKeySpace& ks) {
      comm = ks.comm ;
      rank = ks.rank ;
      np = ks.np ;
      key_ptn = ks.key_ptn ;
      keys = ks.keys ;
      my_keys = ks.my_keys;
      g2f = ks.g2f.clone(); 
    }
  public:
    //constructors
    gKeySpace()
    {
      create_ksm();
      MPI_comm_group(MPI_COMM_WORLD) ;
    }
     
    gKeySpace(const gKeySpace& ks) {create_ksm(); copy_from(ks) ;}
    gKeySpace& operator=(const gKeySpace& ks) {
      create_ksm();
      if(&ks != this)
        copy_from(ks) ;
      return *this ;
    }

    //inspectors for in_vars and out_vars
    variableSet get_in_vars() const{return in_vars;}
    variableSet get_out_vars() const{return out_vars;}

    //mutators for in_vars and out_vars
    void add_in_var(const variable& v){
      in_vars += v;
    }
    void remove_in_var(const variable& v){
      in_vars -= v;
    }
    void add_out_var(const variable& v){
      out_vars += v;
    }
    void remove_out_var(const variable& v){
      out_vars -= v;
    }

    //inspector for key_ptn
    const std::vector<gEntitySet>&
    get_key_ptn() const {return key_ptn ;}
    //mutator for key_ptn
    //whenever key_ptn is modified, keys and my_keys are also modified
    void
    set_key_ptn(const std::vector<gEntitySet>& ptn) {
      key_ptn = ptn ;
      keys = GEMPTY;
      for(unsigned int i = 0; i < ptn.size(); i++){
        keys += ptn[i];
      }
      my_keys = key_ptn[rank] ;
    }
    //inspector for keys
    const gEntitySet&
    get_keys() const {return keys ;}
    //inspector for my_keys
    const gEntitySet&
    get_my_keys() const {return my_keys ;}
    
    //inspectors for send_ptn and revc_ptn 
    const std::vector<gEntitySet>&
    get_send_ptn() const {return send_ptn;}
    const std::vector<gEntitySet>&
    get_recv_ptn() const {return recv_ptn;}
    //mutator for send_ptn and recv_ptn
    //ptn is send_ptn, since recv_ptn is always transposePtn of send_ptn
    //recv_ptn is modified whenever send_ptn is modified 
    void set_send_recv_ptn(const std::vector<gEntitySet>& ptn);
    
    
    static gKeySpaceP get_space(const std::string &spacename, const std::string& casename);
    static std::vector<gKeySpaceP> get_all_spaces();
    void register_space(const std::string &spacename, const std::string& casename) ;


    void set_g2f_map(gStoreRepP rep){g2f.setRep(rep);}
    
   
   
  
    //assume keys and my_keys are given in all processes
    //simple partition keys, and set send_ptn and recv_ptn
    void set_simple_partition(){
      vector<gEntitySet> temp_key_ptn = g_simple_partition<gEntity>(keys, comm);
      vector<gEntitySet> temp_send_ptn(np);
      for(unsigned int pid = 0; pid < temp_key_ptn.size(); pid++){ 
        temp_send_ptn[pid] = my_keys & temp_key_ptn[pid];
      }
      set_send_recv_ptn(temp_send_ptn); 
    }

    //this function is the same as set_simple_partition,
    //except that it outputs a vector specifying the owner of each entity in my_keys
    //assume keys and my_keys are given in all processes
    //simple partition keys, and set send_ptn and recv_ptn
    void set_simple_partition(vector<int>& procmap){
      vector<gEntitySet> temp_key_ptn = g_simple_partition<gEntity>(keys, comm);
      vector<gEntitySet> temp_send_ptn(np);
      procmap.resize(my_keys.size());
      int pid = 0;
      int index = 0;
      for(gEntitySet::const_iterator ei = my_keys.begin(); ei!= my_keys.end(); ei++){
        if(*ei > temp_key_ptn[pid].Max()) pid++; 
        temp_send_ptn[pid] += *ei;
        procmap[index++] = pid;
      }
      set_send_recv_ptn(temp_send_ptn);  
    }

    gKeySpaceP clone() const{
      return new gKeySpace(*this);
    }
    
    virtual ~gKeySpace() {}

    //inspectors for communitor
    MPI_Comm
    get_mpi_comm() const {return comm ;}
    int
    get_comm_rank() const {return rank ;}
    int
    get_np() const {return np ;}
    bool
    is_distributed() const {return np > 1 ;}
    

    //this method is called after recv_ptn is set
    //it generates global keys according to recv_ptn
    //and then it generates g2f_map
    void
    generate_key(gKeyManagerP& km){
      size_t num_keys = 0;
      for(int i = 0; i < np; i++){
        num_keys += recv_ptn[i].size();
      }
      //reset my_keys, 
      my_keys = km->generate_key(num_keys);
      generate_g2f_map();
      scope = GLOBALSPACE; 
    }

    //this method is used to generate ghost cells
    //num_keys: local num of keys need to generate
    //in g2f map, each newly generated key is mapped to itself 
    gEntitySet
    generate_key(gKeyManagerP& km, size_t num_keys){
      gEntitySet temp_keys = km->generate_key(num_keys);
      //add to my_keys
      my_keys += temp_keys;
      GFORALL(temp_keys, k){
        g2f.insert(k,k);
      }ENDGFORALL;
      scope = GLOBALSPACE ;
      return temp_keys;
    }

    
    void generate_g2f_map(){
      g2f.clear();
      g2f.reserve(my_keys.size());
      gEntitySet::const_iterator itr = my_keys.begin();
      for(int i=0;i<np;++i){
        for(gEntitySet::const_iterator itr2 =recv_ptn[i].begin();itr2 !=recv_ptn[i].end(); itr2++){
          gEntity file_num = *itr2; //file number
          gEntity global_num = *itr; //global number
          g2f.insert(global_num, file_num);
          itr++;
        }
      }
      g2f.local_sort();
    }

    //inspectors for maps
    const gStoreRepP
    get_f2g_map() const {return g2f.local_inverse();}
    const gMap&
    get_g2f_map() const {return g2f ;} 
  } ; // end of class gKeySpace
   
}// end of namespace Loci

#endif
