//#############################################################################
//#
//# Copyright 2015-2019, Mississippi State University
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
#include <vector>
#include "defines.h"
#include "dataxferDB.h"
using std::cout;
using std::endl;
using std::string;
using std::ifstream;

namespace Loci {
  storeRepP Local2FileOrder(storeRepP sp, entitySet dom, int &offset,
                            fact_db::distribute_infoP dist, MPI_Comm comm) ;
  void File2LocalOrder(storeRepP &result, entitySet resultSet,
                       storeRepP input, int offset,
                       fact_db::distribute_infoP dist,
                       MPI_Comm comm) ;
}





class init_num_original_nodes : public unit_rule{
  param<int> num_original_nodes;
public:
  init_num_original_nodes(){
    name_store("num_original_nodes", num_original_nodes);
    output("num_original_nodes");
    constraint("UNIVERSE");
    
  }
  //parameter, no loop, 
  virtual void compute(const sequence &seq){
    
    
    *num_original_nodes = 0;
  }
}; 
register_rule<init_num_original_nodes> register_init_num_original_nodes;

class apply_num_original_nodes : public apply_rule<param<int>, Loci::Summation<int> >{
  param<int> num_original_nodes;
  const_store<vect3d> pos;
public:
  apply_num_original_nodes(){
    name_store("num_original_nodes", num_original_nodes);
    name_store("pos", pos);
    input("pos");
    input("num_original_nodes");
    output("num_original_nodes");
   
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
  }
  void calculate(Entity cc){
    *num_original_nodes +=1;
  }
}; 
register_rule<apply_num_original_nodes> register_apply_num_original_nodes;


class get_cellPlan : public pointwise_rule{
  const_param<std::string> planfile_par;
  store<std::vector<char> > cellPlan;
  
  public:
  get_cellPlan(){
    name_store("planfile_par", planfile_par);
    name_store("cellPlan", cellPlan);
    input("planfile_par");
    output("cellPlan");
    constraint("geom_cells");
    disable_threading();
  }
  virtual void compute(const sequence &seq){
       hid_t file_id;
       entitySet dom = entitySet(seq);
       file_id = H5Fopen((*planfile_par).c_str(), H5F_ACC_RDONLY,
                         H5P_DEFAULT);
       
       Loci::readContainer(file_id,"cellPlan",cellPlan.Rep(),dom) ;
       H5Fclose(file_id);
       do_loop(seq, this); 
      
  }
  void calculate(Entity cc){
    if(cellPlan[cc].size() == 1 && cellPlan[cc][0] == 'C') cellPlan[cc].resize(0);
  } 
};
register_rule<get_cellPlan> register_get_cellPlan;

class get_DBcellPlan : public pointwise_rule{
  const_param<std::string> planDB_par;
  store<std::vector<char> > cellPlan;
  
  public:
  get_DBcellPlan(){
    name_store("planDB_par", planDB_par);
    name_store("cellPlan", cellPlan);
    input("planDB_par");
    output("cellPlan");
    constraint("geom_cells");
    disable_threading();
  }
  virtual void compute(const sequence &seq){
    entitySet dom = entitySet(seq);
    store<std::vector<char> > readPlan ;

    readPlan = Loci::DataXFER_DB.getItem((*planDB_par).c_str()) ;
       
    fact_db::distribute_infoP dist = Loci::exec_current_fact_db->get_distribute_info() ;
    if(dist==0) {
      FORALL(dom,ii) {
	cellPlan[ii] = readPlan[ii] ;
      } ENDFORALL ;
    } else {
      store<std::vector<char> > tmp ;
      tmp.allocate(dom) ;
      Loci::storeRepP tmpRep = tmp.Rep() ;
      int offset = 0 ;
      File2LocalOrder(tmpRep,dom,readPlan.Rep(),offset,dist,MPI_COMM_WORLD) ;
      FORALL(dom,ii) {
	cellPlan[ii].swap(tmp[ii]) ;
      } ENDFORALL ;
    }
       
    FORALL(dom,cc) {
      if(cellPlan[cc].size() == 1 && cellPlan[cc][0] == 'C') 
	cellPlan[cc].resize(0);
    } ENDFORALL ;
  }
};
register_rule<get_DBcellPlan> register_get_DBcellPlan;


class get_balanced_cellPlanDB : public pointwise_rule{
  const_param<std::string> planDB_par;
  store<std::vector<char> > cellPlan;
  
  public:
  get_balanced_cellPlanDB(){
    name_store("balanced_planDB_par", planDB_par);
    name_store("priority::refmesh::balancedCellPlan", cellPlan);
    input("balanced_planDB_par");
    output("priority::refmesh::balancedCellPlan");
    constraint("geom_cells");
    disable_threading();
  }
  virtual void compute(const sequence &seq){
   
    entitySet dom = entitySet(seq);
    store<std::vector<char> > readPlan ;

    readPlan = Loci::DataXFER_DB.getItem((*planDB_par).c_str()) ;
       
    fact_db::distribute_infoP dist = Loci::exec_current_fact_db->get_distribute_info() ;
    if(dist==0) {
      FORALL(dom,ii) {
	cellPlan[ii] = readPlan[ii] ;
      } ENDFORALL ;
    } else {
      store<std::vector<char> > tmp ;
      tmp.allocate(dom) ;
      Loci::storeRepP tmpRep = tmp.Rep() ;
      int offset = 0 ;
      File2LocalOrder(tmpRep,dom,readPlan.Rep(),offset,dist,MPI_COMM_WORLD) ;
      FORALL(dom,ii) {
	cellPlan[ii].swap(tmp[ii]) ;
      } ENDFORALL ;
    }
       
    FORALL(dom,cc) {
      if(cellPlan[cc].size() == 1 && cellPlan[cc][0] == 'C') 
	cellPlan[cc].resize(0);
    } ENDFORALL ;
  }

};
register_rule<get_balanced_cellPlanDB> register_get_balanced_cellPlanDB;


class get_balanced_cellPlan : public pointwise_rule{
  const_param<std::string> planfile_par;
  store<std::vector<char> > cellPlan;
  
  public:
  get_balanced_cellPlan(){
    name_store("balanced_planfile_par", planfile_par);
    name_store("priority::refmesh::balancedCellPlan", cellPlan);
    input("balanced_planfile_par");
    output("priority::refmesh::balancedCellPlan");
    constraint("geom_cells");
    disable_threading();
  }
  virtual void compute(const sequence &seq){
   
   
    
       hid_t file_id;
       entitySet dom = entitySet(seq);
       file_id = H5Fopen((*planfile_par).c_str(), H5F_ACC_RDONLY,
                         H5P_DEFAULT);
       
       Loci::readContainer(file_id,"cellPlan",cellPlan.Rep(),dom) ;
       H5Fclose(file_id);
       do_loop(seq, this); 
      
  }
  void calculate(Entity cc){
    if(cellPlan[cc].size() == 1 && cellPlan[cc][0] == 'C') cellPlan[cc].resize(0);
  } 
};
register_rule<get_balanced_cellPlan> register_get_balanced_cellPlan;






//read in plan from former cycle
class get_parentPlan : public pointwise_rule{
  const_param<std::string> parent_planfile_par;
  store<std::vector<char> > parentPlan;
  
  public:
  get_parentPlan(){
    name_store("parent_planfile_par", parent_planfile_par);
    name_store("parentPlan", parentPlan);
    
    input("parent_planfile_par");
    output("parentPlan");
    constraint("geom_cells");
    disable_threading();
  }
  virtual void compute(const sequence &seq){
   
   
    
       hid_t file_id;
       entitySet dom = entitySet(seq);
       file_id = H5Fopen((*parent_planfile_par).c_str(), H5F_ACC_RDONLY,
                         H5P_DEFAULT);
       
       Loci::readContainer(file_id,"cellPlan",parentPlan.Rep(),dom) ;
       H5Fclose(file_id);
       do_loop(seq, this); 
      
  }
  void calculate(Entity cc){
    if(parentPlan[cc].size() == 1 && parentPlan[cc][0] == 'C') parentPlan[cc].resize(0);
  } 
};
register_rule<get_parentPlan> register_get_parentPlan;

//read in plan from former cycle
class get_parentPlanDB : public pointwise_rule{
  const_param<std::string> parent_planDB_par;
  store<std::vector<char> > parentPlan;
  
  public:
  get_parentPlanDB(){
    name_store("parent_planDB_par", parent_planDB_par);
    name_store("parentPlan", parentPlan);
    
    input("parent_planDB_par");
    output("parentPlan");
    constraint("geom_cells");
    disable_threading();
  }
  virtual void compute(const sequence &seq){
    
    entitySet dom = entitySet(seq);
    store<std::vector<char> > readPlan ;

    readPlan = Loci::DataXFER_DB.getItem((*parent_planDB_par).c_str()) ;
       
    fact_db::distribute_infoP dist = Loci::exec_current_fact_db->get_distribute_info() ;
    if(dist==0) {
      FORALL(dom,ii) {
	parentPlan[ii]= readPlan[ii] ;
      } ENDFORALL ;
    } else {
      store<std::vector<char> > tmp ;
      tmp.allocate(dom) ;
      Loci::storeRepP tmpRep = tmp.Rep() ;
      int offset = 0 ;
      File2LocalOrder(tmpRep,dom,readPlan.Rep(),offset,dist,MPI_COMM_WORLD) ;
      FORALL(dom,ii) {
	parentPlan[ii].swap(tmp[ii]) ;
      } ENDFORALL ;
    }
       
    FORALL(dom,cc) {
      if(parentPlan[cc].size() == 1 && parentPlan[cc][0] == 'C') 
	parentPlan[cc].resize(0);
    } ENDFORALL ;
  }
  void calculate(Entity cc){ } 
};
register_rule<get_parentPlanDB> register_get_parentPlanDB;


class reprocess_plan : public pointwise_rule {
  const_store<std::vector<char> > balancedCellPlan1 ;
  const_param<bool> beginWithMarker; //dummy parameter to trick Loci scheduler
  store<std::vector<char> > balancedCellPlan;
 public: 
  reprocess_plan(){
    name_store("balancedCellPlan1", balancedCellPlan1);
    name_store("balancedCellPlan", balancedCellPlan);
    name_store("beginWithMarker", beginWithMarker); 
    input("balancedCellPlan1");
    input("beginWithMarker");
    output("balancedCellPlan");
  }
  virtual void compute(const sequence &seq) {
    do_loop(seq, this);
  }
  
  void calculate(Entity cc){
       balancedCellPlan[cc] = balancedCellPlan1[cc];
  }
    

};
register_rule<reprocess_plan> register_reprocess_plan;


class process_plan : public pointwise_rule {
  const_store<std::vector<char> > balancedCellPlan ;
  store<std::vector<char> > tmpPlan;
 public: 
  process_plan(){
    name_store("balancedCellPlan", balancedCellPlan);
    name_store("tmpPlan", tmpPlan);
    input("balancedCellPlan");
    output("tmpPlan");
  }
  virtual void compute(const sequence &seq) {
    do_loop(seq, this);
  }
  
  void calculate(Entity cc){
    if(balancedCellPlan[cc].size() == 0) tmpPlan[cc].push_back('C');
    else tmpPlan[cc] = balancedCellPlan[cc];
  }
    

};
register_rule<process_plan> register_process_plan;
  
class cellplan_output_file : public pointwise_rule {
  const_store<std::vector<char> > tmpPlan ;
  const_param<string> outfile_par ;
  store<bool> cellplan_output ;
  
public:
  cellplan_output_file(){
    name_store("tmpPlan", tmpPlan);
    name_store("plan_outfile_par", outfile_par);
    name_store("cellplan_output", cellplan_output);
    input("(tmpPlan,plan_outfile_par)");
    output("cellplan_output");
    constraint("geom_cells");
    disable_threading();
  }
  virtual void compute(const sequence &seq) {
    
    hid_t file_id;
    
    file_id = Loci::hdf5CreateFile((*outfile_par).c_str(),H5F_ACC_TRUNC,
                                   H5P_DEFAULT, H5P_DEFAULT) ;
    
    if(Loci::MPI_rank == 0)cout << "writing cellplan into " << *outfile_par << endl;
    Loci::writeContainer(file_id,"cellPlan",tmpPlan.Rep()) ;
    Loci::hdf5CloseFile(file_id) ;
    if(Loci::MPI_rank == 0) cout << "Finish writing cellPlan " << endl;

  
    
  }
 
} ;
register_rule<cellplan_output_file> register_cellplan_output_file;

class cellplan_writeDB : public pointwise_rule {
  const_store<std::vector<char> > tmpPlan ;
  const_param<string> outDB_par ;
  store<bool> cellplan_output ;
  
public:
  cellplan_writeDB(){
    name_store("tmpPlan", tmpPlan);
    name_store("plan_outDB_par", outDB_par);
    name_store("cellplan_output", cellplan_output);
    input("(tmpPlan,plan_outDB_par)");
    output("cellplan_output");
    constraint("geom_cells");
    disable_threading();
  }
  virtual void compute(const sequence &seq) {
    entitySet dom = entitySet(seq);
    
    fact_db::distribute_infoP dist = Loci::exec_current_fact_db->get_distribute_info() ;
    if(dist==0) {
      store<std::vector<char> > pcopy ;
      pcopy.allocate(dom) ;
      FORALL(dom,ii) {
	pcopy[ii] = tmpPlan[ii] ;
      } ENDFORALL ;

      Loci::DataXFER_DB.insertItem((*outDB_par).c_str(),pcopy.Rep()) ;
    } else {
      int offset = 0 ;
      Loci::storeRepP vardist = Loci::Local2FileOrder(tmpPlan.Rep(),dom,offset,
						      dist,
						      MPI_COMM_WORLD) ;
      vardist->shift(offset) ;
      Loci::DataXFER_DB.insertItem((*outDB_par).c_str(),vardist) ;
    }
  }
 
} ;


register_rule<cellplan_writeDB> register_cellplan_writeDB;



// class cellweight_output_file : public pointwise_rule {
//   //this rule is not used by chem, cellweight_writeDB is actually used
//   const_store<int> num_fine_cells;
//   const_param<string> outfile_par ;
//   store<bool> cellweight_output ;
  
// public:
//   cellweight_output_file(){
//     name_store("cellweight_outfile_par", outfile_par);
//     name_store("cellweight_output", cellweight_output);
//     name_store("balanced_num_fine_cells", num_fine_cells);
//     input("(cellweight_outfile_par,balanced_num_fine_cells)");
//     output("cellweight_output");
//     constraint("geom_cells");
//     disable_threading();
//   }
//   virtual void compute(const sequence &seq) {
    
//     hid_t file_id;
//     file_id = Loci::hdf5CreateFile((*outfile_par).c_str(),H5F_ACC_TRUNC,
//                                    H5P_DEFAULT, H5P_DEFAULT) ;
//     Loci::writeContainer(file_id,"cellweights",num_fine_cells.Rep()) ;
//     Loci::hdf5CloseFile(file_id) ;
//   }
// } ;
// register_rule<cellweight_output_file> register_cellweight_output_file;

class cellweight_writeDB : public pointwise_rule {
  const_store<int>  num_fine_cells ;
  const_param<string> outDB_par ;
  store<bool> cellweight_output ;
  
public:
  cellweight_writeDB(){
    name_store("balanced_num_fine_cells", num_fine_cells);
    name_store("cellweight_outDB_par", outDB_par);
    name_store("cellweight_output", cellweight_output);
    input("(balanced_num_fine_cells,cellweight_outDB_par)");
    output("cellweight_output");
    constraint("geom_cells");
    disable_threading();
  }
  virtual void compute(const sequence &seq) {
    entitySet dom = entitySet(seq);
    
    fact_db::distribute_infoP dist = Loci::exec_current_fact_db->get_distribute_info() ;
    if(dist==0) {
      store<int> pcopy ;
      pcopy.allocate(dom) ;
      FORALL(dom,ii) {
	pcopy[ii] = num_fine_cells[ii] ;
      } ENDFORALL ;

      Loci::DataXFER_DB.insertItem((*outDB_par).c_str(),pcopy.Rep()) ;
    } else {
      int offset = 0 ;
      Loci::storeRepP vardist = Loci::Local2FileOrder(num_fine_cells.Rep(),dom,offset,
						      dist,
						      MPI_COMM_WORLD) ;
      vardist->shift(offset) ;
      Loci::DataXFER_DB.insertItem((*outDB_par).c_str(),vardist) ;
    }
  }
 
} ;
register_rule<cellweight_writeDB> register_cellweight_writeDB;

