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
#include <queue>
#include <vector>
#include <utility>
#include <list>
#include <algorithm>
#include <set>
#include <Loci.h>
#include "defines.h"
#include "hex_defines.h"
#include "quadface.h"

using std::queue;
using std::list;
using std::set;
using std::cerr;
using std::endl;

using namespace Loci;

class plannameFile_par : public singleton_rule {
  param<std::string> planname_par  ;
  const_param<std::string> planfile_par ;
public:
  plannameFile_par() {
    name_store("planname_par",planname_par) ;
    name_store("planfile_par",planfile_par) ;
    input("planfile_par") ;
    output("planname_par") ;
  }
  
  virtual void compute(const sequence &seq) {
    *planname_par = *planfile_par ;
  }
} ;

register_rule<plannameFile_par> register_plannameFile_par ;

class plannameDB_par : public singleton_rule {
  param<std::string> planname_par  ;
  const_param<std::string> planDB_par ;
public:
  plannameDB_par() {
    name_store("planname_par",planname_par) ;
    name_store("planDB_par",planDB_par) ;
    input("planDB_par") ;
    output("planname_par") ;
  }
  
  virtual void compute(const sequence &seq) {
    *planname_par = *planDB_par ;
  }
} ;
register_rule<plannameDB_par> register_plannameDB_par ;



class general_edge_points_unit : public unit_rule{
  store<SetLong> pointSet;
  const_param<std::string> planname_par;
public:
  general_edge_points_unit(){
    name_store("pointSet", pointSet);
    name_store("planname_par", planname_par);
    
    input("planname_par");
    output("pointSet");
    constraint("edge2node->pos");
    }
  void calculate(Entity e){
    
  };
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
  }
};



class general_edge_points_apply : public apply_rule<store<SetLong>, SetLongUnion>{
  const_store<std::vector<char> > facePlan;
  const_param<std::string> planname_par;
  const_multiMap face2edge;
  const_multiMap face2node;
  const_MapVec<2> edge2node;
  const_store<vect3d> pos;//dummy
  const_store<bool> is_quadface;
  store<SetLong> pointSet;
public:
  general_edge_points_apply(){
    name_store("facePlan", facePlan);
    name_store("planname_par", planname_par);
    name_store("face2edge", face2edge);
    name_store("pointSet", pointSet);
    name_store("face2node", face2node);
    name_store("edge2node", edge2node);
    name_store("pos", pos);
    name_store("is_quadface", is_quadface);
    
    input("planname_par");
    input("facePlan, is_quadface");
    input("face2node->pos");
    input("face2edge->edge2node->pos");
    input("face2edge->pointSet");
    output("face2edge->pointSet");
    constraint("faces");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
  }
  void calculate(Entity f){
    if(is_quadface[f]){

 if(facePlan[f].size() == 0) return;

    if(face2node.num_elems(f)!=4){
      cerr << "WARNING: Not a hexcell" << endl;
      Loci::Abort();
    }
    
    for( int i = 0; i < face2edge.num_elems(f); i++){

      int64 maxValue = int64(1) << MAXLEVEL;
      std::vector<char> plan;
      extract_quad_edge(facePlan[f], plan, i); 
      if(plan.size() == 0) continue;
      
      bool needReverse = false;
        
      if(i == 0 && edge2node[face2edge[f][0]][0] == face2node[f][1]) needReverse = true;
      else if(i == 1 && edge2node[face2edge[f][1]][0] == face2node[f][2]) needReverse = true;
      else if(i == 2 && edge2node[face2edge[f][2]][0] == face2node[f][2]) needReverse = true;
      else if(i == 3 && edge2node[face2edge[f][3]][0] == face2node[f][3]) needReverse = true;
      
      queue<pair<int64, int64> > Q;
      char code = plan.front();
      unsigned int index = 0;
      pair<int64, int64> pRange = make_pair(int64(0), maxValue);
      Q.push(pRange);
      
      int64 midP;
      
      while(!Q.empty()){
      //read in a code
        if(index >= plan.size()) code = 0;
        else code = plan[index++];
        
      //read in the range from the queue
        pRange = Q.front();
        
        //process the cell, push children ranges into Q
        //and put the middle point into points    
        switch(code){
        case 1:
          midP = (pRange.first + pRange.second)/2;
          Q.push(make_pair(pRange.first, midP));
          Q.push(make_pair(midP, pRange.second));
          
          if(needReverse) (pointSet[face2edge[f][i]].aset).insert(maxValue - midP);
          else (pointSet[face2edge[f][i]].aset).insert(midP);
          break;
          
        case 0:
          break;
          
        case 8:
          Q.push(pRange);
          break;
        default:
          cerr << "WARNING: illegal splitCode in rule edge_points_apply" << endl;
          break;
        }
        Q.pop();
      }
      plan.clear();
    }    
    }//is(is_quadface[f])
    else{
            
    if(facePlan[f].size() == 0) return;
    
    
    int64 maxValue = int64(1) << MAXLEVEL;
    
    //assume facePlan[f][0] ==1 
    int nfnode = face2node.num_elems(f);
    std::vector<bool> needReverse(nfnode);
  
    
    for(int i=0; i<nfnode; i++){
      if(edge2node[face2edge[f][i]][0] == face2node[f][i] &&
         edge2node[face2edge[f][i]][1] == face2node[f][i==(nfnode-1)?0:i+1]) needReverse[i] = false;
      else needReverse[i] = true;
    }
    
    
    //first time split face, put all the edgecenter into pointSet
    
   
    
    std::queue<pair<int, TwoEdge> > Q;
    std::vector<TwoEdge> child(nfnode);
    for(int i =0; i< nfnode; i++){
      if(needReverse[i]){
        child[i].e0 = make_pair(maxValue, maxValue/2);
      }
      else{
        child[i].e0 = make_pair(int64(0),maxValue/2);
      }
      if(needReverse[i==0?nfnode-1:i-1]){
        child[i].e3 = make_pair(maxValue/2, int64(0));
      }
      else{
        child[i].e3 = make_pair(maxValue/2, maxValue);
      }
      Q.push(make_pair(i, child[i]));
      (pointSet[face2edge[f][i]].aset).insert(maxValue/2);
    }
    
   
       
    TwoEdge current;
    unsigned  int index = 1;
    char mySplitCode;
    int edgeID;
    
    while(!Q.empty()){
      current = Q.front().second;
      edgeID = Q.front().first;
      if(index >= facePlan[f].size()){
        mySplitCode = 0;
      }
      else{ 
        //take a code from splitcode
        mySplitCode = facePlan[f][index];
        index++;  
      }
     
      if(mySplitCode == 1){
                
        //define new edgeID and push the children into Q
        //put edgecenter into pointSet if necessary
        int newID[4] = {-1,-1,-1,-1} ;
        //define the 4 children  
        TwoEdge child[4];


        //if edgeID== -1, each newID is -1, no edgecenter need to be put
        //into pointSet
        if(edgeID == -1){
          for(int i=0; i < 3; i++)newID[i] = -1;
        }
        //if edgeID is in [0, nfnode), edge 0 is on edgeID, edge3 is on edgeID-1 
        else if(edgeID >= 0 && edgeID < nfnode){
          newID[2] = -1;
          newID[0] = edgeID; //e0 on edgeID, e3 on edgeID-1
          child[0].e0 = make_pair(current.e0.first, (current.e0.first + current.e0.second)/2);
          child[0].e3 = make_pair((current.e3.first + current.e3.second)/2, current.e3.second);
          
          newID[1] = 2*nfnode+edgeID;//e3 on edgId
          child[1].e3 = make_pair((current.e0.first + current.e0.second)/2, current.e0.second);
          
          newID[3] = nfnode +(edgeID==0?nfnode-1:edgeID-1);//e0 on edgeId
          child[3].e0 = make_pair(current.e3.first, (current.e3.first + current.e3.second)/2);
          
          pointSet[face2edge[f][edgeID]].aset.insert((current.e0.first + current.e0.second)/2);
          
          pointSet[face2edge[f][edgeID==0?nfnode-1:edgeID-1]].aset.insert((current.e3.first + current.e3.second)/2);
          
        }
        //edge 0 on edgeID-nfnode
        else if(edgeID >=nfnode && edgeID <2*nfnode){
          newID[2] = newID[3] = -1;
          newID[0] = edgeID; //e0 on edgeId
          child[0].e0 = make_pair(current.e0.first, (current.e0.first + current.e0.second)/2);
          
          newID[1] = nfnode + edgeID;//e3 on edgeId
          child[1].e3 = make_pair((current.e0.first + current.e0.second)/2, current.e0.second);
          
          pointSet[face2edge[f][edgeID-nfnode]].aset.insert((current.e0.first + current.e0.second)/2);
        }
        //edge 3 on edgeID-2*nfnode
        else if(edgeID >= 2*nfnode && edgeID < 3*nfnode){
          newID[1] = newID[2] = -1;
          newID[0] = edgeID;//e3 on edgeId
          child[0].e3 = make_pair((current.e3.first + current.e3.second)/2, current.e3.second);
          
          newID[3] = edgeID - nfnode; //e0 on edgeId
          child[3].e0 = make_pair(current.e3.first, (current.e3.first + current.e3.second)/2);
          
          pointSet[face2edge[f][edgeID-2*nfnode]].aset.insert((current.e3.first + current.e3.second)/2);
        }
        else{
          cerr << " WARNING: illegal edgeID" << endl;
          Loci::Abort();
        }
        
        for(int i = 0; i < 4; i++){
          Q.push(make_pair(newID[i],child[i])); 
        }
        
        
      }
     
      Q.pop();
    }
    
  }
  }
  
};
register_rule<general_edge_points_unit> register_general_edge_points_unit;
register_rule<general_edge_points_apply> register_general_edge_points_apply;

class make_edge_plan:public pointwise_rule{
  const_store<SetLong> pointSet;
  store<std::vector<char> > edgePlan;
public:
  make_edge_plan(){
    name_store("pointSet", pointSet);
    name_store("edgePlan", edgePlan);
    input("pointSet");
    output("edgePlan");
    constraint("edge2node->pos");
  }
  virtual void compute(const sequence &seq){

    do_loop(seq, this);
  }
  void calculate(Entity e){
     if(pointSet[e].aset.empty())return;
    //define range of the integer coordinates of the edge
  
    int64 maxValue = int64(1) << MAXLEVEL;
    
    //prepare to make up the refinement plan according to list points
    edgePlan[e].clear();
    queue<pair<int64, int64> > Q;
    pair<int64,int64>  pRange = make_pair(int64(0), maxValue);
    Q.push(pRange);
    int64 midP;
    
    while(!Q.empty()){
      //read in the range
      pRange = Q.front();
      
      //define the positions
      midP = (pRange.first + pRange.second) / 2;
      
      //if midP is in pointSet,set the code to 1,
      //else the code is 0
      if(pointSet[e].aset.empty())break;
      else{
        std::set<int64>::iterator p = pointSet[e].aset.find(midP);
        if (p != pointSet[e].aset.end()){
          Q.push(make_pair(pRange.first, midP));
          Q.push(make_pair(midP, pRange.second));
          edgePlan[e].push_back(1);
       
        }
        else{
          edgePlan[e].push_back(0);
        }
      }
      Q.pop();
    }
    while(edgePlan[e].size() != 0 && edgePlan[e].back() == 0 )edgePlan[e].pop_back();
    reduce_vector(edgePlan[e]);
    
   
  }
};

register_rule<make_edge_plan> register_make_edge_plan;


  
