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
//************************************************************************
// this file extract the face refinement plan from the cell refinement
//plan. A queue is used to simulate the tree-building process but no tree is
//actually built.
// Two tables are used. childrenID table indicates the IDs of the children
//whose codes need to be extracted. faceCodeTable transfers from cell code to
//face code.
//***************************************************************************

#include <vector>
#include <queue>
#include "hexcell.h"
#include "tables.h"
using std::queue;
using std::cerr;
using std::endl;
using std::cout;
using std::list;
//using namespace std;


std::vector<char>  extract_hex_face(const  std::vector<char>& cellPlan,  DIRECTION dd){

  //output faceCodeTable
  /*  for(int i=0; i<42; i++){
    if(i%7==0) cout<< endl;
    cout << faceCodeTable[i] <<',';
  }
  cout << endl << endl;
  
  //output faceIDTable
  for(int i=0; i<42; i++){
    if(i%7 == 0)cout <<endl; 
    cout <<'{';
    for(int j=0; j<faceIDTable[i].size(); j++) cout<<faceIDTable[i][j]<<',';
    cout <<'}';
   
  }
  cout << endl;
  
  cout << " cellPlan : " << endl;
  for(unsigned int i =0; i < cellPlan.size(); i++){
    cout << cellPlan[i] << " " ;
  }
  cout << endl;
  */

  std::vector<char> facePlan;
  
  if(cellPlan.size() == 0) {
    facePlan.clear();
    return facePlan;
  }
  unsigned int index = 0;  
  queue<bool> Q;
  Q.push(true);
  bool needExtract; 
  
  char cellCode, faceCode;
  while(!Q.empty()){
    needExtract = Q.front();
    if(index >= cellPlan.size()){
      break;
    }
    else{
      cellCode = cellPlan[index++];
    }
    //    cout << "cellCode: " << cellCode <<endl;
    if(needExtract){
      if(cellCode == 0){
        faceCode = 0;
      }
      else{
        faceCode = faceCodeTable[dd*7+cellCode- 1];
      }
     
      facePlan.push_back(faceCode);
      
      // cout << "cell: " << cellCode <<"  " <<"face: " << faceCode <<endl;
    }
    
    
    if(cellCode != 0){
      std::vector<bool> childrenID = faceIDTable[dd*7+cellCode-1];
      if(!needExtract){
        for(unsigned int i = 0; i < childrenID.size(); i++)childrenID[i] = 0;
      }
      
      for(unsigned int i = 0; i < childrenID.size(); i++){
        Q.push(childrenID[i]);
      }
    }
    Q.pop();
  }
  //delete 0 and 8 in the beginning of faceplan
  while(facePlan.size() != 0 && ((facePlan.front() == 0)||(facePlan.front() == 8))){
    facePlan.erase(facePlan.begin());}
  //delete 0 and 8 at the end of faceplan
  while(facePlan.size() != 0 && ((facePlan.back() == 0)||(facePlan.back() == 8))){
    facePlan.pop_back();
  }
  return facePlan;
}


std::vector<char> merge_quad_face(std::vector<char>& facePlanL, char orientCodeL,
                                  std::vector<char>& facePlanR, char orientCodeR){
  std::vector<char> plan;
  
  list<Edge2d> xEdge;//edges in x direction, will contain edges from both flPlan and frPlan 
  list<Edge2d> yEdge; // edges in y direction, will contain edges from both flPlan and frPlan
    
  int64 maxX = int64(1) << MAXLEVEL;
  int64 maxY = int64(1) << MAXLEVEL;
    
  Range2d rg = Range2d(Point2d(0, 0), Point2d(maxX, maxY));
  
  char orientCode = orientCodeL;
    
  for(int i = 1; i <= 2; i++){
    //first process flPlan, then process frPlan
    
    if( i == 1){
      plan.clear();
      plan = facePlanL;
      if(plan.size() == 0) continue;
      orientCode = orientCodeL;
    }
    if(i == 2){
      plan.clear();
      plan = facePlanR;
      if(plan.size() == 0) continue;
      orientCode = orientCodeR;
      }
    
      
      queue<Range2d> Q;
      char code = plan.front();
      unsigned int index = 0;
      Range2d pRange = rg;
      Q.push(pRange);
    

    //the middle points used for children range
    //       p4     pmax
    //
    //p1     p5     p3
    // 
    //pmin   p2
    
    
    Point2d p1, p2, p3, p4, p5, pmin, pmax;
    
    
    while(!Q.empty()){
      //read in a code
      if(index >= plan.size()) code = 0;
      else code = plan[index++];

      //read in the range from the queue
      pRange = Q.front();
      
       //define positions
      pmin = pRange.minP;
      pmax = pRange.maxP;
      p1.x = pmin.x;
      p1.y = (pmin.y + pmax.y)/2;
      p2.x = (pmin.x + pmax.x)/2;
      p2.y = pmin.y;
      p3.x = pmax.x;
      p3.y = p1.y;
      p4.x = p2.x;
      p4.y = pmax.y;
      p5.x = p2.x;
      p5.y = p1.y;

     
      //process the cell, push children ranges into Q
      //and put the inner edges into xEdges and yEdges    
      switch(code){
        
      case 3:
        Q.push(Range2d(pmin, p5));
        Q.push(Range2d(p1, p4));
        Q.push(Range2d(p2, p3));
        Q.push(Range2d(p5, pmax));
        
        switch(orientCode){
        case 0: //x->x, y->y
          xEdge.push_back(Edge2d(p1.y, p1.x, p3.x));
          yEdge.push_back(Edge2d(p2.x, p2.y, p4.y));
          break;

        case 1: //x->y, maxX-y ->x
          xEdge.push_back(Edge2d(p4.x, maxX - p4.y, maxX - p2.y));
          yEdge.push_back(Edge2d(maxX-p1.y, p1.x, p3.x));
          break;
          
        case 2: //maxX-x -> x, maxY-y -> y
          xEdge.push_back(Edge2d(maxY - p3.y, maxX - p3.x, maxX - p1.x));
          yEdge.push_back(Edge2d(maxX - p4.x, maxY - p4.y, maxY - p2.y));
          break;

        case 3: //y->x, maxY-x ->y
          xEdge.push_back(Edge2d(maxY - p2.x, p2.y, p4.y));
          yEdge.push_back(Edge2d(p3.y, maxY - p3.x, maxY - p1.x));
          break;

        case 4: // y-> x, x -> y
          xEdge.push_back(Edge2d(p2.x, p2.y, p4.y));
          yEdge.push_back(Edge2d(p1.y, p1.x, p3.x));
          break;
          
        case 5: //maxX-x ->x, y->y
          xEdge.push_back(Edge2d(p3.y, maxX - p3.x, maxX - p1.x));
          yEdge.push_back(Edge2d(maxX - p2.x, p2.y, p4.y));
          break;
  
        case 6://maxX-y ->x, maxY-x ->y
          xEdge.push_back(Edge2d(maxY - p4.x, maxX - p4.y, maxX - p2.y));
          yEdge.push_back(Edge2d(maxX - p3.y, maxY - p3.x, maxY - p1.x));
          break;
  
        case 7: //x ->x, maxY -y ->y
          xEdge.push_back(Edge2d(maxY - p1.y, p1.x, p3.x));
          yEdge.push_back(Edge2d(p4.x, maxY - p4.y, maxY - p2.y));
          break;
          
        default:
          cerr << "WARNING: illegal orientCode in rule merge_interior_face" << endl;
          break;
        }
        break;
        
      case 2:
        Q.push(Range2d(pmin, p4));
        Q.push(Range2d(p2, pmax));
        switch(orientCode){
        case 0: //x->x, y->y
          yEdge.push_back(Edge2d(p2.x, p2.y, p4.y));
          break;
          
        case 1: //x->y, maxX-y ->x
          xEdge.push_back(Edge2d(p4.x, maxX - p4.y, maxX - p2.y));
          break;
          
        case 2: //maxX-x -> x, maxY-y -> y
          yEdge.push_back(Edge2d(maxX - p4.x, maxY - p4.y, maxY - p2.y));
          break;
          
        case 3: //y->x, maxY-x ->y
          xEdge.push_back(Edge2d(maxY - p2.x, p2.y, p4.y));
          break;
          
        case 4: // y-> x, x -> y
          xEdge.push_back(Edge2d(p2.x, p2.y, p4.y));
          break;
          
        case 5: //maxX-x ->x, y->y
          yEdge.push_back(Edge2d(maxX - p2.x, p2.y, p4.y));
          break;
          
        case 6://maxX-y ->x, maxY-x ->y
          xEdge.push_back(Edge2d(maxY - p4.x, maxX - p4.y, maxX - p2.y));
          break;
  
        case 7: //x ->x, maxY -y ->y
          yEdge.push_back(Edge2d(p4.x, maxY - p4.y, maxY - p2.y));
          break;

        default:
          cerr << "WARNING: illegal orientCode in rule merge_interior_face" << endl;
          break;
        }
        
        break;
        
      case 1:
        Q.push(Range2d(pmin,p3));
        Q.push(Range2d(p1, pmax));
        switch(orientCode){
        case 0: //x->x, y->y
          xEdge.push_back(Edge2d(p1.y, p1.x, p3.x));
          break;
          
        case 1: //x->y, maxX-y ->x
          yEdge.push_back(Edge2d(maxX-p1.y, p1.x, p3.x));
          break;
          
        case 2: //maxX-x -> x, maxY-y -> y
          xEdge.push_back(Edge2d(maxY - p3.y, maxX - p3.x, maxX - p1.x));
          break;
          
        case 3: //y->x, maxY-x ->y
          yEdge.push_back(Edge2d(p3.y, maxY - p3.x, maxY - p1.x));
          break;

        case 4: // y-> x, x -> y
          yEdge.push_back(Edge2d(p1.y, p1.x, p3.x));
          break;
          
        case 5: //maxX-x ->x, y->y
          xEdge.push_back(Edge2d(p3.y, maxX - p3.x, maxX - p1.x));
          break;
          
        case 6://maxX-y ->x, maxY-x ->y
          yEdge.push_back(Edge2d(maxX - p3.y, maxY - p3.x, maxY - p1.x));
          break;
          
        case 7: //x ->x, maxY -y ->y
          xEdge.push_back(Edge2d(maxY - p1.y, p1.x, p3.x));
          break;

        default:
          cerr << "WARNING: illegal orientCode in rule merge_interior_face" << endl;
          break;
        }
        break;
        
      case 0:
        break;
        
      case 8:
        Q.push(pRange);
        break;

      default:
        cerr << "WARNING: illegal splitCode in rule merge_interior_face" << endl;
        break;
      }
      Q.pop();
    }
  }
  
  //sort and unique the edges
  if(xEdge.size() != 0){
    xEdge.sort();
    xEdge.unique();
  }
  if(yEdge.size() != 0){
    yEdge.sort();
    yEdge.unique();
  }
  
  //process the edges, if two edges overlap, they are replaced by their union
  if(xEdge.size() != 0){
    for(list<Edge2d>::iterator p1 = xEdge.begin(); p1 != xEdge.end(); p1++){
      list<Edge2d>::iterator p2 = p1;
      p2++;
      if(p2 == xEdge.end()) break;
      if(p1->pos == p2->pos){
        if(p2->head <= p1->tail){
          p1->tail = max(p1->tail, p2->tail);
          xEdge.erase(p2);
          p1--;
        }
      }
    }
  }
  if(yEdge.size() != 0){
    for(list<Edge2d>::iterator p1 = yEdge.begin(); p1 != yEdge.end(); p1++){
      list<Edge2d>::iterator p2 = p1;
      p2++;
      if(p2 == yEdge.end()) break;
      if(p1->pos == p2->pos){
        if(p2->head <= p1->tail){
          p1->tail = max(p1->tail, p2->tail);
          yEdge.erase(p2);
          p1--;
        }
      }
    }
  }
  
  //prepare to make up the refinement plan according to the edges
  plan.clear();
  if(xEdge.size()==0 && yEdge.size() == 0) return plan;
  queue<Range2d> Q;
  Range2d pRange = rg;
  Q.push(pRange);
    
  
  //the middle points used for children range
  //       p4     pmax
  //
  //p1     p5     p3
  // 
  //pmin   p2
  Point2d p1, p2, p3, p4, p5, pmin, pmax;
  bitset<2> binaryCode;
  
  while(!Q.empty()){
    //read in the range
    pRange = Q.front();
    
    //define the positions
    pmin = pRange.minP;
    pmax = pRange.maxP;
    p1.x = pmin.x;
    p1.y = (pmin.y + pmax.y)/2;
    p2.x = (pmin.x + pmax.x)/2;
    p2.y = pmin.y;
    p3.x = pmax.x;
    p3.y = p1.y;
    p4.x = p2.x;
    p4.y = pmax.y;
    p5.x = p2.x;
    p5.y = p1.y;
    
    //search the inner edges, if exist, delete the edge or the cut it off,
    //if p1->p3 exists in xEdge, binaryCode.set(0)
    //if p2->p4 exists in yEdge, binaryCode.set(1)
    binaryCode.reset();
    if(xEdge.size() != 0){
      list<Edge2d>::iterator p = std::find(xEdge.begin(), xEdge.end(), Edge2d(p1.y, p1.x, p3.x));
      if(p != xEdge.end()){
        binaryCode.set(0);
        xEdge.erase(p);
      }
      else{
        for(p = xEdge.begin(); p != xEdge.end(); p++){
          if(p->pos == p1.y && p->head <= p1.x && p->tail >= p3.x){
              binaryCode.set(0);
              if((p->head < p1.x) && (p->tail > p3.x)){
                int64 tempHead = p->head;
                p->head = p3.x;
                xEdge.insert(p, Edge2d(p1.y, tempHead, p1.x));
              }
              else if((p->head < p1.x) && (p->tail == p3.x)){
                p->tail = p1.x;}
              else{
                p->head = p3.x;
              }
              
              break;
            }
        }
      }
    }  
    
    
    if(yEdge.size() != 0){
      list<Edge2d>::iterator p = std::find(yEdge.begin(), yEdge.end(), Edge2d(p2.x, p2.y, p4.y));
      if(p != yEdge.end()){
        binaryCode.set(1);
        yEdge.erase(p);
      }
      else{
        for( p = yEdge.begin(); p != yEdge.end(); p++){
          if(p->pos == p2.x && p->head <= p2.y && p->tail >= p4.y){
            binaryCode.set(1);
            if((p->head < p2.y) && (p->tail > p4.y)){
              int64 tempHead = p->head;
              p->head = p4.y;
              yEdge.insert(p, Edge2d(p2.x, tempHead, p2.y));
            }
              else if((p->head < p2.y) && (p->tail == p4.y)){
                p->tail = p2.y;}
            else{
              p->head = p4.y;
            }
            break;
          }
        }
        
      }
    }
    
    //push the children range into Q
    switch(binaryCode.to_ulong()){
    case 3:
      Q.push(Range2d(pmin, p5));
      Q.push(Range2d(p1, p4));
      Q.push(Range2d(p2, p3));
        Q.push(Range2d(p5, pmax));
        plan.push_back(3);
        break;
        
    case 2:
        Q.push(Range2d(pmin, p4));
        Q.push(Range2d(p2, pmax));
        plan.push_back(2);
        break;
        
    case 1:
      Q.push(Range2d(pmin,p3));
      Q.push(Range2d(p1, pmax));
      plan.push_back(1);
      break;
        
    case 0:
      plan.push_back(0);
      break; 
    }
    Q.pop();
  }
  while(plan.size() != 0 && plan.back() == 0 )plan.pop_back();
  return plan;    
  
}


std::vector<char> merge_quad_face(std::vector<char>& facePlanL, char orientCodeL){
   std::vector<char> plan = facePlanL;
   
   
   if(plan.size() == 0)return plan;
    
    list<Edge2d> xEdge;//edges in x direction
    list<Edge2d> yEdge; // edges in y direction
    
    int64 maxX = int64(1) << MAXLEVEL;
    int64 maxY = int64(1) << MAXLEVEL;
    
    Range2d rg = Range2d(Point2d(0, 0), Point2d(maxX, maxY));
  
    char orientCode=orientCodeL;
   
    queue<Range2d> Q;
    char code = plan.front();
    unsigned int index = 0;
    Range2d pRange = rg;
    Q.push(pRange);
    

    //the middle points used for children range
    //       p4     pmax
    //
    //p1     p5     p3
    // 
    //pmin   p2
    
    
    Point2d p1, p2, p3, p4, p5, pmin, pmax;
    
    
    while(!Q.empty()){
      //read in a code
      if(index >= plan.size()) code = 0;
      else code = plan[index++];

      //read in the range from the queue
      pRange = Q.front();
      
       //define positions
      pmin = pRange.minP;
      pmax = pRange.maxP;
      p1.x = pmin.x;
      p1.y = (pmin.y + pmax.y)/2;
      p2.x = (pmin.x + pmax.x)/2;
      p2.y = pmin.y;
      p3.x = pmax.x;
      p3.y = p1.y;
      p4.x = p2.x;
      p4.y = pmax.y;
      p5.x = p2.x;
      p5.y = p1.y;

     
      //process the cell, push children ranges into Q
      //and put the inner edges into xEdges and yEdges    
      switch(code){
        
      case 3:
        Q.push(Range2d(pmin, p5));
        Q.push(Range2d(p1, p4));
        Q.push(Range2d(p2, p3));
        Q.push(Range2d(p5, pmax));
        
        switch(orientCode){
        case 0: //x->x, y->y
          xEdge.push_back(Edge2d(p1.y, p1.x, p3.x));
          yEdge.push_back(Edge2d(p2.x, p2.y, p4.y));
          break;

        case 1: //x->y, maxX-y ->x
          xEdge.push_back(Edge2d(p4.x, maxX - p4.y, maxX - p2.y));
          yEdge.push_back(Edge2d(maxX-p1.y, p1.x, p3.x));
          break;
          
        case 2: //maxX-x -> x, maxY-y -> y
          xEdge.push_back(Edge2d(maxY - p3.y, maxX - p3.x, maxX - p1.x));
          yEdge.push_back(Edge2d(maxX - p4.x, maxY - p4.y, maxY - p2.y));
          break;

        case 3: //y->x, maxY-x ->y
          xEdge.push_back(Edge2d(maxY - p2.x, p2.y, p4.y));
          yEdge.push_back(Edge2d(p3.y, maxY - p3.x, maxY - p1.x));
          break;

        case 4: // y-> x, x -> y
          xEdge.push_back(Edge2d(p2.x, p2.y, p4.y));
          yEdge.push_back(Edge2d(p1.y, p1.x, p3.x));
          break;
          
        case 5: //maxX-x ->x, y->y
          xEdge.push_back(Edge2d(p3.y, maxX - p3.x, maxX - p1.x));
          yEdge.push_back(Edge2d(maxX - p2.x, p2.y, p4.y));
          break;
  
        case 6://maxX-y ->x, maxY-x ->y
          xEdge.push_back(Edge2d(maxY - p4.x, maxX - p4.y, maxX - p2.y));
          yEdge.push_back(Edge2d(maxX - p3.y, maxY - p3.x, maxY - p1.x));
          break;
  
        case 7: //x ->x, maxY -y ->y
          xEdge.push_back(Edge2d(maxY - p1.y, p1.x, p3.x));
          yEdge.push_back(Edge2d(p4.x, maxY - p4.y, maxY - p2.y));
          break;
          
        default:
          cerr << "WARNING: illegal orientCode in rule merge_boundary_face" << endl;
          break;
        }
        break;
        
      case 2:
        Q.push(Range2d(pmin, p4));
        Q.push(Range2d(p2, pmax));
        switch(orientCode){
        case 0: //x->x, y->y
          yEdge.push_back(Edge2d(p2.x, p2.y, p4.y));
          break;
          
        case 1: //x->y, maxX-y ->x
          xEdge.push_back(Edge2d(p4.x, maxX - p4.y, maxX - p2.y));
          break;
          
        case 2: //maxX-x -> x, maxY-y -> y
          yEdge.push_back(Edge2d(maxX - p4.x, maxY - p4.y, maxY - p2.y));
          break;
          
        case 3: //y->x, maxY-x ->y
          xEdge.push_back(Edge2d(maxY - p2.x, p2.y, p4.y));
          break;
          
        case 4: // y-> x, x -> y
          xEdge.push_back(Edge2d(p2.x, p2.y, p4.y));
          break;
          
        case 5: //maxX-x ->x, y->y
          yEdge.push_back(Edge2d(maxX - p2.x, p2.y, p4.y));
          break;
          
        case 6://maxX-y ->x, maxY-x ->y
          xEdge.push_back(Edge2d(maxY - p4.x, maxX - p4.y, maxX - p2.y));
          break;
  
        case 7: //x ->x, maxY -y ->y
          yEdge.push_back(Edge2d(p4.x, maxY - p4.y, maxY - p2.y));
          break;

        default:
          cerr << "WARNING: illegal orientCode in rule merge_boundary_face" << endl;
          break;
        }
        
        break;
        
      case 1:
        Q.push(Range2d(pmin,p3));
        Q.push(Range2d(p1, pmax));
        switch(orientCode){
        case 0: //x->x, y->y
          xEdge.push_back(Edge2d(p1.y, p1.x, p3.x));
          break;
          
        case 1: //x->y, maxX-y ->x
          yEdge.push_back(Edge2d(maxX-p1.y, p1.x, p3.x));
          break;
          
        case 2: //maxX-x -> x, maxY-y -> y
          xEdge.push_back(Edge2d(maxY - p3.y, maxX - p3.x, maxX - p1.x));
          break;
          
        case 3: //y->x, maxY-x ->y
          yEdge.push_back(Edge2d(p3.y, maxY - p3.x, maxY - p1.x));
          break;

        case 4: // y-> x, x -> y
          yEdge.push_back(Edge2d(p1.y, p1.x, p3.x));
          break;
          
        case 5: //maxX-x ->x, y->y
          xEdge.push_back(Edge2d(p3.y, maxX - p3.x, maxX - p1.x));
          break;
          
        case 6://maxX-y ->x, maxY-x ->y
          yEdge.push_back(Edge2d(maxX - p3.y, maxY - p3.x, maxY - p1.x));
          break;
          
        case 7: //x ->x, maxY -y ->y
          xEdge.push_back(Edge2d(maxY - p1.y, p1.x, p3.x));
          break;

        default:
          cerr << "WARNING: illegal orientCode in rule merge_boundary_face" << endl;
          break;
        }
        break;
        
      case 0:
        break;
        
      case 8:
        Q.push(pRange);
        break;

      default:
        cerr << "WARNING: illegal splitCode in rule merge_boundary_face" << endl;
        break;
      }
      Q.pop();
    }
    
    
    //sort and unique the edges
    if(xEdge.size() != 0){
      xEdge.sort();
      xEdge.unique();
    }
    if(yEdge.size() != 0){
      yEdge.sort();
      yEdge.unique();
    }
    
    //process the edges, if two edges overlap, they are replaced by their union
    if(xEdge.size() != 0){
      for(list<Edge2d>::iterator p1 = xEdge.begin(); p1 != xEdge.end(); p1++){
        list<Edge2d>::iterator p2 = p1;
        p2++;
      if(p2 == xEdge.end()) break;
      if(p1->pos == p2->pos){
        if(p2->head <= p1->tail){
          p1->tail = max(p1->tail, p2->tail);
          xEdge.erase(p2);
          p1--;
        }
      }
      }
    }
    if(yEdge.size() != 0){
      for(list<Edge2d>::iterator p1 = yEdge.begin(); p1 != yEdge.end(); p1++){
        list<Edge2d>::iterator p2 = p1;
        p2++;
      if(p2 == yEdge.end()) break;
      if(p1->pos == p2->pos){
        if(p2->head <= p1->tail){
          p1->tail = max(p1->tail, p2->tail);
          yEdge.erase(p2);
          p1--;
        }
      }
      }
    }
    
    //prepare to make up the refinement plan according to the edges
    plan.clear();
    if(xEdge.size()==0 && yEdge.size() == 0) return plan;
    while(!Q.empty())Q.pop();
    pRange = rg;
    Q.push(pRange);
    
    
    //the middle points used for children range
  //       p4     pmax
  //
  //p1     p5     p3
  // 
  //pmin   p2
 
  bitset<2> binaryCode;
  
    while(!Q.empty()){
      //read in the range
      pRange = Q.front();
      
      //define the positions
      pmin = pRange.minP;
      pmax = pRange.maxP;
      p1.x = pmin.x;
      p1.y = (pmin.y + pmax.y)/2;
      p2.x = (pmin.x + pmax.x)/2;
      p2.y = pmin.y;
      p3.x = pmax.x;
      p3.y = p1.y;
      p4.x = p2.x;
      p4.y = pmax.y;
      p5.x = p2.x;
      p5.y = p1.y;

      //search the inner edges, if exist, delete the edge or the cut it off,
      //if p1->p3 exists in xEdge, binaryCode.set(0)
      //if p2->p4 exists in yEdge, binaryCode.set(1)
      binaryCode.reset();
      if(xEdge.size() != 0){
        list<Edge2d>::iterator p = std::find(xEdge.begin(), xEdge.end(), Edge2d(p1.y, p1.x, p3.x));
        if(p != xEdge.end()){
          binaryCode.set(0);
          xEdge.erase(p);
        }
        else{
          for(p = xEdge.begin(); p != xEdge.end(); p++){
            if(p->pos == p1.y && p->head <= p1.x && p->tail >= p3.x){
              binaryCode.set(0);
              if((p->head < p1.x) && (p->tail > p3.x)){
                int64 tempHead = p->head;
                p->head = p3.x;
                xEdge.insert(p, Edge2d(p1.y, tempHead, p1.x));
              }
              else if((p->head < p1.x) && (p->tail == p3.x)){
                p->tail = p1.x;}
              else{
                p->head = p3.x;
              }
            
              break;
            }
          }
        }
      }  
      
      
      if(yEdge.size() != 0){
        list<Edge2d>::iterator p = std::find(yEdge.begin(), yEdge.end(), Edge2d(p2.x, p2.y, p4.y));
        if(p != yEdge.end()){
          binaryCode.set(1);
          yEdge.erase(p);
        }
        else{
          for( p = yEdge.begin(); p != yEdge.end(); p++){
            if(p->pos == p2.x && p->head <= p2.y && p->tail >= p4.y){
              binaryCode.set(1);
              if((p->head < p2.y) && (p->tail > p4.y)){
                int64 tempHead = p->head;
                p->head = p4.y;
                yEdge.insert(p, Edge2d(p2.x, tempHead, p2.y));
              }
              else if((p->head < p2.y) && (p->tail == p4.y)){
                p->tail = p2.y;}
              else{
                p->head = p4.y;
              }
              break;
            }
          }
          
        }
      }
      
      //push the children range into Q
      switch(binaryCode.to_ulong()){
      case 3:
        Q.push(Range2d(pmin, p5));
        Q.push(Range2d(p1, p4));
        Q.push(Range2d(p2, p3));
        Q.push(Range2d(p5, pmax));
        plan.push_back(3);
        break;
        
      case 2:
        Q.push(Range2d(pmin, p4));
        Q.push(Range2d(p2, pmax));
        plan.push_back(2);
        break;
        
      case 1:
        Q.push(Range2d(pmin,p3));
        Q.push(Range2d(p1, pmax));
        plan.push_back(1);
        break;
        
      case 0:
        plan.push_back(0);
        break; 
      }
      Q.pop();
    }
    while(plan.size() != 0 && plan.back() == 0 )plan.pop_back();
    return plan;    
}
