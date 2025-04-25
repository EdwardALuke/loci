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
#include <vector>
#include <queue>

#include "prism.h"

using std::queue;
using std::cerr;
using std::endl;
using std::cout;
using std::list;

std::vector<char>  merge_tri_face_p(const  std::vector<char>& cellPlan1, int dd1, char orientCode1);
  

// void  write_quad_inner_faces(const std::map<QuadFace*, NeibIndex>& faces,
//                                 int cell_offset, int& mxppf, ofstream& ofile){
//   int nfaces = faces.size();
//   if(nfaces == 0) return ;

//   // face2node is store in faceVertex 
//   std::list<int32> faceVertex;
    
//   for(std::map<QuadFace*, NeibIndex>::const_iterator pf = faces.begin();
//       pf != faces.end(); pf++){
//     //for each face, set up face2node 
//     pf->first->set_f2n(faceVertex); 
    
//     //update mxppf
//     mxppf = max(mxppf,static_cast<int>(faceVertex.size()));
    
//     //output the face
//     int  numppf = faceVertex.size();
//     ofile << numppf <<'\t' ;
//     for(std::list<int32>::const_iterator np = faceVertex.begin();
//         np != faceVertex.end(); np++){ofile << *np  <<'\t';} 
//     ofile << pf->second.c1 + cell_offset << '\t' << pf->second.c2 + cell_offset << endl;
    
//   }//finish all faces
// }
        
// void  write_general_inner_faces(const std::map<Face*, NeibIndex>& faces,
//                                 int cell_offset, int& mxppf, ofstream& ofile){
//   int nfaces = faces.size();
//   if(nfaces == 0) return ;

//   // face2node is store in faceVertex 
//   std::list<int32> faceVertex;
    
//   for(std::map<Face*, NeibIndex>::const_iterator pf = faces.begin();
//       pf != faces.end(); pf++){
//     //for each face, set up face2node 
//     pf->first->set_f2n(faceVertex); 
    
//     //update mxppf
//     mxppf = max(mxppf,static_cast<int>(faceVertex.size()));
    
//     //output the face
//     int  numppf = faceVertex.size();
//     ofile << numppf <<'\t' ;
//     for(std::list<int32>::const_iterator np = faceVertex.begin();
//         np != faceVertex.end(); np++){ofile << *np  <<'\t';} 
//     ofile << pf->second.c1 + cell_offset << '\t' << pf->second.c2 + cell_offset << endl;
    
//   }//finish all faces
// }

struct Prism_Quad{
  Prism* c;
  Range2d f;
  int faceID;
  Prism_Quad(Prism* cc, Range2d ff, int d):c(cc), f(ff), faceID(d){}
};
struct Prism_Face{
  Prism* c;
  Face* f;
  Prism_Face(Prism* cc, Face* ff):c(cc), f(ff){}
};

std::vector<int32> get_c1_prism(const std::vector<char>& cellPlan,
                                const std::vector<char>& facePlan,
                                char orientCode, int dd){

  std::vector<int32> c1;
 
  Prism *aCell = new Prism(); // the root, nfold = 3
  aCell->empty_resplit(cellPlan); //only nfold and childCell is defined

  
  if(dd >= 2){
  int64 maxX = int64(1) << MAXLEVEL;
  int64 maxY = int64(1) << MAXLEVEL;
  
  Range2d rg = Range2d(Point2d(0, 0), Point2d(maxX, maxY));

    
  //resplit face
  std::vector<Range2d> leaves;
  std::queue<Range2d> Q;
  Q.push(rg);
   

  Range2d currentF;
  unsigned int index =0;
  char mySplitCode;
  
  while(!Q.empty()){
     currentF = Q.front();
     
     if(index >= facePlan.size()){
       mySplitCode = 0;
     }
     else{ 
       //take a code from splitcode
       mySplitCode = facePlan[index];
       index++;  
     }
     
     
     //the middle points used for children range
     //       p4     pmax
     //
     //p1     p5     p3
     // 
     //pmin   p2
     
     
     Point2d p1, p2, p3, p4, p5, pmin, pmax;
     
     
         
     //define positions
     pmin = currentF.minP;
     pmax = currentF.maxP;
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
     
     
     //push the selected  children in the Q
     
     switch(mySplitCode){
     case 1:
       Q.push( Range2d(pmin,p3));
       Q.push( Range2d(p1, pmax));
       break;
       
     case 2:
       Q.push(Range2d(pmin, p4));
       Q.push(Range2d(p2, pmax));
       break;
       
     case 3:
       Q.push(Range2d(pmin, p5));
       Q.push(Range2d(p1, p4));
       Q.push(Range2d(p2, p3));
       Q.push(Range2d(p5, pmax));
       break;
     case 0:
       
       leaves.push_back(transfer_f2c(currentF,Point2d(maxX, maxY), orientCode));
       
       break;
     default:
       cerr << "illegal face code in resplit(): " << endl; 
       break;
     }
     
     Q.pop();
     
   }//while(!Q.empty())
   

 //resplit cell  
    
    if(cellPlan.size() != 0) {

      std::queue<Prism_Quad> Q;
      Q.push(Prism_Quad(aCell,rg, dd));
      std::vector<pair<Range2d, int32> > faceMap;
      
      Prism* current = 0;
      //   Range2d currentF;
    
      char cellCode, faceCode;
      int faceID;
      int nfold;
      while(!Q.empty()){
        current = Q.front().c;
        currentF = Q.front().f;
        cellCode = current->mySplitCode;
        nfold = current-> getNfold();
        faceID = Q.front().faceID;
      
     
      
        if(cellCode == 0){
          faceMap.push_back(make_pair(currentF,current->getCellIndex()));
        }
        
        //define chidID and child's faceID
        else{

           faceCode = cellCode;
          //the middle points used for children range
          //       p4     pmax
          //
          //p1     p5     p3
          // 
          //pmin   p2

          Point2d p1, p2, p3, p4, p5, pmin, pmax;
          
          //define positions
          pmin = currentF.minP;
          pmax = currentF.maxP;
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
         
      
          switch(faceCode){
          case 1:
            Q.push(Prism_Quad(current->childCell[0],Range2d(pmin, p3), faceID));
            Q.push(Prism_Quad(current->childCell[1], Range2d(p1, pmax), faceID));
            break;
        
          case 2:
            Q.push(Prism_Quad(current->childCell[faceID-2], Range2d(pmin, p4),2));
            Q.push(Prism_Quad(current->childCell[(faceID-2)== nfold-1? 0:(faceID-1)],Range2d(p2, pmax), 5));
            break;
      
          case 3:
            Q.push(Prism_Quad(current->childCell[faceID-2],Range2d(pmin, p5), 2));
            Q.push(Prism_Quad(current->childCell[faceID-2+nfold],Range2d(p1, p4), 2));
            Q.push(Prism_Quad(current->childCell[(faceID-2)== nfold-1? 0:(faceID-1)], Range2d(p2, p3),5));
            Q.push(Prism_Quad(current->childCell[(faceID-2)== nfold-1? nfold:(faceID-1+nfold)],Range2d(p5, pmax), 5));
            break;
          default:
            cerr<<"WARNING: illegal cellCode" << endl;
            break;
          }
        }//else
      
      Q.pop();
      }//while(Q.empty())
      c1 = contain_2d(faceMap, leaves);
     //  cout << "orientCode " << char(orientCode+'0') << "  dd " << dd << " prism facePlan: ";
//       for(int i = 0 ; i < facePlan.size(); i++){
//         cout << char(facePlan[i] + '0') << "  " ;
//       }
//       cout<< endl;

//       cout << "prism cellPlan: ";
//       for(int i = 0 ; i < cellPlan.size(); i++){
//         cout << char(cellPlan[i] + '0') << "  " ;
//       }
//       cout<< endl;

//       cout << "prism c1: ";
//       for(int i = 0 ; i < c1.size(); i++){
//         cout << c1[i] << "  " ;
//       }
//       cout<< endl;
      

      
    }//if(cellPlan.size() != 0)
    else{
  
      for(int i = 0 ; i <int(leaves.size()); i++){
        c1.push_back(int32(1));
      }
    } 
  }//if(dd>=2)
  else{ //dd < 2
    Face* aFace = new Face();
    aFace->numEdge = 3;
    std::vector<Face*> leaves;
    aFace->empty_resplit(facePlan, leaves);

  
    if(cellPlan.size() != 0){
      std::queue<Prism_Face> Q;
      Q.push(Prism_Face(aCell, aFace));
      
    
      Prism* current  =0;
      Face* currentF = 0;
      char cellCode;
     
      
      while(!Q.empty()){
        current = Q.front().c;
        currentF = Q.front().f;
       
        cellCode = current->mySplitCode;
        std::vector<Face*> local_leaves;

        //here if cellCode is put inside switch block, the c1 generated is not in the
        //same order as leaves
        while(cellCode == 1){
          current = current->childCell[dd];
          cellCode = current->mySplitCode;
        }
        switch(cellCode){
        case 0:
          if(currentF->child!=0){
            for(int i = 0; i < currentF->numEdge; i++){
              Q.push(Prism_Face(current, currentF->child[i]));
            }
          }
          else{
            c1.push_back(current->getCellIndex());
           
          }
       
          break;
        case 2:
          for(int i = 0; i <current->nfold; i++){
            Q.push(Prism_Face(current->childCell[general_childID_orient_f2c(i, orientCode, current->nfold)], currentF->child[i]));
          }
          break;
        case 3:
          for(int i = 0; i <current->nfold; i++){
            Q.push(Prism_Face(current->childCell[dd==0?general_childID_orient_f2c(i, orientCode,current->nfold):
                                                 (general_childID_orient_f2c(i, orientCode, current->nfold)+current->nfold)],
                              currentF->child[i]));
          }
          break;
        }
        Q.pop();
      }
      
    }//if(cellPlan.size()!=0)
    else{
      for(unsigned int i = 0; i< leaves.size(); i++){
        c1.push_back(int32(1));
      }
    }
    //clean up
    delete aFace;
    aFace = 0;
  }//else  
 
  //clean up
    delete aCell;
    reduce_vector(c1);
    return c1;
  }
  
