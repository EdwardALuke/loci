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
#include <vector>
#include <set>
#include <fstream>
#include <iostream>
#include <list>
#include <algorithm>
#include <map>
#include "hexcell.h"
#include "tables.h"

using std::cerr;
using std::cout;
using std:: endl;
using std::ofstream;

struct Cell_QuadFace {
  Cell_QuadFace(HexCell* cc, Range2d ff):c(cc), f(ff){};
  HexCell* c;
  Range2d f;
};


std::vector<int32> contain_2d(const std::vector<pair<Range2d, int32> >& faceMap,
                              const std::vector<Range2d>& leaves);


//this function will define face2node for each fine faces and write them out,
//at the same time, mxppf(max num of points per face) will be updated

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
        




std::vector<int32> get_c1_hex(const std::vector<char>& cellPlan,
                              const std::vector<char>& facePlan,
                              char orientCode,
                              char findex){
   
  std::vector<int32> c1;
  int64 maxX = int64(1) << MAXLEVEL;
  int64 maxY = int64(1) << MAXLEVEL;
  
  Range2d rg = Range2d(Point2d(0, 0), Point2d(maxX, maxY));
  
  

  HexCell* aCell = new HexCell();
  // split according to cellPlan,
  //only child is defined
  aCell->empty_resplit(cellPlan);
  
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
 
 if(cellPlan.size() != 0){
   
   std::queue<Cell_QuadFace> Q;
   Q.push(Cell_QuadFace(aCell,rg));
   std::vector<pair<Range2d, int32> > faceMap;
   
   HexCell* current = 0;
   //   Range2d currentF;
   
  
   while(!Q.empty()){
     current = Q.front().c;
  
     currentF = Q.front().f;
     
     if(current->mySplitCode == 0){
      
         faceMap.push_back(make_pair(currentF,current->getCellIndex()));
     }
   
     else{
       //indicate for each child of currentF, the childID of current
         std::vector<int> face2cellID; 
     
         //the code extract from cell
         char faceCode = faceCodeTable[findex*7+(current->mySplitCode)-1];
      
         //indicate how many children current has,and whether each child influence currentF 
         std::vector<bool> childrenID = faceIDTable[findex*7+(current->mySplitCode)-1];

         for(int i = 0; i<int(childrenID.size()); i++){
           if(childrenID[i]){
             face2cellID.push_back(i);
             
           }
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
         
         switch(faceCode){
           
         case 1:
           Q.push(Cell_QuadFace(current->childCell[face2cellID[0]],  Range2d(pmin,p3)));
           Q.push(Cell_QuadFace(current->childCell[face2cellID[1]], Range2d(p1, pmax)));
           break;

         case 2:
           Q.push(Cell_QuadFace(current->childCell[face2cellID[0]], Range2d(pmin, p4)));
           Q.push(Cell_QuadFace(current->childCell[face2cellID[1]], Range2d(p2, pmax)));
           break;
           
         case 3:
           Q.push(Cell_QuadFace(current->childCell[face2cellID[0]], Range2d(pmin, p5)));
           Q.push(Cell_QuadFace(current->childCell[face2cellID[1]], Range2d(p1, p4)));
           Q.push(Cell_QuadFace(current->childCell[face2cellID[2]], Range2d(p2, p3)));
           Q.push(Cell_QuadFace(current->childCell[face2cellID[3]], Range2d(p5, pmax)));
           break;
         case 8:
           if(face2cellID.size()>1){
             cerr<< "WARNING: face2cellID size " <<face2cellID.size()<< endl;
           }
           Q.push(Cell_QuadFace(current->childCell[face2cellID[0]], currentF)); 
           break;
         default:
           cerr << "illegal face code in resplit(): " << endl; 
           break;
         }
     }//else
     Q.pop();
      
   }//while(!Q.empty())
  
 
   
   c1 = contain_2d(faceMap,leaves);
   
 
}//end of if(cellPlan.size() != 0)
 else{
   
   for(int i = 0 ; i <int(leaves.size()); i++){
     c1.push_back(int32(1));
   }
 }
  //clean up
 if(aCell != 0){
   delete aCell;
   aCell = 0;
 }

 
 reduce_vector(c1);

 return c1;
 
}


