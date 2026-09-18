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
#include <queue>
#include "read_par.h"
#include "defines.h"
#include <limits>
#include <fstream>

using std::cout;
using std::endl;
using std::vector;
using std::ifstream;

/**
 * @file read_par.cc
 * @brief Reads source-based refinement parameters and evaluates local spacing
 *        requests.
 */

/**
 * Reads source refinement parameters from a plain text file.
 *
 * The file starts with the number of sources. Each source then provides two
 * endpoints followed by `r0`, `s0`, `r1`, `s1`, and `a`, matching the fields of
 * source_par.
 *
 * @param filename Path to the source-parameter file.
 * @param sources  Output vector replaced with the sources read from the file.
 */
void readPar(string filename, vector<source_par>& sources){
  sources.clear();
 
  ifstream inFile;
  
  inFile.open(filename.c_str());
  if(!inFile){
    cerr <<"can not open " << filename << " for input" << endl;
    Loci::Abort();
  }
  
  int numSources = 0;
  inFile >> numSources;
  if(numSources <= 0){
    cerr<<"can not get numSources " << endl;
    Loci::Abort();
  }
  sources.resize(numSources);
  
 
  
  for(int i = 0; i < numSources; i++){

    source_par s;
    vect3d p1, p2;
    inFile >> p1.x >> p1.y >>p1.z;
    inFile >> p2.x >> p2.y >> p2.z;
    inFile >> s.r0 >> s.s0 >> s.r1 >> s.s1 >> s.a;
    
    s.p1 = p1;
    s.p2 = p2;
    sources[i] = s;
  }
  inFile.close();  
}

/**
 * Computes the shortest distance from a point to a line segment.
 *
 * Degenerate segments are treated as a single endpoint.
 *
 * @param p   Query point.
 * @param p1  First segment endpoint.
 * @param p2  Second segment endpoint.
 * @return Euclidean distance from @p p to the segment with endpoints @p p1 and
 *         @p p2.
 */
double get_distance(const vect3d& p, const vect3d& p1, const vect3d& p2){
  if( norm(p1-p2) < NORMALIZE_ZERO_THRESHOLD) return norm(p - p1);
  double dotp1 = dot(p2-p1, p-p1);
  double dotp2 = dot(p1-p2, p-p2);
  if(dotp1 > 0 && dotp2 > 0) return norm(cross(p2-p1, p-p1))/norm(p2-p1);
  if(dotp1 <= 0) return norm(p-p1);
  if(dotp2 <= 0) return norm(p-p2);
  cerr<<"WARNING: get_distance reach dummy code" << endl;
  return norm(p-0.5*(p1+p2));
}

/**
 * Return the spacing requested by source s at point p. Let r be the distance
 * to the segment from p1 to p2.
 *
 * Use s0 when r <= r0 and the projection of p lies strictly between the
 * endpoints. Otherwise use s1 for r <= r1, or s1 * pow(r/r1, a) outside that
 * distance.
 */
double get_spacing(const vect3d& p, const source_par& s){
 
  double r = get_distance(p, s.p1, s.p2);
  double dotp1 = dot(s.p2-s.p1, p-s.p1);
  double dotp2 = dot(s.p1-s.p2, p-s.p2);
  if(dotp1 > 0 && dotp2 > 0 && r <= s.r0) return s.s0;
  if(r <= s.r1) return s.s1;
  return s.s1*pow(r/s.r1, s.a);
}

/**
 * Return the mean of nodes[i]->p. An empty nodes vector reports an error and
 * returns the zero vector.
 */
vect3d get_center(const vector<Node*>& nodes){
  vect3d center = vect3d(0.0, 0.0, 0.0);
  int num_nodes = nodes.size();
  if(num_nodes == 0){
    cerr << "ERROR: num of nodes of a cell is zero" << endl;
    return center;
  }
  for(int i = 0; i< num_nodes; i++)
    center = center + nodes[i]->p;
  center = (1.0/ num_nodes)*center;
  return center;
}

/**
 * Return the minimum source spacing sampled at the cell vertices and their
 * mean position. If ss is empty, return numeric_limits<double>::max().
 */
double get_min_spacing(const vector<Node*>& nodes, const vector<source_par>& ss){
  vect3d center = get_center(nodes);
  double spacing = std::numeric_limits<double>::max();
  for(unsigned int i = 0; i < ss.size(); i++){
    for(unsigned int np = 0; np < nodes.size(); np++){
      spacing = min(spacing, get_spacing(nodes[np]->p, ss[i]));
    }
    spacing = min(spacing,get_spacing(center, ss[i]));
  }
  return spacing;
}


/**
 * Return 1 if a source-spacing check requests refinement, otherwise 0.
 *
 * For a source with min_edge_len > r1/2, refine if either endpoint is inside
 * the cell bounding box or the distance from the mean vertex position to the
 * segment midpoint is less than min_edge_len + r1. Otherwise, refine if
 * min_edge_len is at least twice the minimum spacing sampled at the vertices
 * and their mean position.
 */
int tag_cell(const vector<Node*>& nodes, const vector<source_par>& sources, double min_edge_len){


  //get bounding box of the cell
  double minx, miny, minz, maxx, maxy, maxz;
  minx=miny=minz=std::numeric_limits<float>::max();
  maxx=maxy=maxz=std::numeric_limits<float>::min();
  for(unsigned int i = 0; i < nodes.size(); i++){
    minx = min(minx, (nodes[i]->p).x);
    miny = min(miny, (nodes[i]->p).y);
    minz = min(minz, (nodes[i]->p).z);
    maxx = max(maxx, (nodes[i]->p).x);
    maxy = max(maxy, (nodes[i]->p).y);
    maxz = max(maxz, (nodes[i]->p).z);
  }
  vect3d center = get_center(nodes); 
  
  for(unsigned int i = 0; i<sources.size(); i++){
    //if any feature is in the bounding box,
    //and min_edge_length is larger than 0.5r1,
    //then split the cell
    vect3d p1 = sources[i].p1;
    vect3d p2 = sources[i].p2;
    if((p1.x <= maxx && p1.x >= minx &&
        p1.y <= maxy && p1.y >= miny &&
        p1.z <= maxz && p1.z >= minz
        )||(p2.x <= maxx && p2.x >= minx &&
            p2.y <= maxy && p2.y >= miny &&
            p2.z <= maxz && p2.z >= minz)){
      
      if(min_edge_len>0.5*sources[i].r1)return 1;
    }
    //if anyany feature sphere intersect with the cell, 
    //and min_edge_length of the cell is larger than 0.5r1,
    //then split the cell
    
    
    if((norm(center-0.5*(p1+p2))< (min_edge_len+sources[i].r1)) &&
       (min_edge_len>0.5*sources[i].r1))return 1;
    
  }
      
  double min_spacing = get_min_spacing(nodes, sources);
  if(min_edge_len>=2*min_spacing)return 1;
  else return 0;
  
}
