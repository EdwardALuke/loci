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
#include "vogtools.h"
#include <map>
#include <unistd.h>

#include <list>
using std::list ;
#include <string>
using std::string ;
#include <vector>
using std::vector ;
using std::pair ;
using std::cout ;
using std::endl ;
using std::cerr ;
using std::ifstream ;
using std::ios ;

//using Loci::fact_db ;
using Loci::debugout ;

namespace Loci {

  extern void distributed_inverseMap(multiMap &result,
                                     vector<pair<Entity,Entity> > &input,
                                     entitySet input_image,
                                     entitySet input_preimage,
                                     const std::vector<entitySet> &init_ptn) ;
}
namespace VOG {
  
  //#define MEMDIAG
  
#ifdef MEMDIAG
  void *memtop =0;

  class beginexec {
  public:
    beginexec() {
      memtop = sbrk(0) ;
    }
  } ;

  beginexec hackit;
#endif
  void memSpace(string s) {
#ifdef MEMDIAG

    unsigned long memsize = (char *)sbrk(0)-(char *)memtop ;
    debugout << s << ": malloc = " << double(memsize)/(1024.*1024) << endl ;
    //#define MEMINFO
#ifdef MEMINFO
    struct mallinfo info = mallinfo() ;
    debugout << s << ": minfo, arena=" << info.arena
             << ", ordblks=" << info.ordblks
             << ", hblks="<< info.hblks
             << ", hblkhd="<< info.hblkhd
             << ", uordblks=" << info.uordblks
             << ", fordblks="<<info.fordblks
             << ", keepcost="<<info.keepcost << endl ;
#endif
    debugout.flush() ;
#endif
  }

  vector<BC_descriptor> readTags(string filename) {
    vector<BC_descriptor> bcs ;
    ifstream file(filename.c_str(),ios::in) ;
    if(file.fail()) {
      return bcs ;
    }
    char buf[1024] ;
    buf[1023] = '\0' ;
    if(file.peek() != '#') {
      while(file.peek() != EOF) {
        int id = -1;
        file >> id ;
        if(id == -1)
          break ;
        string name ;
        file >> name ;
        if(name == "")
          break ;
        string tmp = name ;
        size_t nsz = name.size() ;
        if(!(name[0] >= 'a' && name[0] <= 'z') &&
           !(name[0] >= 'A' && name[0] <= 'Z'))
          name[0] = '_' ;
        for(size_t i=1;i<nsz;++i) 
          if(!(name[i] >= 'a' && name[i] <= 'z') &&
             !(name[i] >= 'A' && name[i] <= 'Z') &&
             !(name[i] >= '0' && name[i] <= '9'))
            name[i] = '_' ;
        if(tmp != name) 
          cerr << "Renaming tag '" << tmp << "' to '" << name << "'!" << endl ;
        
        BC_descriptor BC ;
        BC.id = id ;
        BC.name = name ;
        BC.Trans = false ;
        file.getline(buf,1023) ;
        bcs.push_back(BC) ;
        while(file.peek() != EOF && (isspace(file.peek())))
          file.get() ;
        while(file.peek() == '#') { // eat comments
          file.getline(buf,1023) ;
          while(file.peek() != EOF && (isspace(file.peek())))
            file.get() ;
        }
      }
      return bcs ;
    }
    file.getline(buf,1023) ;
    if(strncmp(buf,"# Generated",11) == 0) {
      cout << "Tagfile: " << buf << endl ;
      // New tag file, dump first two lines
      file.getline(buf,1023) ;
      file.getline(buf,1023) ;
      int nsp = 0 ;
      for(int i=0;i<1023;++i) {
        if(buf[i] == '\0')
          break ;
        if(isspace(buf[i])) {
          nsp++ ;
          while(isspace(buf[i]))
            i++ ;
        }
        if(strncmp(&buf[i],"Trans",5) == 0)
          break ;
      }
          
      int n2trans = nsp-1 ;
      while(file.peek() != EOF) {
        int id = -1;
        file >> id ;
        if(id == -1)
          break ;
        string name ;
        file >> name ;
        if(name == "")
          break ;
        string tmp = name ;
        size_t nsz = name.size() ;
        if(!(name[0] >= 'a' && name[0] <= 'z') &&
           !(name[0] >= 'A' && name[0] <= 'Z'))
          name[0] = '_' ;
        for(size_t i=1;i<nsz;++i) 
          if(!(name[i] >= 'a' && name[i] <= 'z') &&
             !(name[i] >= 'A' && name[i] <= 'Z') &&
             !(name[i] >= '0' && name[i] <= '9'))
            name[i] = '_' ;
        if(tmp != name) 
          cerr << "Renaming tag '" << tmp << "' to '" << name << "'!" << endl ;
        
        bool trans ;
        for(int i=0;i<n2trans;++i)
          file >> trans ;
        BC_descriptor BC ;
        BC.id = id ;
        BC.name = name ;
        BC.Trans = trans ;
        file.getline(buf,1023) ;
        bcs.push_back(BC) ;
        while(file.peek() != EOF && (isspace(file.peek())))
          file.get() ;
        while(file.peek() == '#') { // eat comments
          file.getline(buf,1023) ;
          while(file.peek() != EOF && (isspace(file.peek())))
            file.get() ;
        }
      }
    } else {
      if(strncmp(buf,"#ID:",4) != 0) {
        cerr << "unable to interpret tag file format" << endl ;
        cerr << "first line: " << buf << endl ;
        return bcs ;
      }
      int nsp = 0 ;
      for(int i=0;i<120;++i) {
        if(buf[i] == '\0')
          break ;
        if(isspace(buf[i])) {
          nsp++ ;
          while(isspace(buf[i]))
            i++ ;
        }
        if(strncmp(&buf[i],"Trans",5) == 0)
          break ;
      }
      int n2trans = nsp-1 ;
      while(file.peek() != EOF) {
        while(file.peek() != EOF&&file.peek() != '#')
          file.get() ;

        char c ;
        c = file.get() ;
        if(c != '#') {
          cerr << "expected to get a hash while reading '" << filename << "'"
               << endl ;
          bcs.clear() ;
          return bcs ;
        }
        int id = 0 ;
        file >> id ;
        c = file.get() ;
        if(c != ':') {
          cerr << "expected to get a colon while reading '" << filename << "'"
               << endl ;
          bcs.clear() ;
          return bcs ;
        }
        string name ;
        file >> name ;
        string tmp = name ;
        size_t nsz = name.size() ;
        if(!(name[0] >= 'a' && name[0] <= 'z') &&
           !(name[0] >= 'A' && name[0] <= 'Z'))
          name[0] = '_' ;
        for(size_t i=1;i<nsz;++i) 
          if(!(name[i] >= 'a' && name[i] <= 'z') &&
             !(name[i] >= 'A' && name[i] <= 'Z') &&
             !(name[i] >= '0' && name[i] <= '9'))
            name[i] = '_' ;
        if(tmp != name) 
          cerr << "Renaming tag '" << tmp << "' to '" << name << "'!" << endl ;

        bool trans ;
        for(int i=0;i<n2trans;++i)
          file >> trans ;
        BC_descriptor BC ;
        BC.id = id ;
        BC.name = name ;
        BC.Trans = trans ;
        file.getline(buf,1024) ;
        bcs.push_back(BC) ;
        while(file.peek() != EOF && (isspace(file.peek())))
          file.get() ;
      }
    }
    return bcs ;
  }
  
  using Loci::MapRepP ;
  using Loci::storeRepP ;
  using Loci::MPI_processes ;
  using Loci::MPI_rank ;
  using Loci::create_intervalSet ;
  using Loci::UNIVERSE_MIN ;
  using Loci::UNIVERSE_MAX ;

  struct nodeSort {
    int key ;
    int orig_id ;
    vector3d<double> pos ;
  } ;

  inline bool operator<(const nodeSort &n1, const nodeSort &n2) {
    return n1.key < n2.key ;
  }

  vector<int> simplePartitionVec(int mn, int mx, int p) {
    vector<int> nums(p+1) ;
    int n = mx-mn+1 ;
    int dn = n/p ; // divisor
    int rn = n%p ; // remainder
    int start = mn ;
    nums[0] = start ;
    for(int i=0;i<p;++i) {
      start += dn+((i<rn)?1:0) ;
      nums[i+1] = start ;
    }
    FATAL(start != mx+1) ;
    return nums ;
  }

  vector<entitySet> getNodesDist(entitySet &nodes, store<vector3d<double> > &pos) {
    // First establish current distribution of entities across processors
    vector<entitySet> ptn(MPI_processes) ; // entity Partition
    // Get entity distributions
    nodes = pos.domain() ;
    entitySet allNodes = Loci::all_collect_entitySet(nodes) ;
    vector<int> nodesizes(MPI_processes) ;
    int size = nodes.size() ;
    MPI_Allgather(&size,1,MPI_INT,&nodesizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
    int cnt = allNodes.Min() ;
    for(int i=0;i<MPI_processes;++i) {
      ptn[i] = interval(cnt,cnt+nodesizes[i]-1) ;
      cnt += nodesizes[i] ;
    }
    nodes = allNodes ;
    return ptn ;
  }

  vector<entitySet> getFacesDist(entitySet &faces,multiMap &face2node) {
    // First establish current distribution of entities across processors
    vector<entitySet> ptn(MPI_processes) ; // entity Partition
    
    faces = face2node.domain() ;
    entitySet allFaces = Loci::all_collect_entitySet(faces) ;
    vector<int> facesizes(MPI_processes) ;
    int size = faces.size() ;
    MPI_Allgather(&size,1,MPI_INT,&facesizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
    int cnt = allFaces.Min() ;
    for(int i=0;i<MPI_processes;++i) {
      ptn[i] += interval(cnt,cnt+facesizes[i]-1) ;
      cnt += facesizes[i] ;
    }
    faces = allFaces ;
    return ptn ;
  }

  vector<entitySet> getCellsDist(entitySet &cells,Map &cl, Map &cr) {
    // First establish current distribution of entities across processors
    vector<entitySet> ptn(MPI_processes) ; // entity Partition

    entitySet tmp_cells = cl.image(cl.domain())+cr.image(cr.domain()) ;
    entitySet loc_geom_cells = tmp_cells & interval(0,Loci::UNIVERSE_MAX) ;
    entitySet geom_cells = Loci::all_collect_entitySet(loc_geom_cells) ;
    int mn = geom_cells.Min() ;
    int mx = geom_cells.Max() ;
    vector<int> pl = simplePartitionVec(mn,mx,MPI_processes) ;
    for(int i=0;i<MPI_processes;++i)
      ptn[i] += interval(pl[i],pl[i+1]-1) ;
    cells = geom_cells ;

    return ptn ;
  }

  vector<entitySet> getDist(entitySet &nodes, entitySet &faces,
                            entitySet &cells,
                            store<vector3d<double> > &pos,
                            Map &cl, Map &cr, multiMap &face2node) {
    // First establish current distribution of entities across processors
    vector<entitySet> ptn(MPI_processes) ; // entity Partition

    // Get entity distributions
    nodes = pos.domain() ;
    entitySet allNodes = Loci::all_collect_entitySet(nodes) ;
    vector<int> nodesizes(MPI_processes) ;
    int size = nodes.size() ;
    MPI_Allgather(&size,1,MPI_INT,&nodesizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
    int cnt = allNodes.Min() ;
    for(int i=0;i<MPI_processes;++i) {
      ptn[i] = interval(cnt,cnt+nodesizes[i]-1) ;
      cnt += nodesizes[i] ;
    }

    faces = face2node.domain() ;
    entitySet allFaces = Loci::all_collect_entitySet(faces) ;
    vector<int> facesizes(MPI_processes) ;
    size = faces.size() ;
    MPI_Allgather(&size,1,MPI_INT,&facesizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
    cnt = allFaces.Min() ;
    for(int i=0;i<MPI_processes;++i) {
      ptn[i] += interval(cnt,cnt+facesizes[i]-1) ;
      cnt += facesizes[i] ;
    }
    
    entitySet tmp_cells = cl.image(cl.domain())+cr.image(cr.domain()) ;
    entitySet loc_geom_cells = tmp_cells & interval(0,Loci::UNIVERSE_MAX) ;
    entitySet geom_cells = Loci::all_collect_entitySet(loc_geom_cells) ;
    int mn = geom_cells.Min() ;
    int mx = geom_cells.Max() ;
    vector<int> pl = simplePartitionVec(mn,mx,MPI_processes) ;
    for(int i=0;i<MPI_processes;++i)
      ptn[i] += interval(pl[i],pl[i+1]-1) ;
    nodes = allNodes ;
    faces = allFaces ;
    cells = geom_cells ;
    return ptn ;
  }

  // Establish geometrically consistent face orientation
  void orientFaces(store<vector3d<double> > &pos,
                   Map &cl, Map &cr, multiMap &face2node) {

    entitySet nodes, faces,cells ;
    vector<entitySet> init_ptn = getDist(nodes,faces,cells,
                                         pos,cl,cr,face2node);
    
    store<vector3d<double> > fpos ;
    store<vector3d<double> > area ;

    entitySet loc_faces = faces & init_ptn[MPI_rank] ;
    entitySet geom_cells = cells ;
    entitySet negs = interval(UNIVERSE_MIN,-1) ;
    entitySet boundary_faces = cr.preimage(negs).first ;
    entitySet interior_faces = loc_faces - boundary_faces ;

    entitySet myNodes = nodes & init_ptn[MPI_rank] ;
    entitySet total_dom =
      Loci::MapRepP(face2node.Rep())->image(loc_faces) + myNodes ;


    dstore<vector3d<double> > tmp_pos ;
    FORALL(pos.domain(), pi) {
      tmp_pos[pi] = pos[pi] ;
    } ENDFORALL ;
    
    if(MPI_processes > 1) {
      Loci::storeRepP sp = tmp_pos.Rep() ;
      Loci::fill_clone(sp, total_dom, init_ptn) ;
    }
    
    entitySet face_dom = face2node.domain() ;
    fpos.allocate(face_dom) ;
    area.allocate(face_dom) ;

    FORALL(face_dom,fc) {
      int nnodes = face2node.end(fc) - face2node.begin(fc) ;
      vector3d<double> fp(0,0,0) ;
      double w = 0 ;
      for(int i=0;i<nnodes;++i) {
        vector3d<double> p1 = (tmp_pos[face2node[fc][i]]) ;
        vector3d<double> p2 = (tmp_pos[face2node[fc][(i+1)%nnodes]]) ;

        double len = norm(p1-p2) ;

        fp += len*(p1+p2) ;
        w += len ;
      }
      fpos[fc] = fp/(2.*w) ;
      vector3d<double> a(0,0,0) ;
      for(int i=0;i<nnodes;++i) {
        vector3d<double> p1 = (tmp_pos[face2node[fc][i]]) ;
        vector3d<double> p2 = (tmp_pos[face2node[fc][(i+1)%nnodes]]) ;
        a += cross(p1-fpos[fc],p2-fpos[fc]) ;
      }
      area[fc] = .5*a ;
    } ENDFORALL ;


    dstore<vector3d<double> > cpos ;
    dstore<double> cnum ;

    entitySet tmp_cells =  cl.image(loc_faces) | cr.image(interior_faces) ;
    // Add cells owned by this processor!
    tmp_cells += cells & init_ptn[MPI_rank] ; 
    cpos.allocate(tmp_cells) ;
    cnum.allocate(tmp_cells) ;
    FORALL(tmp_cells,cc) {
      cpos[cc] = vector3d<double>(0,0,0) ;
      cnum[cc] = 0 ;
    } ENDFORALL ;
    FORALL(loc_faces,fc) {
      double A = norm(area[fc]) ;
      cpos[cl[fc]] += A*fpos[fc] ;
      cnum[cl[fc]] += A ;
    } ENDFORALL ;
    FORALL(interior_faces,fc) {
      double A = norm(area[fc]) ;
      cpos[cr[fc]] += A*fpos[fc] ;
      cnum[cr[fc]] += A ;
    } ENDFORALL ;
    Loci::storeRepP cp_sp = cpos.Rep() ;
    Loci::storeRepP cn_sp = cnum.Rep() ;
    entitySet clone_cells = tmp_cells - (cells&init_ptn[MPI_rank]) ;
    std::vector<Loci::storeRepP> v_cpos = Loci::send_global_clone_non(cp_sp, clone_cells, init_ptn) ;
    std::vector<Loci::storeRepP> v_cnum = Loci::send_global_clone_non(cn_sp, clone_cells, init_ptn) ;
    for(int i = 0; i < Loci::MPI_processes; ++i) {
      entitySet dom = v_cpos[i]->domain() & cpos.domain() ;
      dstore<vector3d<double> > tmp_cpos(v_cpos[i]) ;
      dstore<double> tmp_cnum(v_cnum[i]) ;
      FORALL(dom, di) {
	cpos[di] += tmp_cpos[di] ;
	cnum[di] += tmp_cnum[di] ;
      } ENDFORALL ;
    }
    Loci::fill_clone(cp_sp, clone_cells, init_ptn) ;
    Loci::fill_clone(cn_sp, clone_cells, init_ptn) ;   
    FORALL(tmp_cells,cc) {
      cpos[cc] = cpos[cc]/cnum[cc] ;
    } ENDFORALL ;

    vector<int> broken_faces ;

    FORALL(interior_faces,fc) {
      vector3d<double> dv = cpos[cr[fc]]-cpos[cl[fc]] ; 
      vector3d<double> dv2 = fpos[fc]-cpos[cl[fc]] ; 
      vector3d<double> dv3 = cpos[cr[fc]]-fpos[fc] ; 

      int t1 = (dot(area[fc],dv) <0.0)?1:0 ;
      int t2 = (dot(area[fc],dv2) <0.0)?1:0 ;
      int t3 = (dot(area[fc],dv3) <0.0)?1:0 ;
      int test = t1+t2+t3 ;
      if(test != 3 && test != 0) {
        debugout << "problem with face located at " << fpos[fc]
                 << endl ;
        broken_faces.push_back(fc) ;
      }

      
      else if(t1 == 1) { // Face oriented incorrectly
	int i = 0 ;
	int j = face2node.end(fc) - face2node.begin(fc) -1 ;
	while(i < j) {
          std::swap(face2node[fc][i],face2node[fc][j]) ;
	  i++ ;
	  j-- ;
	} 
      }
    } ENDFORALL ;


    FORALL(boundary_faces,fc) {
      const vector3d<double> center = fpos[fc] ;

      vector3d<double> ccenter ;
      ccenter = cpos[cl[fc]] ;
      vector3d<double> dv = center-ccenter ;
      if(dot(area[fc],dv) < 0.0) {
	int i = 0 ;
	int j = face2node.end(fc) - face2node.begin(fc) -1 ;
	while(i < j) {
          std::swap(face2node[fc][i],face2node[fc][j]) ;
	  i++ ;
	  j-- ;
	} 
      }
      
    } ENDFORALL ;

    int rsize = 0 ;
    int size = broken_faces.size() ;
    
    MPI_Allreduce(&size,&rsize,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD) ;
    if(rsize!=0) {
      if(MPI_rank == 0) {
        cerr << "Bad Grid: Non-Convex Cell (centroid outside cell bounds)" << endl ;
      }
      //      Loci::Abort() ;
    }
  }


  void colorMatrix(store<vector3d<double> > &pos,
                   Map &cl, Map &cr, multiMap &face2node) {
    
    entitySet nodes, faces,cells ;
    //    vector<entitySet> ptn = getDist(nodes,faces,cells,
    //                                    pos,cl,cr,face2node);
    //    vector<entitySet> nptn = getNodesDist(nodes,pos) ;
    vector<entitySet> fptn = getFacesDist(faces,face2node) ;
    vector<entitySet> cptn = getCellsDist(cells,cl,cr) ;
    entitySet loc_faces = faces & fptn[MPI_rank] ;
    entitySet geom_cells = cells & cptn[MPI_rank] ;
    entitySet negs = interval(UNIVERSE_MIN,-1) ;
    entitySet boundary_faces = cr.preimage(negs).first ;
    entitySet interior_faces = loc_faces - boundary_faces ;

    using std::pair ;
    vector<pair<Entity,Entity> > cellmap(interior_faces.size()*2) ;
    int cnt = 0 ;
    FORALL(interior_faces,fc) {
      cellmap[cnt++] = pair<Entity,Entity>(cl[fc],cr[fc]) ;
      cellmap[cnt++] = pair<Entity,Entity>(cr[fc],cl[fc]) ;
    } ENDFORALL ;
    multiMap c2c ;
    Loci::distributed_inverseMap(c2c,cellmap,cells,cells,cptn) ;

    store<int> ctmp ;
    ctmp.allocate(geom_cells) ;
    FORALL(geom_cells,cc) {
      ctmp[cc] = -1 ;
    } ENDFORALL ;

    int ncells_local = geom_cells.size() ;
    int ncells = 0 ; 
    MPI_Allreduce(&ncells_local,&ncells, 1, MPI_INT, MPI_MAX,MPI_COMM_WORLD) ;
    int col = ncells*Loci::MPI_rank ;
    
    vector<int> visited ;
    entitySet left_out = geom_cells ;
    int lo_p = geom_cells.Min() ;
    while(left_out != EMPTY) {
      vector<int> work ;
      work.push_back(left_out.Min()) ;
      while(work.size() != 0) {
	vector<int> working ;
	for(size_t i=0;i<work.size();++i) {
          int cc = work[i] ;
          if(ctmp[cc] == -1) {
            ctmp[cc] = col++ ;
            visited.push_back(cc) ;
            for(const int *pi = c2c.begin(cc);pi!=c2c.end(cc);++pi)
              if(geom_cells.inSet(*pi) && ctmp[*pi] == -1) {
                working.push_back(*pi) ;
              }
          }
	}
        work.swap(working) ;
      }
      left_out = EMPTY ;
      entitySet candidates = geom_cells & interval(lo_p,UNIVERSE_MAX) ;
      FORALL(candidates,cc) {
        if(ctmp[cc] == -1) {
          left_out += cc ;
          break ;
        }
      } ENDFORALL ;
      if(left_out != EMPTY)
        lo_p = left_out.Min() ;
    }

    dstore<int> color ;
    FORALL(geom_cells,cc) {
      color[cc] = ctmp[cc];
    } ENDFORALL ;

    entitySet clone_cells = cl.image(interior_faces)
      + cr.image(interior_faces) ;
    clone_cells -= geom_cells ;
    Loci::storeRepP cp_sp = color.Rep() ;
    Loci::fill_clone(cp_sp, clone_cells, cptn) ;

    FORALL(interior_faces,fc) {
      int color_l = color[cl[fc]] ;
      int color_r = color[cr[fc]] ;
      if(color_l == color_r) 
        cerr << "color equal == " << color_l << endl ;
      if(color_l == -1 || color_r == -1)
        cerr << "matrix coloring internal error" << endl ;
                                                              
      if(color_l > color_r) {
        // change face orientation to match matrix coloring
        std::swap(cl[fc],cr[fc]) ;
        int i = 0 ;
        int j = face2node[fc].size() - 1;
       	while(i < j) {
          std::swap(face2node[fc][i],face2node[fc][j]) ;
          i++ ;
          j-- ;
        } 
      }
    } ENDFORALL ;
    
  }


  void getCellCenters(store<vector3d<double> > &cellcenter,
                      store<vector3d<double> > &pos,
                      Map &cl, Map &cr, multiMap &face2node,
                      vector<entitySet> &nptn,
		      vector<entitySet> &cptn,
		      vector<entitySet> &fptn) {
    entitySet faces = face2node.domain() ;
    entitySet tmp_cells = cl.image(cl.domain())+cr.image(cr.domain()) ;
    entitySet loc_geom_cells = tmp_cells & interval(0,Loci::UNIVERSE_MAX) ;
    entitySet geom_cells = Loci::all_collect_entitySet(loc_geom_cells) ;

    // Compute cell centers
    entitySet total_dom =
      Loci::MapRepP(face2node.Rep())->image(faces) + pos.domain() ;

    dstore<vector3d<double> > tmp_pos ;
    FORALL(pos.domain(), pi) {
      tmp_pos[pi] = pos[pi] ;
    } ENDFORALL ;
    
    if(MPI_processes > 1) {
      Loci::storeRepP sp = tmp_pos.Rep() ;
      Loci::fill_clone(sp, total_dom, nptn) ;
    }
    
    store<vector3d<double> > fpos, area ;
    fpos.allocate(faces) ;
    area.allocate(faces) ;

    FORALL(faces,fc) {
      int nnodes = face2node[fc].size() ;
      vector3d<double> fp(0,0,0) ;
      double w = 0 ;
      for(int i=0;i<nnodes;++i) {
        vector3d<double> p1 = (tmp_pos[face2node[fc][i]]) ;
        vector3d<double> p2 = (tmp_pos[face2node[fc][(i+1)%nnodes]]) ;

        double len = norm(p1-p2) ;

        fp += len*(p1+p2) ;
        w += len ;
      }
      fpos[fc] = fp/(2.*w) ;
      vector3d<double> a(0,0,0) ;
      for(int i=0;i<nnodes;++i) {
        vector3d<double> p1 = (tmp_pos[face2node[fc][i]]) ;
        vector3d<double> p2 = (tmp_pos[face2node[fc][(i+1)%nnodes]]) ;
        a += cross(p1-fpos[fc],p2-fpos[fc]) ;
      }
      area[fc] = .5*a ;
    } ENDFORALL ;

    dstore<vector3d<double> > cpos ;
    dstore<double> cnum ;

    // Add cells owned by this processor!
    tmp_cells &= interval(0,UNIVERSE_MAX) ;
    tmp_cells += geom_cells & cptn[MPI_rank] ; 
    cpos.allocate(tmp_cells) ;
    cnum.allocate(tmp_cells) ;
    FORALL(tmp_cells,cc) {
      cpos[cc] = vector3d<double>(0,0,0) ;
      cnum[cc] = 0 ;
    } ENDFORALL ;
    FORALL(faces,fc) {
      double A = norm(area[fc]) ;
      cpos[cl[fc]] += A*fpos[fc] ;
      cnum[cl[fc]] += A ;
    } ENDFORALL ;
    FORALL(faces,fc) {
      if(cr[fc] >= 0) {
        double A = norm(area[fc]) ;
        cpos[cr[fc]] += A*fpos[fc] ;
        cnum[cr[fc]] += A ;
      }
    } ENDFORALL ;

    entitySet clone_cells = tmp_cells - cptn[MPI_rank] ;

    Loci::storeRepP cp_sp = cpos.Rep() ;
    Loci::storeRepP cn_sp = cnum.Rep() ;
    std::vector<Loci::storeRepP> v_cpos =
      Loci::send_global_clone_non(cp_sp, clone_cells, cptn) ;
    std::vector<Loci::storeRepP> v_cnum =
      Loci::send_global_clone_non(cn_sp, clone_cells, cptn) ;

    for(int i = 0; i < Loci::MPI_processes; ++i) {
      entitySet dom = v_cpos[i]->domain() & cpos.domain() ;
      dstore<vector3d<double> > tmp_cpos(v_cpos[i]) ;
      dstore<double> tmp_cnum(v_cnum[i]) ;
      FORALL(dom, di) {
	cpos[di] += tmp_cpos[di] ;
	cnum[di] += tmp_cnum[di] ;
      } ENDFORALL ;
    }
    Loci::fill_clone(cp_sp, clone_cells, cptn) ;
    Loci::fill_clone(cn_sp, clone_cells, cptn) ;   

    tmp_cells = geom_cells & cptn[MPI_rank] ;
    cellcenter.allocate(tmp_cells) ;
    FORALL(tmp_cells,cc) {
      cellcenter[cc] = cpos[cc]/cnum[cc] ;
    } ENDFORALL ;
  }


  // Hilbert key
  typedef Loci::Array<unsigned int,3> Hcode ;

  // Integer coordinate
  typedef Loci::Array<unsigned int,3> Point ;

  struct Key {
    Hcode key ; 
    int key_id ;
  } ;

  inline bool operator<(const Key &k1, const Key &k2) {
    return (
            (k1.key[2]<k2.key[2]) ||
            (k1.key[2]==k2.key[2]&&k1.key[1]<k2.key[1]) ||
            (k1.key[2]==k2.key[2]&&k1.key[1]==k2.key[1]&&k1.key[0]<k2.key[0])
            ) ;
  }

  // Morton Space filling curve encoding
  Hcode M_encode(Point p) {
    int bit = 0 ;
    Hcode m ;
    m[0] = 0 ;
    m[1] = 0 ;
    m[2] = 0 ;
    for(int i=0;i<32;++i) 
      for(int j=0;j<3;++j) {
        if((p[j] & (1 << i)) != 0) {
          int indx = bit>>5 ;
          int b = bit&31 ;
          m[indx] |= (1 << b) ;
        }
        bit++ ;
      }
    return m ;
  }
  // Mask for 3-Dimensions
  const unsigned int g_mask[] = {4,2,1};

  // given the coordinates of a point find the sequence number of the point
  // on the Hilbert Curve
  Hcode H_encode(Point p) {
    const int DIM = 3 ;
    const int WORDBITS = 32 ;
    const int NUMBITS = 32 ;
    unsigned int mask = (unsigned long)1 << (WORDBITS - 1) ;
    unsigned int element, temp1, temp2, A, W = 0, S, tS, T, tT, J, P = 0, xJ;
    
    Hcode	h;
    h[0] = 0;
    h[1] = 0;
    h[2] = 0;
    
    int	i = NUMBITS * DIM - DIM, j;
    
    for (j = A = 0; j < DIM; j++)
      if (p[j] & mask)
        A |= g_mask[j];
    
    S = tS = A;
    
    P |= S & g_mask[0];
    for (j = 1; j < DIM; j++)
      if( (S & g_mask[j]) ^ ((P >> 1) & g_mask[j]))
        P |= g_mask[j];
    
    /* add in DIM bits to hcode */
    element = i / WORDBITS;
    if (i % WORDBITS > WORDBITS - DIM) {
        h[element] |= P << i % WORDBITS;
	h[element + 1] |= P >> (WORDBITS - i % WORDBITS) ;
    } else
      h[element] |= P << (i - element * WORDBITS) ;

    J = DIM;
    for (j = 1; j < DIM; j++)
      if ((P >> j & 1) == (P & 1))
        continue;
      else
        break;
    if (j != DIM)
      J -= j;
    xJ = J - 1;
    
    if (P < 3)
      T = 0;
    else
      if (P % 2)
        T = (P - 1) ^ (P - 1) / 2;
      else
        T = (P - 2) ^ (P - 2) / 2;
    tT = T;
    
    for (i -= DIM, mask >>= 1; i >=0; i -= DIM, mask >>= 1) {
      for (j = A = 0; j < DIM; j++)
        if (p[j] & mask)
          A |= g_mask[j];
      
      W ^= tT;
      tS = A ^ W;
      if (xJ % DIM != 0) {
        temp1 = tS << xJ % DIM;
        temp2 = tS >> (DIM - xJ % DIM) ;
        S = temp1 | temp2;
        S &= ((unsigned int)1 << DIM) - 1;
      } else
        S = tS;

      P = S & g_mask[0];
      for (j = 1; j < DIM; j++)
        if( (S & g_mask[j]) ^ ((P >> 1) & g_mask[j]))
          P |= g_mask[j];
      
      /* add in DIM bits to hcode */
      element = i / WORDBITS;
      if (i % WORDBITS > WORDBITS - DIM) {
        h[element] |= P << i % WORDBITS;
        h[element + 1] |= P >> (WORDBITS - i % WORDBITS) ;
      } else
        h[element] |= P << (i - element * WORDBITS) ;

      if (i > 0) {
        if (P < 3)
          T = 0;
        else
          if (P % 2)
            T = (P - 1) ^ (P - 1) / 2;
          else
            T = (P - 2) ^ (P - 2) / 2;
        
        if (xJ % DIM != 0) {
          temp1 = T >> xJ % DIM;
          temp2 = T << (DIM - xJ % DIM) ;
          tT = temp1 | temp2;
          tT &= ((unsigned int)1 << DIM) - 1;
        } else
          tT = T;
        
        J = DIM;
        for (j = 1; j < DIM; j++)
          if ((P >> j & 1) == (P & 1))
            continue;
          else
            break;
        if (j != DIM)
          J -= j;
        
        xJ += J - 1;

      }
    }
    return h;
  }


  // Optimize indicies of mesh to increase locality
  void optimizeMesh(store<vector3d<double> > &pos,
                    Map &cl, Map &cr, multiMap &face2node) {
    
    memSpace("start optimizeMesh") ;
    // First establish current distribution of entities across processors
    vector<entitySet> nptn(MPI_processes) ; // entity Partition
    vector<entitySet> cptn(MPI_processes) ; // entity Partition
    vector<entitySet> fptn(MPI_processes) ; // entity Partition

    // Get entity distributions
    entitySet nodes = pos.domain() ;
    entitySet allNodes = Loci::all_collect_entitySet(nodes) ;
    vector<int> nodesizes(MPI_processes) ;
    int size = nodes.size() ;
    MPI_Allgather(&size,1,MPI_INT,&nodesizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
    int cnt = allNodes.Min() ;
    for(int i=0;i<MPI_processes;++i) {
      nptn[i] = interval(cnt,cnt+nodesizes[i]-1) ;
      cnt += nodesizes[i] ;
    }
      
    entitySet faces = face2node.domain() ;
    entitySet allFaces = Loci::all_collect_entitySet(faces) ;
    vector<int> facesizes(MPI_processes) ;
    size = faces.size() ;
    MPI_Allgather(&size,1,MPI_INT,&facesizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
    cnt = allFaces.Min() ;
    for(int i=0;i<MPI_processes;++i) {
      fptn[i] += interval(cnt,cnt+facesizes[i]-1) ;
      cnt += facesizes[i] ;
    }
    
    entitySet tmp_cells = cl.image(cl.domain())+cr.image(cr.domain()) ;
    entitySet loc_geom_cells = tmp_cells & interval(0,Loci::UNIVERSE_MAX) ;
    entitySet geom_cells = Loci::all_collect_entitySet(loc_geom_cells) ;
    int mn = geom_cells.Min() ;
    int mx = geom_cells.Max() ;
    vector<int> pl = simplePartitionVec(mn,mx,MPI_processes) ;
    for(int i=0;i<MPI_processes;++i)
      cptn[i] += interval(pl[i],pl[i+1]-1) ;


    memSpace("collect partition Info") ;
    // Compute distribution of cells based on space filling curve
    store<vector3d<double> > cellcenter ;
    getCellCenters(cellcenter, pos, cl,  cr, face2node, nptn, cptn, fptn) ;
    vector3d<double> maxVec, minVec,tmaxVec, tminVec ;
    
    memSpace("get cell centers") ;
    loc_geom_cells = geom_cells & cptn[MPI_rank] ;

    
    maxVec = cellcenter[loc_geom_cells.Min()] ;
    minVec = maxVec ;
    FORALL(loc_geom_cells,cc) {
      maxVec = vector3d<double>(max(maxVec.x,cellcenter[cc].x),
                                max(maxVec.y,cellcenter[cc].y),
                                max(maxVec.z,cellcenter[cc].z)) ;
      minVec = vector3d<double>(min(minVec.x,cellcenter[cc].x),
                                min(minVec.y,cellcenter[cc].y),
                                min(minVec.z,cellcenter[cc].z)) ;
    } ENDFORALL ;
    
    MPI_Allreduce(&maxVec.x,&tmaxVec.x,3,MPI_DOUBLE,
                  MPI_MAX,MPI_COMM_WORLD) ;
    MPI_Allreduce(&minVec.x,&tminVec.x,3,MPI_DOUBLE,
                  MPI_MIN,MPI_COMM_WORLD) ;

    maxVec = tmaxVec ;
    minVec = tminVec ;
    vector<Key> keyList(loc_geom_cells.size()) ;
    cnt = 0 ;

    vector3d<double> s = maxVec-minVec ;
    double scale = 4e9/max(s.x,max(s.y,s.z)) ;

    FORALL(loc_geom_cells,cc) {
      Point p ;
      vector3d<double> pbase = scale*(cellcenter[cc]-minVec) ;
      p[0] = (unsigned int)(pbase.x) ;
      p[1] = (unsigned int)(pbase.y) ;
      p[2] = (unsigned int)(pbase.z) ;
      keyList[cnt].key = H_encode(p) ;
      keyList[cnt++].key_id = cc ;
    } ENDFORALL ;

    Loci::parSampleSort(keyList,MPI_COMM_WORLD) ;
    
    memSpace("sorted Keys") ;
    vector<int> keysizes(MPI_processes) ;
    size = keyList.size() ;
    MPI_Allgather(&size,1,MPI_INT,&keysizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
    for(int i=1;i<MPI_processes;++i)
      keysizes[i] += keysizes[i-1] ;
    cnt = geom_cells.Min() ;
    if(MPI_rank > 0)
      cnt += keysizes[MPI_rank-1] ;

    vector<pair<Entity, Entity> > keypair(keyList.size()) ;
    for(size_t i=0;i<keyList.size();++i) {
      keypair[i] = pair<Entity,Entity>(cnt++,keyList[i].key_id) ;
    }
    multiMap mapping ;

    Loci::distributed_inverseMap(mapping,keypair,geom_cells,geom_cells,cptn) ;
    memSpace("inverseMap") ;
    // renumber cells 
    dMap cell2cell ;
    FORALL(loc_geom_cells,cc) {
      if(mapping[cc].size() != 1)
        cerr << "problem generating cell2cell map" << endl ;
      cell2cell[cc] = mapping[cc][0] ;
    } ENDFORALL ;
                               
    entitySet cimage = tmp_cells & interval(0,Loci::UNIVERSE_MAX) ;

    cell2cell.setRep(MapRepP(cell2cell.Rep())->expand(cimage,cptn)) ;
    
    memSpace("expand cell2cell map") ;
    FORALL(faces,fc) {
      cl[fc] = cell2cell[cl[fc]] ;
      if(cr[fc] >= 0)
        cr[fc] = cell2cell[cr[fc]] ;
    } ENDFORALL ;

    // Now renumber faces to match cell numbering
    multiMap nf2n ;
    Map ncl,ncr ;
    store<int> count ;
    entitySet newFaces ;
    const int p = MPI_processes ;
    if(p == 1) {
      newFaces = face2node.domain() ;
      nf2n.setRep(face2node.Rep()) ;
      ncl.setRep(cl.Rep()) ;
      ncr.setRep(cr.Rep()) ;
      count.allocate(newFaces) ;
    } else {
      vector<entitySet> send_sets(p) ;
      for(int i=0;i<MPI_processes;++i) 
	send_sets[i] = cl.preimage(cptn[i]&geom_cells).first ;

      vector<int> scounts(p,0) ;
      for(size_t i=0;i<send_sets.size();++i) {
	int nod_tot = 0 ;
	FORALL(send_sets[i],fc) {
	  nod_tot += face2node[fc].size() ;
	} ENDFORALL ;
	scounts[i] = 1 + send_sets[i].size()*3+nod_tot ;
      }
      vector<int> sdispls(p) ;
      sdispls[0] = 0 ;
      for(int i=1;i<p;++i)
	sdispls[i] = sdispls[i-1]+scounts[i-1] ;
      
      int send_size = sdispls[p-1]+scounts[p-1] ;
      vector<int> sbuffer(send_size) ;
      // Fill in send buffers
      for(size_t i=0;i<send_sets.size();++i) {
	int k = sdispls[i] ;
	sbuffer[k++] = send_sets[i].size() ;
	FORALL(send_sets[i],fc) {
	  sbuffer[k++] = cl[fc] ;
	  sbuffer[k++] = cr[fc] ;
	  sbuffer[k++] = face2node[fc].size() ;
	} ENDFORALL ;
	FORALL(send_sets[i],fc) {
	  for(int j=0;j<face2node[fc].size();++j)
	    sbuffer[k++] = face2node[fc][j] ;
	} ENDFORALL ;
      }
      
      vector<int> rcounts(p) ;
      MPI_Alltoall(&scounts[0],1,MPI_INT,&rcounts[0],1,MPI_INT,MPI_COMM_WORLD) ;
      
      vector<int> rdispls(p) ;
      rdispls[0] = 0 ;
      for(int i=1;i<p;++i) {
	rdispls[i] = rdispls[i-1]+rcounts[i-1] ;
      }
      int recv_size = rdispls[p-1]+rcounts[p-1] ;
      vector<int> rbuffer(recv_size) ;
      
      MPI_Alltoallv(&sbuffer[0],&scounts[0],&sdispls[0],MPI_INT,
		    &rbuffer[0],&rcounts[0],&rdispls[0],MPI_INT,
		    MPI_COMM_WORLD) ;
      
      int num_recv_faces = 0 ;
      for(int i=0;i<p;++i) {
	num_recv_faces += rbuffer[rdispls[i]] ;
      }
      vector<int> face_dist(p) ;
      MPI_Allgather(&num_recv_faces,1,MPI_INT,&face_dist[0],1,
		    MPI_INT,MPI_COMM_WORLD) ;
      int face_off = allFaces.Min() ;
      for(int i=0;i<MPI_rank-1;++i)
	face_off += face_dist[i] ;
      newFaces = interval(face_off,face_off+num_recv_faces-1) ;
      ncl.allocate(newFaces) ;
      ncr.allocate(newFaces) ;
      count.allocate(newFaces) ;
    
      memSpace("allocate new faces") ;
      int off = face_off ;
      for(int i=0;i<p;++i) {
	int k = rdispls[i] ;
	int nf = rbuffer[k++] ;
	for(int j=0;j<nf;++j) {
	  int fc = off+j ;
	  ncl[fc] = rbuffer[k++] ;
	  ncr[fc] = rbuffer[k++] ;
	  count[fc] = rbuffer[k++] ;
	}
	off += nf ;
      }
      nf2n.allocate(count) ;
      off = face_off ;
      for(int i=0;i<p;++i) {
	int k = rdispls[i] ;
	int nf = rbuffer[k++] ;
	k += nf*3 ;
	for(int j=0;j<nf;++j) {
	  int fc = off+j ;
	  for(int c=0;c<count[fc];++c)
	    nf2n[fc][c] = rbuffer[k++] ;
	}
	off += nf ;
      }
    }

    vector<pair<pair<int,int>,int> > sort_list(newFaces.size()) ;
    cnt = 0 ;
    FORALL(newFaces,fc) {
      sort_list[cnt].first.first = ncl[fc] ;
      sort_list[cnt].first.second = ncr[fc] ;
      sort_list[cnt].second = fc ;
      cnt++ ;
    } ENDFORALL ;
    std::sort(sort_list.begin(),sort_list.end()) ;
    Map scl,scr ;
    scl.allocate(newFaces) ;
    scr.allocate(newFaces) ;
    cnt = 0 ;
    FORALL(newFaces,fc) {
      scl[fc] = ncl[sort_list[cnt].second] ;
      scr[fc] = ncr[sort_list[cnt].second] ;
      count[fc] = nf2n[sort_list[cnt].second].size() ;
      cnt++ ;
    } ENDFORALL ;

    multiMap sf2n ;
    sf2n.allocate(count) ;
    cnt = 0 ;
    FORALL(newFaces,fc) {
      for(int j=0;j<count[fc];++j)
        sf2n[fc][j] = nf2n[sort_list[cnt].second][j] ;
      cnt++ ;
    } ENDFORALL ; 

    cl = scl.Rep() ;
    cr = scr.Rep() ;
    face2node = sf2n.Rep() ;
    memSpace("created new face ordering") ;
    // Now order nodes to match faces
    entitySet loc_faces = face2node.domain() ;
    entitySet node_access = MapRepP(face2node.Rep())->image(loc_faces)+pos.domain() ;
    dstore<int> node_key ;
    FORALL(node_access,nd) {
      node_key[nd] = std::numeric_limits<int>::max() ;
    } ENDFORALL ;

    FORALL(loc_faces,fc) {
      for(int i=0;i<face2node[fc].size();++i) {
        int nd = face2node[fc][i] ;
        node_key[nd] = min(node_key[nd],fc) ;
      }
    } ENDFORALL ;
    entitySet clone_nodes = node_access - pos.domain() ;
    Loci::storeRepP nk_sp = node_key.Rep() ;
    std::vector<Loci::storeRepP> v_nk =
      Loci::send_global_clone_non(nk_sp,clone_nodes,nptn) ;
    for(int i=0;i<MPI_processes;++i) {
      entitySet dom = v_nk[i]->domain() & pos.domain() ;
      dstore<int> tmp_nk(v_nk[i]) ;
      FORALL(dom,di) {
        node_key[di] = min(node_key[di],tmp_nk[di]) ;
      } ENDFORALL ;
    }

    vector<nodeSort> node_data(nodes.size()) ;
    cnt = 0 ;
    FORALL(nodes,nd) {
      node_data[cnt].key = node_key[nd] ;
      node_data[cnt].orig_id = nd ;
      node_data[cnt].pos = pos[nd] ;
      cnt++ ;
    } ENDFORALL ;

    Loci::parSampleSort(node_data,MPI_COMM_WORLD) ;

    memSpace("ordering Nodes") ;

    size = node_data.size() ;
    MPI_Allgather(&size,1,MPI_INT,&keysizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
    for(int i=1;i<MPI_processes;++i)
      keysizes[i] += keysizes[i-1] ;
    cnt = allNodes.Min() ;
    if(MPI_rank > 0)
      cnt += keysizes[MPI_rank-1] ;


    int base = cnt ;
    entitySet newnodes = interval(base,base+node_data.size()-1) ;
    store<vector3d<double> > npos ;
    npos.allocate(newnodes) ;

    vector<pair<Entity, Entity> > nodepair(node_data.size()) ;
    for(size_t i=0;i<node_data.size();++i) {
      nodepair[i] = pair<Entity,Entity>(cnt++,node_data[i].orig_id) ;
      npos[base+i] = node_data[i].pos ;
    }

    pos = npos.Rep() ;
    
    multiMap nmapping ;

    Loci::distributed_inverseMap(nmapping,nodepair,allNodes,allNodes,nptn) ;
    memSpace("inverseMap nmapping") ;
    // renumber nodes
    dMap node2node ;
    FORALL(nodes,nd) {
      if(nmapping[nd].size() != 1)
        cerr << "problem generating node2node map" << endl ;
      node2node[nd] = nmapping[nd][0] ;
    } ENDFORALL ;

    node2node.setRep(MapRepP(node2node.Rep())->expand(node_access,nptn)) ;
    memSpace("expanding node2node") ;
    FORALL(face2node.domain(),fc) {
      for(int i=0;i<face2node[fc].size();++i)
        face2node[fc][i] = node2node[face2node[fc][i]] ;
    } ENDFORALL ;
    
  }
}
