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
#ifndef MPI_CONTAINERIO_H
#define MPI_CONTAINERIO_H
#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <dist_internal.h>
#include <distribute_container.h>
#include <distribute.h>
#include <distribute_io.h>
#include <mpi.h>
#include <execute.h>
#include <multiStoreIO.h>
#include <mpi_multiStoreIO.h>

namespace Loci {
  /*
    mpi version to read/write stores. parallel io only
    Notice: each variable is read in from /write into a separate file
    argument xfer_type is used to control collective or independent data transfer
    global variable PHDF5_MPI_Info is used.
  */

  typedef struct {
    int ordered;
    int vec_size;
    int dom_size; //not entitySet.size(), but enritySet.num_intervals()*2
  } store_header;

  extern entitySet BcastEntitySet(entitySet set, int root, MPI_Comm comm);
  extern std::vector<int> simplePartitionVec(int mn, int mx, int p);
  extern std::vector<Loci::entitySet> simplePartition(int mn, int mx, MPI_Comm comm);
  extern  storeRepP Local2FileOrder(storeRepP sp, entitySet dom, int &offset,
                                    fact_db::distribute_infoP dist, MPI_Comm comm);
  

  extern void File2LocalOrder(storeRepP &result, entitySet resultSet,
                              storeRepP input, int offset,
                              fact_db::distribute_infoP dist,
                              MPI_Comm comm) ;

  extern int getMinFileNumberFromLocal(entitySet read_set,
                                       fact_db::distribute_infoP dist,
                                       MPI_Comm comm LOCI_DEFAULT_COMM) ;
  
  //only process 0 will call this function to write header, if ordered, also write domain eset
  void pmpi_WriteHeader(MPI_File fh, const entitySet &eset, const store_header& header );

  //only process 0 will call this function to read in header, if ordered, also read in domain eset
  void pmpi_ReadHeader( MPI_File fh,  store_header& header,  entitySet &eset);

  template< class T >
  inline void pmpi_writeStoreP(std::string& filename,                              
                               store<T> &var, const entitySet& write_set,
                               fact_db &facts, int iotype, bool ordered) {
    
    const_store<T> var_const ;
    var_const = var.Rep() ;
#ifndef MPI_STUBB
    pmpi_writeStoreP(filename, var_const, write_set, facts, iotype, ordered) ;
#else
    cerr << "MPI/IO is not supported with this compile of Loci!" << endl ;
    Loci::Abort() ;
#endif
  }
  
  template< class T >
  inline void pmpi_writeStoreVecP( std::string& filename,
                                   storeVec<T> &var,
                                   const entitySet& write_set, fact_db &facts, int xfer_type, bool ordered){
    const_storeVec<T> var_const;
    var_const = var.Rep() ;
    
#ifndef MPI_STUBB
    pmpi_writeStoreVecP(filename, var_const, write_set, facts, xfer_type, ordered) ;
#else
    cerr << "MPI/IO is not supported with this compile of Loci!" << endl ;
    Loci::Abort() ;
#endif
  }

  template< class T >
  inline void pmpi_write_ordered_store(std::string& filename, 
                                       const const_store<T> &var,int offset,
                                       fact_db &facts, int xfer_type) {
#ifndef MPI_STUBB
    /*
      write out the store that already in file numbering
    */
    MPI_Comm comm = facts.get_comm() ;
    int prank = 0 ;
    int np = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    MPI_Comm_size(comm,&np) ;

    
    //open for write only
    MPI_File fh = 0;
    MPI_File_open( comm, filename.c_str(),
                   MPI_MODE_WRONLY | MPI_MODE_CREATE, PHDF5_MPI_Info, &fh) ; 
   
 
    // Shift domain by offset
    entitySet dom  = var.domain();
    entitySet dom_file = dom >> offset ;

    // Compute overall domain across processors
    std::vector<entitySet> dom_vector = all_collect_vectors(dom_file,comm);
    entitySet q_dom;
    for(int i = 0; i < np; i++)
      q_dom += dom_vector[i];

    // If nothing to write, don't proceed
    if(q_dom == EMPTY){
      MPI_File_close(&fh);
      return ;
    }
    
   
    store_header header;
    header.ordered = 1;
    header.vec_size = 1;
    header.dom_size = 2*q_dom.num_intervals(); //how many int

    
    if(prank == 0)
      pmpi_WriteHeader(fh, q_dom, header);

    //need set moffset
    MPI_Offset moffset;
    moffset =(MPI_Offset) sizeof(header); //skip header
    moffset += (MPI_Offset) sizeof(int)*header.dom_size; //skip domain
   

    //pack data
    std::vector<T> data(dom.size()) ;
    int cnt = 0;
    FORALL(dom,ii) {
      data[cnt] = var[ii];
      cnt++ ;
    } ENDFORALL ;

    //write data
    pmpi_writeUnorderedVectorP(fh, moffset, data, comm, xfer_type) ;

    //close file
    MPI_File_close(&fh);
#else
    cerr << "MPI/IO is not supported with this compile of Loci!" << endl ;
    Loci::Abort() ;
#endif

  }

  template< class T >
  inline void pmpi_write_ordered_storeVec(std::string& filename, 
                                          const const_storeVec<T> &var,int offset,
                                          fact_db &facts, int xfer_type) {
#ifndef MPI_STUBB
    /*
      write out the storeVec that already in file numbering
    */
    MPI_Comm comm = facts.get_comm() ;
    int prank = 0 ;
    int np = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    MPI_Comm_size(comm,&np) ;

    
    //open for write only
    MPI_File fh = 0;
    MPI_File_open( comm, filename.c_str(),
                   MPI_MODE_WRONLY | MPI_MODE_CREATE, PHDF5_MPI_Info, &fh) ; 

   
 
    // Shift domain by offset
    entitySet dom  = var.domain();
    entitySet dom_file = dom >> offset ;

    // Compute overall domain across processors
    std::vector<entitySet> dom_vector = all_collect_vectors(dom_file,comm);
    entitySet q_dom;
    for(int i = 0; i < np; i++)
      q_dom += dom_vector[i];

    // If nothing to write, don't proceed
    if(q_dom == EMPTY){
      MPI_File_close(&fh);
      return ;
    }
    
    int vec_size = var.vecSize();
    
    store_header header;
    header.ordered = 1;
    header.vec_size = vec_size;
    header.dom_size = 2*q_dom.num_intervals(); //how many int

    
    if(prank == 0)
      pmpi_WriteHeader(fh, q_dom, header);

    //need set moffset
    MPI_Offset moffset;
    moffset =(MPI_Offset) sizeof(header); //skip header
    moffset += (MPI_Offset) sizeof(int)*header.dom_size; //skip domain
   

    //pack data
    std::vector<T> data(dom.size()*vec_size) ;
    int cnt = 0;
    FORALL(dom,ii) {
      for(int i = 0; i < vec_size; i++)
        data[cnt++] = var[ii][i];
    } ENDFORALL ;

    //write data
    pmpi_writeUnorderedVectorP(fh, moffset, data, comm, xfer_type) ;

    //close file
    MPI_File_close(&fh);
#else
    cerr << "MPI/IO is not supported with this compile of Loci!" << endl ;
    Loci::Abort() ;
#endif
  }

  
  template< class T >
  inline void pmpi_writeStoreP(std::string& filename,
                               const const_store<T> &var,
                               const entitySet& write_set, fact_db &facts, int xfer_type, bool ordered) {
#ifndef MPI_STUBB
    MPI_Comm comm = facts.get_comm() ;
    int prank = 0 ;
    int np = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    MPI_Comm_size(comm,&np) ;

      
    if(!ordered){
#ifdef io_performance
      MPI_Barrier(comm);
      stopWatch s;
      s.start();
#endif   

     
      
      //open file for write only
      MPI_File fh = 0;
      MPI_Offset moffset = 0;
      MPI_File_open( comm, filename.c_str(),
                     MPI_MODE_WRONLY | MPI_MODE_CREATE, PHDF5_MPI_Info, &fh) ; 
          
      entitySet dom = write_set ;
      fact_db::distribute_infoP df = facts.get_distribute_info() ;
      //just write out the entities belong to me
      if(df != 0)
        dom  = write_set & df->my_entities ;
   
      //pack fileids and data
     
      std::vector<int> fileids(dom.size()) ;
      std::vector<T> data(dom.size()) ;

      if(df == 0) {
        int cnt = 0 ;  
        FORALL(dom,ii) {
          int file_no = ii ;
          fileids[cnt] = file_no ;
          data[cnt] = var[ii];
          cnt++ ;
        } ENDFORALL ;

      } else {
        dMap l2f ;
        l2f = df->l2f.Rep() ;

        int cnt = 0 ;
        FORALL(dom,ii) {
          int file_no = l2f[ii] ;
          fileids[cnt] = file_no ;
          data[cnt] = var[ii];
          cnt++ ;
        } ENDFORALL ;
    
      }

   
      //process 0 write out header
      store_header header;
      header.ordered = 0;
      header.vec_size = 1;
      header.dom_size = 0; //how many int
      entitySet q_dom = EMPTY;
      if(prank == 0)
        pmpi_WriteHeader(fh, q_dom, header);

      //all processes update moffset 
      moffset =(MPI_Offset) sizeof(header); //skip header
      
      // write out two vectors fileids and data
      pmpi_writeUnorderedVectorP(fh, moffset, fileids, comm, xfer_type) ;
      pmpi_writeUnorderedVectorP(fh, moffset, data, comm, xfer_type) ;

      //close file
      MPI_File_close(&fh);
#ifdef io_performance  
      MPI_Barrier(comm);
      double wall_time = s.stop();
      if(prank == 0) std::cerr << "                                                    mpi parallel time to write unordered store : "  << wall_time << std::endl; 
#endif
    }else{//ordered writing
#ifdef io_performance
      MPI_Barrier(comm);
      stopWatch s;
      s.start();
#endif   
      int offset = 0 ;         
      Loci::fact_db::distribute_infoP dist = facts.get_distribute_info() ;
      if(dist == 0 || np == 1 ) { //no need to reorder, write directly
        pmpi_write_ordered_store(filename, var, offset, facts, xfer_type);
        return;
      }
     
      // parallel store write.. reorder to file numbering then write out

      // reorder from local to file numbering
      entitySet dom = var.domain() ;
      // Create container vardist that is ordered across processors in the
      // file numbering, the domain of this container shifted by offset
      // is the actual file numbering.
      // Loci::storeRepP vardist = Local2FileOrder(var.Rep(),dom,offset,dist,MPI_COMM_WORLD) ;
      const_store<T> out;
      out = Local2FileOrder(var.Rep(),dom,offset,dist,comm);
      // Write out container that has been distributed in the file numbering
      pmpi_write_ordered_store(filename, out, offset,facts, xfer_type) ;//group_id is replaced by fh
  
#ifdef io_performance  
      MPI_Barrier(comm);
      double wall_time = s.stop();
      if(prank == 0) std::cerr << "                                                    mpi parallel time to write   ordered store : "  << wall_time << std::endl; 
#endif
    }
#else
    cerr << "MPI/IO is not supported with this compile of Loci!" << endl ;
    Loci::Abort() ;
#endif

  }

  template< class T >
  inline void pmpi_writeStoreVecP( std::string& filename,
                                   const const_storeVec<T> &var,
                                   const entitySet& write_set, fact_db &facts, int xfer_type, bool ordered) {
#ifndef MPI_STUBB

    MPI_Comm comm = facts.get_comm() ;
    int prank = 0 ;
    int np = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    MPI_Comm_size(comm,&np) ;
    
    if(!ordered){
#ifdef io_performance
      MPI_Barrier(comm);
      Loci::stopWatch s ;
      s.start() ;
#endif

      //open file for write only
      MPI_File fh = 0;
      MPI_Offset moffset = 0;
      MPI_File_open( comm, filename.c_str(),
                     MPI_MODE_WRONLY | MPI_MODE_CREATE, PHDF5_MPI_Info, &fh) ; 
          
      entitySet dom = write_set ;
      fact_db::distribute_infoP df = facts.get_distribute_info() ;
      //just write out the entities belong to me
      if(df != 0)
        dom  = write_set & df->my_entities ;
    
      int vec_size = var.vecSize();

      // pack fileids and data 
      std::vector<int> fileids(dom.size()) ;
      std::vector<T> data(dom.size()*vec_size);
    
      if(df == 0) {
        int cnt = 0 ;
        FORALL(dom,ii) {
          int file_no = ii ;
          fileids[cnt] = file_no ;
          cnt++ ;
        } ENDFORALL ;

      } else {
        dMap l2f ;
        l2f = df->l2f.Rep() ;

        int cnt = 0 ;
        FORALL(dom,ii) {
          int file_no = l2f[ii] ;
          fileids[cnt] = file_no ;
          cnt++ ;
        } ENDFORALL ;
      
      }
    
      int cnt = 0;
      FORALL(dom,ii) {
        for(int i = 0; i < vec_size; i++)
          data[cnt++] = var[ii][i];  
      } ENDFORALL ;  
    
      //process 0 write out header
      store_header header;
      header.ordered = 0;
      header.vec_size = vec_size;
      header.dom_size = 0; //how many int
      entitySet q_dom = EMPTY;
      if(prank == 0)
        pmpi_WriteHeader(fh, q_dom, header);
      
      //all processes update moffset 
      moffset += (MPI_Offset)sizeof(header); //skip header

      //write out two vectors fileids and data
      pmpi_writeUnorderedVectorP(fh, moffset, fileids, comm, xfer_type) ;
      pmpi_writeUnorderedVectorP(fh, moffset, data, comm, xfer_type) ;

      //close file
      MPI_File_close(&fh);

#ifdef io_performance  
      MPI_Barrier(comm);
      double wall_time = s.stop();
      if(prank == 0) std::cerr << "                                                    mpi parallel time to write unordered storeVec : "  << wall_time << std::endl; 
#endif
      
    }else{//ordered writing
#ifdef io_performance
      MPI_Barrier(comm);
      stopWatch s;
      s.start();
#endif

      int offset = 0 ;         
      Loci::fact_db::distribute_infoP dist = facts.get_distribute_info() ;
      if(dist == 0 || np == 1 ) { //no need to reorder, write directly
        pmpi_write_ordered_storeVec(filename, var, offset, facts, xfer_type);
        return;
      }
      // parallel store write.. reorder to file numbering then write out

      // reorder from local to file numbering
      entitySet dom = var.domain() ;
      // Create container vardist that is ordered across processors in the
      // file numbering, the domain of this container shifted by offset
      // is the actual file numbering.
      
      const_storeVec<T> out;
      out = Local2FileOrder(var.Rep(),dom,offset,dist,comm);
      // Write out container that has been distributed in the file numbering
      pmpi_write_ordered_storeVec(filename, out, offset,facts, xfer_type) ;//group_id is replaced by fh
  
#ifdef io_performance  
      MPI_Barrier(comm);
      double wall_time = s.stop();
      if(prank == 0) std::cerr << "                                                    mpi parallel time to write   ordered storeVec : "  << wall_time << std::endl; 
#endif
    }
#else
    cerr << "MPI/IO is not supported with this compile of Loci!" << endl ;
    Loci::Abort() ;
#endif
  }


  template< class T >
  void pmpi_readStoreP(std::string& filename,
                       store<T> &var,
                       entitySet read_set,
                       fact_db &facts, int xfer_type) {
#ifndef MPI_STUBB
    MPI_Comm comm = facts.get_comm() ;
    int prank = 0 ;
    int np = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    MPI_Comm_size(comm,&np) ;
    
    //open file for read only
    MPI_File fh = 0;
    MPI_Offset moffset = 0;
    MPI_File_open( comm, filename.c_str(),
                   MPI_MODE_RDONLY, PHDF5_MPI_Info, &fh) ;
    
    //process 0 read in the header, if ordered, also read in domain
    entitySet q_dom;
    store_header header;
    if(prank == 0)pmpi_ReadHeader(fh, header, q_dom);
   
    //process 0 broadcast header
    MPI_Bcast(&header, sizeof(header), MPI_BYTE, 0, comm);

    //update moffset
    moffset =(MPI_Offset) sizeof(header); //skip header
    if(header.ordered) moffset += (MPI_Offset) sizeof(int)*header.dom_size; //skip domain

    if(!header.ordered){ //if store is unordered

      //read in fileID and data 
      std::vector<T> data ;
      std::vector<int> fileID ;
      pmpi_readUnorderedVectorP(fh, moffset, fileID, comm, xfer_type) ;
      pmpi_readUnorderedVectorP(fh, moffset, data, comm, xfer_type) ;

      //and close file
      MPI_File_close(&fh); 

      
      // Now we need to create a mapping from fileID to local and processor
      // number
      std::vector<int> local_num(fileID.size()) ;
      std::vector<int> procID(fileID.size()) ;
      findMapping(local_num,procID,read_set,fileID,facts) ;
     
      // const int p = MPI_processes ;
      std::vector<int> send_sz(np,0) ;
      std::vector<int> recv_sz(np,0) ;
      std::vector<int> recv_local_num ;
      distributeMapMultiStore(send_sz,recv_sz,recv_local_num,local_num,procID,comm) ;
     
      std::vector<T> recv_data ;
      sendData(recv_data,send_sz,recv_sz,recv_local_num,data,procID,comm) ;
  
     
      std::vector<int> alloc_set = recv_local_num ;
      std::sort(alloc_set.begin(),alloc_set.end()) ;
      entitySet alloc ;
      for(size_t i=0;i<alloc_set.size();++i)
        alloc += alloc_set[i] ;
     
    
      // if(alloc - read_set != EMPTY){
      //   cerr << " WARNING: var is not allocated on read_set" << endl;
      // };
     
      var.allocate(alloc) ;
      for(size_t i=0;i<recv_local_num.size();++i) {
        var[recv_local_num[i]] = recv_data[i] ;
      }
     
      return;
    }else{//if the store is ordered in filenumbering

      //read in data 
      std::vector<T> data ;
      pmpi_readUnorderedVectorP(fh, moffset, data, comm, xfer_type) ;
    
      //and close file
      MPI_File_close(&fh); 
      
      //process 0 broadcast domain
      q_dom = BcastEntitySet(q_dom,0,comm) ;
     
      //each process do partition and compute dom
      //std::vector<int> interval_sizes ;
      entitySet dom ;
      if(q_dom==EMPTY) {
        var.allocate(q_dom) ;
        return ;
      }
      
      {//ATTENTION: the data partition in pmpi_readUnorderedVectorP()
        //should be exactly the same as that in simplePartition
        std::vector<entitySet> ptn = simplePartition(q_dom.Min(),q_dom.Max(),comm) ;
        dom = ptn[prank] &q_dom ;
      }

     
      int offset = dom.Min() ;
      dom <<= offset ;
  
      //allocate qrep
      store<T> qvar ;
      qvar.allocate(dom) ;
      
      //unpack data, 
      int cnt = 0;
      FORALL(dom,ii) {
        qvar[ii] = data[cnt];
        cnt++ ;
      } ENDFORALL ;

      Loci::storeRepP new_store = qvar.Rep();
      // map from file number to local numbering
      fact_db::distribute_infoP dist = facts.get_distribute_info() ;
      if(dist != 0) {
        // Correct offset if file numbering changes.  Assume read_set is being
        // read in over the same set
        int minID = offset ;
        MPI_Bcast(&minID,1,MPI_INT,0,comm) ;
        const int minIDf = getMinFileNumberFromLocal(read_set,dist,comm) ;
        const int correct = minIDf - minID ;
        offset += correct  ;

        // Allocate space for reordered container
        Loci::storeRepP result = var.Rep()->new_store(read_set) ;
        Loci::File2LocalOrder(result,read_set,new_store,offset,dist,comm) ;
        // Copy results into container
        if(read_set == EMPTY) {
          read_set = result->domain() ;
          var.allocate(read_set) ;
        }
        var.Rep()->copy(result,read_set) ;
      } else {
        if(read_set == EMPTY) {
          read_set = new_store->domain() ;
          var.allocate(read_set) ;
        } else {
          offset = read_set.Min() - new_store->domain().Min() ;
          if(offset != 0) {
            // shift new store by offset to correct alignment
            new_store->shift(offset) ;
          }
        }
        var.Rep()->copy(new_store,read_set) ;
      }
    }
#else
    cerr << "MPI/IO is not supported with this compile of Loci!" << endl ;
    Loci::Abort() ;
#endif

  }


  template< class T >
  void pmpi_readStoreVecP(std::string& filename,
                          storeVec<T> &var,
                          entitySet read_set,
                          fact_db &facts, int xfer_type) {
#ifndef MPI_STUBB
    MPI_Comm comm = facts.get_comm() ;
    int prank = 0 ;
    int np = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    MPI_Comm_size(comm,&np) ;

    //open file for read only
    MPI_File fh = 0;
    MPI_Offset moffset = 0;
    MPI_File_open( comm, filename.c_str(),
                   MPI_MODE_RDONLY, PHDF5_MPI_Info, &fh) ; 
   


    //process 0 read in the header, if ordered, also read in domain
    
    entitySet q_dom;
    store_header header;
    if(prank == 0)pmpi_ReadHeader(fh, header, q_dom);
   
    //process 0 broadcast header
    MPI_Bcast(&header, sizeof(header), MPI_BYTE, 0, comm);

    int vec_size = header.vec_size;
    
    //update moffset
    moffset =(MPI_Offset) sizeof(header); //skip header
    if(header.ordered) moffset += (MPI_Offset) sizeof(int)*header.dom_size; //skip domain

    if(!header.ordered){ //if store is unordered
    
      //read in fileID and data
      std::vector<T> data ;
      std::vector<int> fileID ;
      pmpi_readUnorderedVectorP(fh, moffset, fileID, comm, xfer_type) ;
      data.resize(fileID.size()*vec_size);//must resize the vector to read in the data
      pmpi_readUnorderedDataVectorP(fh, moffset, data, comm, xfer_type) ;

      //and close file
      MPI_File_close(&fh); 
    
      // Now we need to create a mapping from fileID to local and processor
      // number
      std::vector<int> local_num(fileID.size()) ;
      std::vector<int> procID(fileID.size()) ;
      findMapping(local_num,procID,read_set,fileID,facts) ;

    
      std::vector<int> send_sz(np,0) ;
      std::vector<int> recv_sz(np,0) ;
      std::vector<int> recv_local_num ;
      distributeMapMultiStore(send_sz,recv_sz,recv_local_num,local_num,procID,comm) ;
    
      std::vector<T> recv_data ;
      sendMultiData(recv_data,send_sz,recv_sz,recv_local_num,data,procID, vec_size,comm) ;

   
      std::vector<int> alloc_set = recv_local_num ;
      std::sort(alloc_set.begin(),alloc_set.end()) ;
      entitySet alloc ;
      for(size_t i=0;i<alloc_set.size();++i)
        alloc += alloc_set[i] ;

  
      var.allocate(alloc) ;
      var.setVecSize(vec_size);
      int cnt = 0;
      for(size_t i=0;i<recv_local_num.size();++i) {
        for(int j = 0; j < vec_size; j++){
          var[recv_local_num[i]][j] = recv_data[cnt++] ;
        }
      }
      return;
    }else{//if the store is ordered in filenumbering

      //process 0 broadcast domain
      q_dom = BcastEntitySet(q_dom,0,comm) ;
     
      //each process do partition and compute dom
      //std::vector<int> interval_sizes ;
      entitySet dom ;
      if(q_dom==EMPTY) {
        var.allocate(q_dom) ;
        return ;
      }
      
      {
        std::vector<entitySet> ptn = simplePartition(q_dom.Min(),q_dom.Max(),comm) ;
        dom = ptn[prank] &q_dom ;
      }

     
      int offset = dom.Min() ;
      dom <<= offset ;

    
      //read in data, the size of vector data must be set first
      std::vector<T> data(dom.size()*vec_size) ;
      pmpi_readUnorderedDataVectorP(fh, moffset, data, comm, xfer_type) ;
    
      //and close file
      MPI_File_close(&fh); 
      
      //allocate qrep
      storeVec<T> qvar ;
      qvar.allocate(dom) ;
      qvar.setVecSize(vec_size);
      
      //unpack data, 
      int cnt = 0;
      FORALL(dom,ii) {
        for(int i = 0; i < vec_size; i++)
          qvar[ii][i] = data[cnt++];
      } ENDFORALL ;

      Loci::storeRepP new_store = qvar.Rep();
      // map from file number to local numbering
      fact_db::distribute_infoP dist = facts.get_distribute_info() ;
      if(dist != 0) {
        // Correct offset if file numbering changes.  Assume read_set is being
        // read in over the same set
        int minID = offset ;
        MPI_Bcast(&minID,1,MPI_INT,0,comm) ;
        const int minIDf = getMinFileNumberFromLocal(read_set,dist,comm) ;
        const int correct = minIDf - minID ;
        offset += correct  ;

        // Allocate space for reordered container
        Loci::storeRepP result = var.Rep()->new_store(read_set) ;
        Loci::File2LocalOrder(result,read_set,new_store,offset,dist,comm) ;
        // Copy results into container
        if(read_set == EMPTY) {
          read_set = result->domain() ;
          var.allocate(read_set) ;
          var.setVecSize(vec_size);
        }
        var.Rep()->copy(result,read_set) ;
      } else {
        if(read_set == EMPTY) {
          read_set = new_store->domain() ;
          var.allocate(read_set) ;
          var.setVecSize(vec_size);
        } else {
          offset = read_set.Min() - new_store->domain().Min() ;
          if(offset != 0) {
            // shift new store by offset to correct alignment
            new_store->shift(offset) ;
          }
        }
        var.Rep()->copy(new_store,read_set) ;
      }
    }
#else
    cerr << "MPI/IO is not supported with this compile of Loci!" << endl ;
    Loci::Abort() ;
#endif
  }
  
}
#endif
