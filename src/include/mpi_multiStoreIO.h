/*
  hdfio is replaced with mpiio

*/

#ifndef MPI_MULTISTOREIO_H
#define MPI_MULTISTOREIO_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif

//#ifdef H5_HAVE_PARALLEL
#include <dist_internal.h>
#include <distribute_container.h>
#include <distribute.h>
#include <distribute_io.h>
#include <multiStore.h>
#include "execute.h"
#include <multiStoreIO.h>
#include <mpi.h>

#ifndef MPI_STUBB
/*
  When writing in collective mode, the MPI library carries out a
  number of optimisations
  –  It uses fewer processes to actually do the writing
  •  Typically one per node
  –  It aggregates data in appropriate chunks before writing
*/


/*
  mpi version to read/write multiStores. parallel io only
  Notice: each variable is read in from /write into a separate file
  argument xfer_type is used to control collective or independent data transfer
  global variable PHDF5_MPI_Info is used.
*/
namespace Loci {   
  
  
  template<class T> void pmpi_writeUnorderedVectorP(MPI_File fh,
                                                    MPI_Offset &moffset, 
                                                    std::vector<T> &v,
                                                    MPI_Comm prime_comm,
                                                    int xfer_type) {
    /*this function write out the total size of v (one unsignd long) and the data in v
      in the order of process 0 , process 1, ...
      parallel io,
      moffset is used to specify the global location before writing ,
      and it is updated after writing
    */
    int mpi_size;
    int mpi_rank;
    MPI_Comm_rank(prime_comm, &mpi_rank);
    MPI_Comm_size(prime_comm, &mpi_size);

    //each process figure out the array_size_combined
    unsigned long local_size = v.size() ; 
    unsigned long rsize = local_size;
    std::vector<unsigned long> prime_count(mpi_size);  
    MPI_Allgather(&rsize,sizeof(unsigned long),MPI_BYTE,
                  &prime_count[0],sizeof(unsigned long),MPI_BYTE,prime_comm) ;

    unsigned long array_size_combined = 0;
    for(int i = 0; i < mpi_size; i++) array_size_combined += prime_count[i];
    
    if(array_size_combined == 0)
      return ;
  
  
    //process 0 write out array_size_combined
    if(mpi_rank==0){
      MPI_File_write_at(fh, moffset, &array_size_combined, sizeof(unsigned long), MPI_BYTE, MPI_STATUS_IGNORE);
    }
    //each process update moffset after write array_size_combined 
    moffset += (MPI_Offset)sizeof(unsigned long); 
    
   
    unsigned long start = 0;
    for(int i = 0; i < mpi_rank; i++)  start += prime_count[i];

    //to avoid communication between processes, when each process has a different offset, use  tmp_moffset
    MPI_Offset tmp_moffset = moffset;
    tmp_moffset += (MPI_Offset)start*sizeof(T);
    
    if(xfer_type == DXFER_COLLECTIVE_IO)
      MPI_File_write_at_all(fh, tmp_moffset, &v[0], rsize*sizeof(T), MPI_BYTE, MPI_STATUS_IGNORE);
    else
      MPI_File_write_at(fh, tmp_moffset, &v[0], rsize*sizeof(T), MPI_BYTE, MPI_STATUS_IGNORE);
  
    //update moffset to the end of whole dataset
    moffset += array_size_combined*sizeof(T);
  }

 
 
  template<class T> void pmpi_writeVectorSerialP(MPI_File fh,
                                                 MPI_Offset &moffset,
                                                 std::vector<T> &v ,
                                                 MPI_Comm prime_comm) {
    
    /*only process 0 write out data, independent only
      but v is available on all processes here
    */
    unsigned long array_size_combined = v.size() ; //if array_size_combined is not known to processes other than process 0, use MPI_Bcast instead of else statement 

    if(array_size_combined == 0)
      return ;

    if(MPI_COMM_NULL != prime_comm){
      int prime_rank = -1;
      MPI_Comm_rank(prime_comm, &prime_rank);
      if(prime_rank==0){
        MPI_File_write_at(fh,  moffset, &array_size_combined, sizeof(unsigned long), MPI_BYTE, MPI_STATUS_IGNORE);
        moffset += (MPI_Offset) sizeof(unsigned long);
        MPI_File_write_at(fh, moffset,  &v[0], array_size_combined*sizeof(T), MPI_BYTE, MPI_STATUS_IGNORE);
        moffset += (MPI_Offset) array_size_combined*sizeof(T); 
      }else{
        moffset +=  (MPI_Offset) (sizeof(unsigned long) + array_size_combined*sizeof(T)); 
      }
      
      // MPI_Bcast(&moffset, 1, MPI_OFFSET, 0, prime_comm); //process 0 is at the end of file now

    }
  }

  

  
  template<class T> void pmpi_readVectorSerialP(MPI_File fh,
                                                MPI_Offset& moffset,
                                                std::vector<T> &v,
                                                MPI_Comm prime_comm){
    /*
      prime_rank 0 read all the data (in this application, the block_schedule and block_sets,
      these are vector<int>, consider small data, thus process 0 read it in and broadcast it
      Notice: caller's responsibility to broadcast these data to other processes
    */
    
    if (MPI_COMM_NULL != prime_comm){
      int prime_rank = -1;
      MPI_Comm_rank(prime_comm, &prime_rank);
      if(prime_rank==0){
        unsigned long array_size_combined = 0;
        MPI_File_read_at(fh, moffset,  &array_size_combined, sizeof(unsigned long), MPI_BYTE, MPI_STATUS_IGNORE);
        moffset += (MPI_Offset)sizeof(unsigned long);
        if(array_size_combined > 0){
          v.resize(array_size_combined);
          MPI_File_read_at(fh, moffset, &v[0], array_size_combined*sizeof(T), MPI_BYTE, MPI_STATUS_IGNORE);
        }else {
          std::vector<T> tmp ;
          v.swap(tmp) ;
        }
        moffset += (MPI_Offset) array_size_combined*sizeof(T);
      }
      //need broadcast moffset to  all processes
      MPI_Bcast(&moffset, 1, MPI_OFFSET, 0, prime_comm); //process 0 is at the end of file now
    }
  }

  template<class T>
  void pmpi_readUnorderedVectorP(MPI_File fh, MPI_Offset& moffset, std::vector<T> &v, MPI_Comm prime_comm, int xfer_type) {
    /*
      each process read in different data in parallel
    */
    int mpi_rank = -1;
    int mpi_size = -1;
    
    if (MPI_COMM_NULL != prime_comm){
      MPI_Comm_rank(prime_comm, &mpi_rank);
      MPI_Comm_size(prime_comm, &mpi_size);

      //process 0 read in array_size_combined and broadcast it 
      unsigned long array_size_combined = 0;
      if(mpi_rank == 0)
        MPI_File_read_at(fh, moffset, &array_size_combined, sizeof(unsigned long), MPI_BYTE, MPI_STATUS_IGNORE);
      MPI_Bcast(&array_size_combined, sizeof(unsigned long), MPI_BYTE, 0, prime_comm); 
      

      
      moffset +=  (MPI_Offset)sizeof(unsigned long); //all the process update moffset
  
    
      //partition the data
      unsigned long delta = array_size_combined/mpi_size ;
      unsigned long rem = array_size_combined%mpi_size ;
      unsigned long loc_array_size = delta + ((mpi_rank<int(rem))?1:0);
  
  
      std::vector<T> tmp(loc_array_size) ;
      v.swap(tmp) ;
  
  
      //each process read in its own part in parallel here
      unsigned long start = (mpi_rank<int(rem))?mpi_rank*(delta+1):(delta+1)*rem+delta*(mpi_rank-rem) ;
      MPI_Offset tmp_moffset = moffset;
      tmp_moffset += (MPI_Offset) start*sizeof(T);
  
      unsigned long count = loc_array_size  ;
      
      if(xfer_type == DXFER_COLLECTIVE_IO)
        MPI_File_read_at_all(fh, tmp_moffset, &v[0], count*sizeof(T), MPI_BYTE, MPI_STATUS_IGNORE);
      else
        MPI_File_read_at(fh, tmp_moffset, &v[0], count*sizeof(T), MPI_BYTE, MPI_STATUS_IGNORE);

      //update moffset
      moffset += array_size_combined*sizeof(T);
    }
  }
  template<class T>
  void pmpi_readUnorderedDataVectorP(MPI_File fh, MPI_Offset& moffset, std::vector<T> &v, MPI_Comm prime_comm,  int xfer_type) {
    /*
      each process read in different data in parallel,
      this function is used when reading in storeVec data,
      v is already resized according to its domain,
      the data is partitioned according to the size of v
    */
    int mpi_rank = -1;
    int mpi_size = -1;
    
    if (MPI_COMM_NULL != prime_comm){
      MPI_Comm_rank(prime_comm, &mpi_rank);
      MPI_Comm_size(prime_comm, &mpi_size);

      //process 0 read in array_size_combined and broadcast it 
      unsigned long array_size_combined = 0;
      if(mpi_rank == 0)
        MPI_File_read_at(fh, moffset, &array_size_combined, sizeof(unsigned long), MPI_BYTE, MPI_STATUS_IGNORE);
     
      

      
      moffset +=  (MPI_Offset)sizeof(unsigned long); //all the process update moffset
  
    
     
  
      unsigned long loc_array_size = v.size(); 
      std::vector<unsigned long> all_local_sizes(mpi_size) ;
      MPI_Allgather(&loc_array_size,1,MPI_UNSIGNED_LONG,
                    &all_local_sizes[0],1,MPI_UNSIGNED_LONG,
                    prime_comm) ;

      //process 0 perform sanity check
      if(mpi_rank == 0) {
        unsigned long total = 0;
        for(int i = 0; i < mpi_size; i++) total += all_local_sizes[i];
        if(total != array_size_combined){
          cerr << "ERROR: total size of data vector doesn't match the size in file" << endl;
          Loci::debugout << "ERROR: total size of data vector doesn't match the size in file "<<endl;
          Loci::Abort();
        }
      }
  
      //each process read in its own part in parallel here
      unsigned long start = 0;
      for(int i = 0; i < mpi_rank; i++) start += all_local_sizes[i];
      MPI_Offset tmp_moffset = moffset;
      tmp_moffset += (MPI_Offset) start*sizeof(T);
  
      unsigned long count = loc_array_size  ;
      
      if(xfer_type == DXFER_COLLECTIVE_IO)
        MPI_File_read_at_all(fh, tmp_moffset, &v[0], count*sizeof(T), MPI_BYTE, MPI_STATUS_IGNORE);
      else
        MPI_File_read_at(fh, tmp_moffset, &v[0], count*sizeof(T), MPI_BYTE, MPI_STATUS_IGNORE);

      //update moffset
      moffset += array_size_combined*sizeof(T);
    }
  }

  
    
  template< class T >
  inline void pmpi_writeMultiStoreP(std::string filename,                              
                                    multiStore<T> &var, entitySet write_set,
                                    fact_db &facts, int iotype) {
    const_multiStore<T> var_const ;
    var_const = var.Rep() ;
    pmpi_writeMultiStoreP(filename, var_const,write_set, facts, iotype) ;
  }

  template< class T >
  void pmpi_writeMultiStoreP(std::string filename,
                             const const_multiStore<T> &var,
                             entitySet write_set, fact_db &facts, int xfer_type) {
    MPI_File fh = 0;
    MPI_Offset moffset = 0;
    MPI_File_open( MPI_COMM_WORLD, filename.c_str(),
                   MPI_MODE_WRONLY | MPI_MODE_CREATE, PHDF5_MPI_Info, &fh) ; 
      
    std::vector<int> sizes_local ;//unique sizes on this process, such as 0, 35, 105,...
      
    //    double time_write = 0 ;
    //    double pre_time = 0 ;
    
    entitySet dom = write_set ;
    fact_db::distribute_infoP df = facts.get_distribute_info() ;
    //just write out the entities belong to me???
    if(df != 0)
      dom  = write_set & df->my_entities ;

    Loci::stopWatch sp ;

    sp.start() ;
    { // Write out book keeping stuff.
      std::vector<int> counts(dom.size()) ;
      int cnt = 0 ;
      std::vector<int> fileids(dom.size()) ;
      if(df == 0) {
        FORALL(dom,ii) {

          counts[cnt] = var.vec_size(ii) ;
          int file_no = ii ;
          fileids[cnt] = file_no ;
          cnt++ ;
        } ENDFORALL ;

      } else {
        dMap l2f ;
        l2f = df->l2f.Rep() ;

        FORALL(dom,ii) {
          counts[cnt] = var.vec_size(ii) ;

          int file_no = l2f[ii] ;
          fileids[cnt] = file_no ;
          cnt++ ;
        } ENDFORALL ;

      }

      
      //no groups, write out two vectors counts and fileids
      pmpi_writeUnorderedVectorP(fh, moffset, counts, MPI_COMM_WORLD, xfer_type) ;
      pmpi_writeUnorderedVectorP(fh, moffset, fileids, MPI_COMM_WORLD, xfer_type) ;
            
      std::sort(counts.begin(),counts.end()) ;
      std::vector<int>::const_iterator lu = std::unique(counts.begin(),counts.end()) ;
      int nunique = lu-counts.begin() ;
      std::vector<int> tmp(nunique) ;
      sizes_local.swap(tmp) ;
      for(int i=0;i<nunique;++i)
        sizes_local[i] = counts[i] ;
    }

    // schedule blocking for writing

    // Find global sizes array
    sizes_local.push_back(std::numeric_limits<int>::max()) ;
    int iloc = 0 ;
    std::vector<int> sizes ; //unique sizes for all processors
    int cmin = -1 ;
    do {
      MPI_Allreduce(&sizes_local[iloc],&cmin,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD) ;
      if(cmin != std::numeric_limits<int>::max())
        sizes.push_back(cmin) ;
      if(sizes_local[iloc] == cmin)
        iloc++ ;
    } while (cmin != std::numeric_limits<int>::max()) ;

    std::vector<int> block_sizes ;//for all processes
    size_t cnt = 0 ;
    if(sizes[0] == 0)
      cnt++ ;
    const int blk_factor = 256 ; //actual block size is Min(internal block size, 256)
    int cloc = 0 ;
    int cloc_prev = 0 ;
    while(cnt < sizes.size()) {
      cloc += blk_factor ;
      cloc = min(cloc,sizes[cnt]) ;
      block_sizes.push_back(cloc-cloc_prev) ;
      cloc_prev = cloc ;
      if(cloc == sizes[cnt])
        cnt++ ;
    }
    if(block_sizes.size() == 0) {
      block_sizes.push_back(0) ;
    }
    std::vector<int> num_local_xfers(block_sizes.size(),0) ; //local::for each size, how many elements need xfer

    FORALL(dom,ii) {
      int lsz = var.vec_size(ii) ;
      int tot = 0 ;
      for(size_t i=0;i<block_sizes.size();++i) {
        tot += block_sizes[i] ;
        if(lsz >= tot && block_sizes[i] != 0)
          num_local_xfers[i]++ ;
      }
    } ENDFORALL ;

    std::vector<int> block_data_elems(block_sizes.size(),0) ;//global::for each block, how many elements need xfer
    MPI_Allreduce(&num_local_xfers[0],&block_data_elems[0],block_sizes.size(),
                  MPI_INT,MPI_SUM,MPI_COMM_WORLD) ;

    // write out write schedule information
   
    pmpi_writeVectorSerialP(fh, moffset, block_sizes, MPI_COMM_WORLD) ;
    pmpi_writeVectorSerialP(fh, moffset, block_data_elems, MPI_COMM_WORLD) ;
   
    //    pre_time += sp.stop() ;
    // Now write out the main data block


    unsigned long total_size = 0 ;//for all processes, all blocks
    for(size_t i=0;i<block_sizes.size();++i) {
      total_size += block_sizes[i]*block_data_elems[i] ;
    }

    int p = 0;
    int r = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &r);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    if(r == 0)MPI_File_write_at(fh, moffset, &total_size, sizeof(unsigned long), MPI_BYTE, MPI_STATUS_IGNORE);
    moffset += (MPI_Offset)sizeof(unsigned long);
    


    int max_local_size = 0 ;
    for(size_t i=0;i<block_sizes.size();++i) {
      int local_size = block_sizes[i]*num_local_xfers[i] ;
      max_local_size = max(max_local_size,local_size) ;
    }
   
    //allgather local sizes over
    std::vector<int> all_local_sizes(block_sizes.size()*p) ;
    MPI_Allgather(&num_local_xfers[0],block_sizes.size(),MPI_INT,
                  &all_local_sizes[0],block_sizes.size(),MPI_INT,
                  MPI_COMM_WORLD) ;

 
    /* Create a large dataset for all processes  */
  
    //send_buffer  over all blocks
    std::vector<T> send_buffer(max_local_size) ;

    int loc = 0 ;
    unsigned long start = 0 ;
 
    for(size_t i=0;i<block_sizes.size();++i) { //write out data block by block
      int bs = block_sizes[i] ;

      //prepare send_buffer
      if(num_local_xfers[i] != 0) {
        int bcnt = 0 ;
        FORALL(dom,ii) {
          int lsz = var.vec_size(ii) ;
          if(lsz >= loc+bs) {
            const T *base = var.begin(ii) ;
            for(int j=0;j<bs;++j) {
              send_buffer[bcnt++] = base[loc+j] ;
            }
          }
        } ENDFORALL ;
      }
      loc += bs ;
      {
        int local_size =  all_local_sizes[r*block_sizes.size()+i]*bs;
        //start depend on prime_rank
        int offset = 0;
        for(int pr = 0; pr < r ; pr++){
          offset +=  all_local_sizes[pr*block_sizes.size()+i]*bs;
        }
        unsigned long total_block_size = 0; //to update moffset
        for(int pr = 0; pr < p ; pr++)total_block_size +=  all_local_sizes[pr*block_sizes.size()+i]*bs;
        
        Loci::stopWatch s ;
        s.start() ;
        unsigned long prime_start = start + offset;
      
        {
          Loci::stopWatch s ;
          s.start() ;
 
          unsigned long count = local_size *sizeof(T);
          // moffset += (MPI_Offset) prime_start * sizeof(T);
          MPI_Offset tmp_moffset = moffset + (MPI_Offset) (prime_start * sizeof(T));
  
          int lsz = (count==0)?0:1 ;
          int gsz  = lsz ;

          MPI_Allreduce(&lsz,&gsz,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
          if(gsz != 0) {

            if(xfer_type == DXFER_COLLECTIVE_IO)
              MPI_File_write_at_all(fh, tmp_moffset, &send_buffer[0], count, MPI_BYTE, MPI_STATUS_IGNORE);
            else
              MPI_File_write_at(fh, tmp_moffset, &send_buffer[0], count, MPI_BYTE, MPI_STATUS_IGNORE);

            //moffset += (MPI_Offset)count;
            //MPI_Bcast(&moffset, 1, MPI_OFFSET, p-1, MPI_COMM_WORLD); //process mpi_size-1 is at the end of file now
            moffset += (MPI_Offset) total_block_size*sizeof(T);
          }
        }
	//        time_write += s.stop() ;
        start += block_sizes[i]*block_data_elems[i]; //how many total elements written for this block
      }
    }
    MPI_File_close(&fh);
  }
  
 
  template< class T >
  void pmpi_readMultiStoreP(std::string& filename,
                            multiStore<T> &var,
                            entitySet read_set,
                            fact_db &facts, int xfer_type) {

    //group_id is known to the world
    MPI_File fh = 0;
    MPI_Offset moffset = 0;
    MPI_File_open( MPI_COMM_WORLD, filename.c_str(),
                   MPI_MODE_RDONLY, PHDF5_MPI_Info, &fh) ; 
    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    int r = MPI_rank ;
    
   
    std::vector<int> counts ;
    std::vector<int> fileID ;
    pmpi_readUnorderedVectorP(fh, moffset, counts, MPI_COMM_WORLD, xfer_type) ;
    pmpi_readUnorderedVectorP(fh, moffset, fileID, MPI_COMM_WORLD, xfer_type) ;
   
    // Now we need to create a mapping from fileID to local and processor
    // number
    std::vector<int> local_num(fileID.size()) ;
    std::vector<int> procID(fileID.size()) ;
    findMapping(local_num,procID,read_set,fileID,facts) ;

    const int p = MPI_processes ;
    std::vector<int> send_sz(p,0) ;
    std::vector<int> recv_sz(p,0) ;
    std::vector<int> recv_local_num ;
    distributeMapMultiStore(send_sz,recv_sz,recv_local_num,local_num,procID) ;
    std::vector<int> recv_count ;
    sendData(recv_count,send_sz,recv_sz,recv_local_num,counts,procID) ;
  

    std::vector<int> alloc_set = recv_local_num ;
    std::sort(alloc_set.begin(),alloc_set.end()) ;
    entitySet alloc ;
    for(size_t i=0;i<alloc_set.size();++i)
      alloc += alloc_set[i] ;

    store<int> count_store ;
    count_store.allocate(alloc) ;
    for(size_t i=0;i<recv_local_num.size();++i) {
      count_store[recv_local_num[i]] = recv_count[i] ;
    }
    var.allocate(count_store) ;

    // Now with the background work done we can begin reading in the
    // multiStore in segments, but we need to get the schedule ready
    // schedule of the blocking sizes
    std::vector<int> block_schedule ;
    // sizes of the subarrays for each block
    std::vector<int> block_sets ;


    pmpi_readVectorSerialP(fh, moffset, block_schedule, MPI_COMM_WORLD) ;
    pmpi_readVectorSerialP(fh, moffset, block_sets, MPI_COMM_WORLD) ;
    int bsize = block_schedule.size() ;
    MPI_Bcast(&bsize,1,MPI_INT,0,MPI_COMM_WORLD) ;
    if(0 != r) {
      std::vector<int> bsched(bsize) ;
      block_schedule.swap(bsched) ;
      std::vector<int> bsets(bsize) ;
      block_sets.swap(bsets) ;
    }
    MPI_Bcast(&block_schedule[0],bsize,MPI_INT,0,MPI_COMM_WORLD) ;
    MPI_Bcast(&block_sets[0],bsize,MPI_INT,0,MPI_COMM_WORLD) ;

    //prime processes Open data array
  
    unsigned long dimension = 0 ;
    //read in dimension here
    if(r == 0) MPI_File_read_at(fh, moffset, &dimension, sizeof(unsigned long), MPI_BYTE, MPI_STATUS_IGNORE);
    MPI_Bcast(&dimension,sizeof(unsigned long), MPI_BYTE,0,MPI_COMM_WORLD) ;
    moffset += (MPI_Offset)sizeof(unsigned long);

    

    // Loop over block schedule, communicate each block
    
    unsigned long start = 0 ;
 
 
    int blksz = 0 ;
    std::vector<int> file_read_sizes(p,0) ;
    for(size_t i=0;i<block_sets.size();++i) {
      //for each block
      blksz += block_schedule[i] ;
      int lbsz = 0 ;
      std::vector<int> send_sz_blk(p,0) ;
      for(size_t j=0;j<counts.size();++j)
        if(counts[j] >=blksz && blksz > 0) {
          lbsz++ ;
          send_sz_blk[procID[j]]++ ;
        }

      MPI_Allgather(&lbsz,1,MPI_INT,&file_read_sizes[0],1,MPI_INT,MPI_COMM_WORLD) ;

      int tot = 0 ;
      for(int j=0;j<p;++j)
        tot+= file_read_sizes[j] ;
      if(r == 0 && tot != block_sets[i]) {
        std::cerr << "inconsistent, tot = "
                  << tot << " block_sets[" << i << "]=" << block_sets[i]
                  << endl ;
      }

      int msgsz = lbsz*block_schedule[i] ;

      std::vector<T> read_buffer ;
      if(msgsz > 0) read_buffer.resize(msgsz);

      {
        unsigned long offset = 0;
        for(int pr = 0; pr < r ; pr++){
          offset += file_read_sizes[pr]*block_schedule[i] ;
        }
        unsigned long total_block_size = 0;//for update moffset
        for(int pr = 0; pr < p ; pr++)total_block_size += file_read_sizes[pr]*block_schedule[i];

        unsigned long prime_start = start + offset;

        unsigned long count = msgsz;
     

        /* set up the collective transfer properties list */

     

        MPI_Offset tmp_moffset = moffset + (MPI_Offset) prime_start*sizeof(T);
        if(xfer_type == DXFER_COLLECTIVE_IO)
          MPI_File_read_at_all(fh, tmp_moffset, &read_buffer[0], count*sizeof(T), MPI_BYTE, MPI_STATUS_IGNORE);
        else
          MPI_File_read_at(fh, tmp_moffset, &read_buffer[0], count*sizeof(T), MPI_BYTE, MPI_STATUS_IGNORE);
        moffset += (MPI_Offset)total_block_size *sizeof(T);
        start += block_sets[i]*block_schedule[i];////how many total elements read for this block
      }

      // After reading in the data, send it to the appropriate processor
      //First figure out what needs to be sent

      // copy data to send ordered on processor
      std::vector<T> xmit_buffer ; //give it a default capacity
      if(msgsz > 0)xmit_buffer.resize(msgsz);

      std::vector<int> cnt(p,0) ;
      std::vector<int> xoffset(p+1,0) ;
      for(int j=0;j<p;++j)
        xoffset[j+1] = xoffset[j]+send_sz_blk[j] ;

      int rcnt = 0 ;
      for(size_t j=0;j<counts.size();++j)
        if(counts[j] >=blksz) {
          int ps = procID[j] ;
          int loc = xoffset[ps]+cnt[ps] ;
          int bs = block_schedule[i] ;
          for(int k=0;k<bs;++k)
            xmit_buffer[loc*bs+k] = read_buffer[rcnt++] ;
          cnt[ps]++ ;
        }

      std::vector<int> recv_sz_blk(p,0) ;
      int indx = 0 ;
      for(int j=0;j<p;++j)
        for(int k=0;k<recv_sz[j];++k,indx++) {
          int c = count_store[recv_local_num[indx]] ;
          if(c >= blksz && c != 0) {
            recv_sz_blk[j]++ ;
          }
        }
      std::vector<int> roffset(p+1,0) ;
      for(int j=0;j<p;++j)
        roffset[j+1] = roffset[j]+recv_sz_blk[j] ;
      std::vector<T> rdata(roffset[p]*block_schedule[i]) ;
      int nreq = 0 ;
      for(int j=0;j<p;++j) {
        if(send_sz_blk[j] > 0)
          nreq++ ;
        if(recv_sz_blk[j] > 0)
          nreq++ ;
      }

      std::vector<MPI_Request> recv_Requests(nreq) ;
      int req = 0 ;
      for(int j=0;j<p;++j)
        if(recv_sz_blk[j] > 0) {
          MPI_Irecv(&rdata[roffset[j]*block_schedule[i]],
                    recv_sz_blk[j]*block_schedule[i]*sizeof(T),
                    MPI_BYTE,j,3,
                    MPI_COMM_WORLD,&recv_Requests[req]) ;
          req++ ;
        }
      for(int j=0;j<p;++j)
        if(send_sz_blk[j] > 0) {
          MPI_Isend(&xmit_buffer[xoffset[j]*block_schedule[i]],
                    send_sz_blk[j]*block_schedule[i]*sizeof(T),
                    MPI_BYTE,j,3,
                    MPI_COMM_WORLD,&recv_Requests[req]) ;
          req++ ;
        }
      std::vector<MPI_Status> statuslist(nreq) ;
      MPI_Waitall(nreq,&recv_Requests[0],&statuslist[0]) ;

      // Now copy data into multistore
      indx = 0 ;
      int l = 0 ;
      int bs = block_schedule[i] ;
      for(int j=0;j<p;++j)
        for(int k=0;k<recv_sz[j];++k,indx++) {
          int le = recv_local_num[indx] ;
          if(count_store[le] >= blksz) {
            for(int m=0;m<bs;++m) {
              var[le][blksz+m-bs] = rdata[l++] ;
            }
          }
        }
    }
    MPI_File_close(&fh);
  }
}
#endif

#endif
