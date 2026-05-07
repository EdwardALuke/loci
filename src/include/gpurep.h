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
#ifndef GPUREP_H
#define GPUREP_H

// This header file contains the class definition of
// gpustoreRepI, store, and const_store.
// Their corresponding template implementation is in
// store_impl.h
// This separation is necessary to resolve some dependency
// problems in the class hierarchy.
// The same design applies to all other container classes.
#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <DMap.h>
#include <store_rep.h>
#include <data_traits.h>
#include <sstream>
#include <hdf5_readwrite.h>
#include <mpi.h>
#include <string.h>
#include <dist_internal.h>

#ifdef USE_CUDA_RT
#include <cuda_runtime_api.h>
#endif

namespace Loci {

#ifdef USE_CUDA_RT

  extern int MAXGPUStreamAlloc ;
  extern int GPUStreamAlloc ;

  inline int getGPUStreamID() {
    int retval = GPUStreamAlloc ;
    GPUStreamAlloc = (retval + 1) & (MAXGPUStreamAlloc-1) ;
    return retval ;
  }
  
  extern cudaStream_t streamSet[] ;

  inline cudaStream_t &getGPUStream() {
    return streamSet[getGPUStreamID()] ;
  }
  
#endif
  
  class gpuRep : public storeRep {
  public:
    gpuRep() { } ;
    ~gpuRep() {} ;
    // Code to copy from cpu container to gpu container
    virtual void copyFrom(const storeRepP &p, entitySet set) = 0 ;
    // code to copy from gpu container to cpu container
    virtual void copyTo(storeRepP &p, entitySet set) const = 0 ;
  } ;

  typedef NPTR<gpuRep> gpuRepP ;


  // GPU store allocation object
  struct GPUstoreAllocateInfo {
    void *alloc_ptr1 ;
    void *alloc_ptr2 ;
    void *base_ptr ;
    int base_offset ;
    int size ;
    size_t allocated_size ;
    bool allocated ;
    entitySet allocset ;

    GPUstoreAllocateInfo()
      : alloc_ptr1(0), alloc_ptr2(0), base_ptr(0),
        base_offset(0), size(0), allocated_size(0),
        allocated(false) {
      allocset = EMPTY ;
    }

    template<class T> void release() {
      if(alloc_ptr1 != 0) {
        // Call placement delete
        if(!std::is_trivially_default_constructible<T>::value) {
          cerr << "gpu container types should be trivially constructable"
               << endl ;
          Loci::Abort() ;
        }

#ifdef USE_CUDA_RT
        cudaError_t err = cudaFree(alloc_ptr1) ;
        if(err!= cudaSuccess) {
          cerr << "cudaFree failed" << endl ;
          cerr << "error = " << cudaGetErrorString(err) << endl;
          Loci::Abort() ;
        }
#else
        free(alloc_ptr1) ;
#endif

        if(alloc_ptr2 != 0) {
#ifdef USE_CUDA_RT
          err = cudaFree(alloc_ptr2) ;
          if(err!= cudaSuccess) {
            cerr << "cudaFree failed" << endl ;
            cerr << "error = " << cudaGetErrorString(err) << endl;
            Loci::Abort() ;
          }
#else
          free(alloc_ptr2) ;
#endif
        }

        alloc_ptr1 = 0 ;
        alloc_ptr2 = 0 ;
        base_ptr = 0 ;
        base_offset = 0 ;
        size = 0 ;
        allocated_size = 0 ;
        allocated = false ;
        allocset = EMPTY ;
      }
    }

    template<class T>
    void allocBasic(const entitySet &eset, int sz) {
      // if the pass in is EMPTY, we delete the previous allocated memory
      // this equals to free the memory
      if( eset == EMPTY ) {
        // just return instead of releasing memory
        return ;

        if(alloc_ptr1 != 0) {
          // Call placement delete
          if(!std::is_trivially_default_constructible<T>::value) {
            cerr << "gpu container types should be trivially constructible"
                 << endl ;
            Loci::Abort() ;
          }

#ifdef USE_CUDA_RT
          cudaError_t err = cudaFree(alloc_ptr1) ;
          if(err!= cudaSuccess) {
            cerr << "cudaFree failed" << endl ;
            cerr << "error = " << cudaGetErrorString(err) << endl;
            Loci::Abort() ;
          }
#else
          free(alloc_ptr1) ;
#endif
        }

        if(alloc_ptr2 != 0) {
#ifdef USE_CUDA_RT
          cudaError_t err = cudaFree(alloc_ptr2) ;
          if(err != cudaSuccess) {
            cerr << "cudaFree failed" << endl ;
            cerr << "error = " << cudaGetErrorString(err) << endl ;
            Loci::Abort() ;
          }
#else
          free(alloc_ptr2) ;
#endif
        }

        alloc_ptr1 = 0 ;
        alloc_ptr2 = 0 ;
        base_ptr = 0 ;
        size = sz ;
        base_offset = 0 ;
        allocated = false ;
        allocated_size = 0 ;
        allocset = eset ;
      } else {
        int_type old_range_min = allocset.Min() ;
        int_type old_range_max = allocset.Max() ;
        int_type new_range_min = eset.Min() ;
        int_type new_range_max = eset.Max() ;

        // if the old range and the new range are equal, nothing
        // needs to be done, just return
        if( (old_range_min == new_range_min) &&
            (old_range_max == new_range_max) &&
            size == sz) {
          allocset = eset ;
          //	cerr << "doing nothing... allocset=" << allocset << ", eset=" << eset  << endl ;
          return ;
        }

        size_t alloc_size = (new_range_max-new_range_min+1)*sz ;

        // TODO: copying overlapped entity data is currently disabled for GPU
        // containers.

        // is there any overlap between the old and the new domain?
        // we copy the contents in the overlap region to the new
        // allocated storage
        //entitySet ecommon = allocset & eset ;

#ifdef USE_CUDA_RT
        T * tmp_alloc_pointer ;
        cudaError_t err = cudaMalloc((void **) &tmp_alloc_pointer,sizeof(T)*alloc_size) ;
        if(err != cudaSuccess) {
          cerr << "cudaMalloc  failed" << endl ;
          cerr << "error = " << cudaGetErrorString(err) << endl;
          Loci::Abort() ;
        }
        cudaDeviceSynchronize() ;
#else
        T * tmp_alloc_pointer = (T *) malloc(sizeof(T)*(alloc_size)+(STORE_ALIGN_SIZE)) ;
#endif

        T * tmp_base_ptr = tmp_alloc_pointer ;
        T * tmp_base_algn = (T *) ((uintptr_t) tmp_base_ptr ) ;
        if(tmp_base_ptr != tmp_base_algn) {
          tmp_base_ptr = (T *) ( tmp_base_algn) ;
        }

        // Call placement new
        if(!std::is_trivially_default_constructible<T>::value) {
          cerr << "gpu container types should be trivially constructable" << endl ;
          Loci::Abort() ;
        }

        if(sz == size) {
          //Not copying on realloc
          //	T *p = ((T *) base_ptr)  ;
          // if ecommon == EMPTY, then nothing is done in the loop
          //	FORALL(ecommon,ii) {
          //	  for(int i=0;i<sz;++i)
          //	    tmp_base_ptr[ii-new_range_min*sz+i] = p[ii-base_offset*sz+i] ;
          //	} ENDFORALL ;
        }

        // Not calling placement new/delete, or copying on realloc
        if(!std::is_trivially_default_constructible<T>::value) {
          cerr << "gpu container types should be trivially constructable" << endl ;
          Loci::Abort() ;
        }

        if(alloc_ptr1 != 0) {
#ifdef USE_CUDA_RT
          cudaError_t err = cudaFree(alloc_ptr1) ;
          if(err != cudaSuccess) {
            cerr << "cudaFree  failed" << endl ;
            cerr << "error = " << cudaGetErrorString(err) << endl;
            Loci::Abort() ;
          }
#else
          free(alloc_ptr1) ;
#endif
        }

        if(alloc_ptr2 != 0) {
#ifdef USE_CUDA_RT
          cudaError_t err = cudaFree(alloc_ptr2) ;
          if(err != cudaSuccess) {
            cerr << "cudaFree failed" << endl ;
            cerr << "error = " << cudaGetErrorString(err) << endl ;
            Loci::Abort() ;
          }
#else
          free(alloc_ptr2) ;
#endif
        }

        alloc_ptr1 = tmp_alloc_pointer ;
        alloc_ptr2 = 0 ;
        base_ptr = tmp_base_ptr ;
        base_offset = new_range_min ;
        size = sz;
        allocated_size = alloc_size ;
        allocated = true ;
        allocset = eset ;
      }
    }

    // elist is the list of entities (in increasing order)
    // clist is the list of size for each of these entities
    // ptn is the final allocated size
    template<class T> void allocMulti(const storeAllocateInfo &count,
				      entitySet ptn) {
      if(count.base_ptr ==0)
	return ; // count must be valid
      size_t sum = 0 ;
      entitySet context = ptn & count.allocset ;
      FORALL(context,ii) {
	sum += ((int *)count.base_ptr)[ii-count.base_offset] ;
      } ENDFORALL ;

      size_t npntrs = 0 ;
      if(ptn != EMPTY)
	npntrs = (ptn.Max()-ptn.Min()+2) ; // Always allocate one more so
      // we can get the size of the last element.
      // Ok, alloc_ptr1 is the pointer allocated for the data
      // alloc_ptr2 is the pointer allocated for the pointer array
      // this will be used for accessing data in the multiStore
      // The base_ptr will be used to call the in place constructor/destructor
      // for the type allocation
      T * tmp_alloc_ptr1 = 0 ;
      T * tmp_base_ptr = 0 ;
      T ** tmp_alloc_ptr2 = 0 ;
      size_t tmp_allocated_sz = sum + 1 ; 
      
      tmp_alloc_ptr1 = (T *) malloc(sizeof(T)*tmp_allocated_sz +(STORE_ALIGN_SIZE)) ;
      tmp_base_ptr = tmp_alloc_ptr1 ;
      T * tmp_base_algn = (T *) ((uintptr_t) tmp_base_ptr & ~(uintptr_t)(STORE_ALIGN_SIZE-1)) ;
      if(tmp_base_ptr !=tmp_base_algn) 
	tmp_base_ptr = (T *) ((uintptr_t) tmp_base_algn+(uintptr_t)STORE_ALIGN_SIZE) ;
      // Call placement new
      if(!std::is_trivially_default_constructible<T>::value) {
	cerr << "gpu container types should be trivially constructable"
	     << endl ;
	Loci::Abort() ;
      }
      tmp_alloc_ptr2 = (T **) malloc(sizeof(T*)*npntrs) ;
      int tmp_base_offset = ptn.Min() ;
      if(ptn==EMPTY)
	tmp_base_offset = 0 ;
      size_t loc = 0 ;
      FORALL(ptn,ii) {
	tmp_alloc_ptr2[ii-tmp_base_offset] = tmp_base_ptr+loc ;
	if(count.allocset.inSet(ii))
	  loc += ((int *)count.base_ptr)[ii-count.base_offset] ;
	tmp_alloc_ptr2[ii-tmp_base_offset+1] = tmp_base_ptr+loc ;
      } ENDFORALL ;

      // release existing memory
      release<T>() ;
      
      alloc_ptr1 = tmp_alloc_ptr1 ;
      alloc_ptr2 = tmp_alloc_ptr2 ;
      base_ptr = tmp_base_ptr ;
      base_offset = tmp_base_offset ;
      allocated_size = tmp_allocated_sz ;
      size = -1 ;
      allocset = ptn ;
      return ;
    }
    
  } ;


  extern std::vector<GPUstoreAllocateInfo> GPUstoreAllocateData ;

  extern int getGPUStoreAllocateID() ;
  extern void releaseGPUStoreAllocateID(int id) ;
  

}
#endif
