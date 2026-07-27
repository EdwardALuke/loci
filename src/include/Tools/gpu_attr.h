#ifndef LOCI_GPU_ATTR_H
#define LOCI_GPU_ATTR_H

#ifdef __CUDACC__
#define GPU_DECL __host__ __device__
#else
#define GPU_DECL
#endif

#ifdef __CUDACC__
#define GPU_ONLY_DECL __device__
#else
#define GPU_ONLY_DECL
#endif

#endif
