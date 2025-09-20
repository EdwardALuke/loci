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
#ifndef MALLOC_ALLOC_H
#define MALLOC_ALLOC_H
namespace Loci {
  template <class T> class malloc_alloc ;

  template<> class malloc_alloc<void>
  {
  public:
    typedef void        value_type;
    typedef void*       pointer;
    typedef const void* const_pointer;
    
    template <class U> 
    struct rebind { typedef malloc_alloc<U> other; };
  };
  

  template <class T> class malloc_alloc
  {
  public:
    typedef T                 value_type;
    typedef value_type*       pointer;
    typedef const value_type* const_pointer;
    typedef value_type&       reference;
    typedef const value_type& const_reference;
#if defined(__sgi)
    typedef size_t       size_type;
    typedef ptrdiff_t    difference_type;
#else
    typedef std::size_t       size_type;
    typedef std::ptrdiff_t    difference_type;
#endif
    
    template <class U> 
    struct rebind { typedef malloc_alloc<U> other; };
    
    malloc_alloc() {}
    malloc_alloc(const malloc_alloc&) {}
    template <class U> 
    malloc_alloc(const malloc_alloc<U>&) {}
    ~malloc_alloc() {}
    
    pointer address(reference x) const { return &x; }
    const_pointer address(const_reference x) const { 
      return x;
    }

    pointer allocate(size_type n) 
    {
#if defined(__sgi)
      void* p = malloc(n * sizeof(T));
      if (!p)
        throw std::bad_alloc();
#else
      void* p = std::malloc(n * sizeof(T));
      if (!p)
        throw std::bad_alloc();
#endif
      return static_cast<pointer>(p);
    }

    pointer allocate(size_type n, typename Loci::malloc_alloc<void>::const_pointer  hint)  
    {
#if defined(__sgi)
      void* p = malloc(n * sizeof(T));
      if (!p)
        throw std::bad_alloc();
#else
      void* p = std::malloc(n * sizeof(T));
      if (!p)
        throw std::bad_alloc();
#endif
      return static_cast<pointer>(p);
    }
    
#if defined(__sgi)
    void deallocate(pointer p, size_type) { free(p); }
#else
    void deallocate(pointer p, size_type) { std::free(p); }
#endif
    
    size_type max_size() const { 
      return static_cast<size_type>(-1) / sizeof(T);
    }
    
    void construct(pointer p, const value_type& x) { 
      new(p) value_type(x); 
    }
    void destroy(pointer p) { p->~value_type(); }
    
  private:
    void operator=(const malloc_alloc&);
  };
  
  
  template <class T>
  inline bool operator==(const malloc_alloc<T>&, 
                         const malloc_alloc<T>&) {
    return true;
  }
  
  template <class T>
  inline bool operator!=(const malloc_alloc<T>&, 
                         const malloc_alloc<T>&) {
    return false;
  }

    
  
}
#endif
