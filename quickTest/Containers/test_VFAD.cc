#include <Loci.h>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>

#include <iostream>
#include <random>
#include <vector>
#include <algorithm>

using std::max ;
using std::min ;
using std::fabs ;
using std::abs ;
using namespace Loci ;

#define TOLCHECK(a,b,tol) (a == doctest::Approx(b).epsilon(tol))
const double tol = 1e-7 ;

TEST_CASE("VFAD arithmetic operations") {
  std::random_device rd;
  std::default_random_engine re(rd());
  std::uniform_real_distribution<double> rdist(-10, 10);
  
  typedef VFAD type;
  
  unsigned int const N = 100;
  
  std::vector<type> a(N), b(N), c(N);
  std::vector<double> ad(N), bd(N) ;
  for(unsigned int e = 0; e < N; ++e) {
    a[e] = rdist(re);
    b[e] = rdist(re);
    if(b[e].data.value == 0) // Don't let b be zero so it can test a/b
      b[e].data.value = 1.0 ;
    ad[e] = realToDouble(a[e]) ;
    bd[e] = realToDouble(b[e]) ;
    for(unsigned int i=0;i<type::maxN;++i) {
      a[e].data.grad[i] = rdist(re) ;
      b[e].data.grad[i] = rdist(re) ;
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = a[e]+b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == realToDouble(a[e])+realToDouble(b[e])) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      FADd B(b[e].data.value,b[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(A + B).grad,tol)) ;
    }
  }

  for(unsigned int e = 0; e < N; ++e) {
    c[e] = a[e]-b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == realToDouble(a[e])-realToDouble(b[e])) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      FADd B(b[e].data.value,b[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(A - B).grad,tol)) ;
    }
  }

  for(unsigned int e = 0; e < N; ++e) {
    c[e] = a[e]*b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == realToDouble(a[e])*realToDouble(b[e])) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      FADd B(b[e].data.value,b[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(A * B).grad,tol)) ;
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = a[e]/b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == realToDouble(a[e])/realToDouble(b[e])) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      FADd B(b[e].data.value,b[e].data.grad[i]) ;
      // Check derivative within tight tolerance
      CHECK(TOLCHECK(c[e].data.grad[i],(A / B).grad,tol)) ;
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = a[e] ;
    c[e] += b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == realToDouble(a[e])+realToDouble(b[e])) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      FADd B(b[e].data.value,b[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(A + B).grad,tol)) ;
    }
  }

  for(unsigned int e = 0; e < N; ++e) {
    c[e] = a[e] ;
    c[e] -= b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == realToDouble(a[e])-realToDouble(b[e])) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      FADd B(b[e].data.value,b[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(A - B).grad,tol)) ;
    }
  }

  for(unsigned int e = 0; e < N; ++e) {
    c[e] = a[e] ;
    c[e] *= b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == realToDouble(a[e])*realToDouble(b[e])) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      FADd B(b[e].data.value,b[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(A * B).grad,tol)) ;
    }
  }

  for(unsigned int e = 0; e < N; ++e) {
    c[e] = a[e] ;
    c[e] /= b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == realToDouble(a[e])/realToDouble(b[e])) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      FADd B(b[e].data.value,b[e].data.grad[i]) ;
      // Check derivative within tight tolerance
      CHECK(TOLCHECK(c[e].data.grad[i],(A / B).grad,tol)) ;
    }
  }

  // Checking with a as a double
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = ad[e]+b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == realToDouble(ad[e])+realToDouble(b[e])) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A = ad[e] ;
      FADd B(b[e].data.value,b[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(A + B).grad,tol)) ;
    }
  }

  for(unsigned int e = 0; e < N; ++e) {
    c[e] = ad[e]-b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == realToDouble(ad[e])-realToDouble(b[e])) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A = ad[e] ;
      FADd B(b[e].data.value,b[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(A - B).grad,tol)) ;
    }
  }

  for(unsigned int e = 0; e < N; ++e) {
    c[e] = ad[e]*b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == realToDouble(ad[e])*realToDouble(b[e])) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A= ad[e] ;
      FADd B(b[e].data.value,b[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(A * B).grad,tol)) ;
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = ad[e]/b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == realToDouble(ad[e])/realToDouble(b[e])) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A = ad[e] ;
      FADd B(b[e].data.value,b[e].data.grad[i]) ;
      // Check derivative within tight tolerance
      CHECK(TOLCHECK(c[e].data.grad[i], (A / B).grad,tol)) ;
    }
  }

  // Check b is double
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = a[e]+bd[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == realToDouble(a[e])+realToDouble(bd[e])) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      FADd B=bd[e] ;
      CHECK(TOLCHECK(c[e].data.grad[i],(A + B).grad,tol)) ;
    }
  }

  for(unsigned int e = 0; e < N; ++e) {
    c[e] = a[e]-bd[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == realToDouble(a[e])-realToDouble(bd[e])) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      FADd B=bd[e] ;
      CHECK(TOLCHECK(c[e].data.grad[i],(A - B).grad,tol)) ;
    }
  }

  for(unsigned int e = 0; e < N; ++e) {
    c[e] = a[e]*bd[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == realToDouble(a[e])*realToDouble(bd[e])) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      FADd B=bd[e] ;
      CHECK(TOLCHECK(c[e].data.grad[i],(A * B).grad,tol)) ;
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = a[e]/bd[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == realToDouble(a[e])/realToDouble(bd[e])) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      FADd B=bd[e] ;
      // Check derivative within tight tolerance
      CHECK(TOLCHECK(c[e].data.grad[i],(A / B).grad,tol)) ;
    }
  }

  // Unary operator checks
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = -a[e] ;
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == -realToDouble(a[e])) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(-A).grad,tol)) ;
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = +a[e] ;
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == +realToDouble(a[e])) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(+A).grad,tol)) ;
    }
  }
  
}

TEST_CASE("VFAD arithmetic functions") {
  std::random_device rd;
  std::default_random_engine re(rd());
  std::uniform_real_distribution<double> rdist(-10, 10);
  
  typedef VFAD type;
  
  unsigned int const N = 100;
  const double tol = 1e-6 ;
  
  std::vector<type> a(N), b(N), c(N);
  std::vector<double> ad(N), bd(N) ;
  for(unsigned int e = 0; e < N; ++e) {
    a[e] = rdist(re);
    b[e] = rdist(re);
    if(b[e].data.value == 0) // Don't let b be zero so it can test a/b
      b[e].data.value = 1.0 ;
    ad[e] = realToDouble(a[e]) ;
    bd[e] = realToDouble(b[e]) ;
    for(unsigned int i=0;i<type::maxN;++i) {
      a[e].data.grad[i] = rdist(re) ;
      b[e].data.grad[i] = rdist(re) ;
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = sqrt(abs(a[e])) ;
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == sqrt(abs(realToDouble(a[e])))) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(sqrt(abs(A))).grad,tol)) ;
    }
  }

  for(unsigned int e = 0; e < N; ++e) {
    c[e] = erf(a[e]) ;
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == erf(realToDouble(a[e]))) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(erf(A)).grad,tol)) ;
    }
  }

  for(unsigned int e = 0; e < N; ++e) {
    c[e] = sin(a[e]) ;
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == sin(realToDouble(a[e]))) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(sin(A)).grad,tol)) ;
    }
  }
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = cos(a[e]) ;
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == cos(realToDouble(a[e]))) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(cos(A)).grad,tol)) ;
    }
  }

  for(unsigned int e = 0; e < N; ++e) {
    c[e] = tan(a[e]) ;
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == tan(realToDouble(a[e]))) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(tan(A)).grad,tol)) ;
    }
  }

  for(unsigned int e = 0; e < N; ++e) {
    c[e] = exp(a[e]) ;
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == exp(realToDouble(a[e]))) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(exp(A)).grad,tol)) ;
    }
  }

  for(unsigned int e = 0; e < N; ++e) {
    c[e] = log(fabs(b[e])) ;
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == log(fabs(realToDouble(b[e])))) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(b[e].data.value,b[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(log(fabs(A))).grad,tol)) ;
    }
  }
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = log10(fabs(b[e])) ;
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(TOLCHECK(realToDouble(c[e]),log10(fabs(realToDouble(b[e]))),tol)) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(b[e].data.value,b[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(log10(fabs(A))).grad,tol)) ;
    }
  }

  for(unsigned int e = 0; e < N; ++e) {
    c[e] = pow(a[e],3) ;
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == pow(realToDouble(a[e]),3)) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(pow(A,3)).grad,tol)) ;
    }
  }
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = pow(a[e],3.0f) ;
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == pow(realToDouble(a[e]),3.0f)) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(pow(A,3.0f)).grad,tol)) ;
    }
  }
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = pow(a[e],3.0) ;
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == pow(realToDouble(a[e]),3.0)) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(pow(A,3.0)).grad,tol)) ;
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = pow(abs(b[e]),a[e]) ;
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == pow(abs(realToDouble(b[e])),realToDouble(a[e]))) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      FADd B(b[e].data.value,b[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(pow(abs(B),A)).grad,tol)) ;
    }
  }

  for(unsigned int e = 0; e < N; ++e) {
    c[e] = pow(2,a[e]) ;
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == pow(2,realToDouble(a[e]))) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(pow(2,A)).grad,tol)) ;
    }
  }

  for(unsigned int e = 0; e < N; ++e) {
    c[e] = pow(2.0,a[e]) ;
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == pow(2.0,realToDouble(a[e]))) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(pow(2.0,A)).grad,tol)) ;
    }
  }
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = pow(2.0f,a[e]) ;
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == pow(2.0f,realToDouble(a[e]))) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(pow(2.0f,A)).grad,tol)) ;
    }
  }
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = sinh(a[e]) ;
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == sinh(realToDouble(a[e]))) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(sinh(A)).grad,tol)) ;
    }
  }

  for(unsigned int e = 0; e < N; ++e) {
    c[e] = cosh(a[e]) ;
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == cosh(realToDouble(a[e]))) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(cosh(A)).grad,tol)) ;
    }
  }

  for(unsigned int e = 0; e < N; ++e) {
    c[e] = tanh(a[e]) ;
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == tanh(realToDouble(a[e]))) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(tanh(A)).grad,tol)) ;
    }
  }

  for(unsigned int e = 0; e < N; ++e) {
    c[e] = asin(a[e]*0.1) ;
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == asin(realToDouble(a[e])*0.1)) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(asin(A*0.1)).grad,tol)) ;
    }
  }
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = acos(a[e]*0.1) ;
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == acos(realToDouble(a[e])*0.1)) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(acos(A*0.1)).grad,tol)) ;
    }
  }
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = atan(a[e]*0.1) ;
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == atan(realToDouble(a[e])*0.1)) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(atan(A*0.1)).grad,tol)) ;
    }
  }

  for(unsigned int e = 0; e < N; ++e) {
    c[e] = asinh(a[e]) ;
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == asinh(realToDouble(a[e]))) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(asinh(A)).grad,tol)) ;
    }
  }
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = acosh(1.0+fabs(b[e])) ;
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == acosh(1.0+fabs(realToDouble(b[e])))) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd B(b[e].data.value,b[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(acosh(1.0+fabs(B))).grad,tol)) ;
    }
  }
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = atanh(a[e]*0.1) ;
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    CHECK(realToDouble(c[e]) == atanh(realToDouble(a[e])*0.1)) ;
    for(unsigned int i = 0; i < type::maxN; ++i) {
      FADd A(a[e].data.value,a[e].data.grad[i]) ;
      CHECK(TOLCHECK(c[e].data.grad[i],(atanh(A*0.1)).grad,tol)) ;
    }
  }
}
