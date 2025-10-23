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

TEST_CASE("vtype<int, 8> arithmetic operations") {
  std::random_device rd;
  std::default_random_engine re(rd());
  std::uniform_int_distribution<int> rdist(-10, 10);
  
  typedef vtype<int, 8> type;
  
  unsigned int const N = 100;
  
  std::vector<type> a(N), b(N), c(N);
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      a[e][i] = rdist(re);
      b[e][i] = rdist(re);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = a[e]+b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i] + b[e][i]);
    }
  }

  for(unsigned int e = 0; e < N; ++e) {
    a[e] += b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i]);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = c[e] - b[e];
    a[e] -= b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i]);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = a[e]*b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i]*b[e][i]);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    a[e] *= b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i]);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = c[e]/b[e];
    a[e] /= b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i]);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = a[e]*type(10);
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i]*10);
    }
  }
}

TEST_CASE("vtype<float, 4> arithmetic operations") {
  std::random_device rd;
  std::default_random_engine re(rd());
  std::uniform_real_distribution<float> rdist(-10.0, 10.0);
  
  typedef vtype<float, 4> type;
  
  unsigned int const N = 100;
  
  std::vector<type> a(N), b(N), c(N);
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      a[e][i] = rdist(re);
      b[e][i] = rdist(re);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = a[e]+b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i] + b[e][i]);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    a[e] += b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i]);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = c[e] - b[e];
    a[e] -= b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i]);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = a[e]*b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i]*b[e][i]);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    a[e] *= b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i]);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = c[e]/b[e];
    a[e] /= b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i]);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = a[e]*type(10.0f);
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i]*10.0f);
    }
  }
}

TEST_CASE("vtype<float, 8> arithmetic operations") {
  std::random_device rd;
  std::default_random_engine re(rd());
  std::uniform_real_distribution<float> rdist(-10.0, 10.0);
  
  typedef vtype<float, 8> type;
  
  unsigned int const N = 100;
  
  std::vector<type> a(N), b(N), c(N);
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      a[e][i] = rdist(re);
      b[e][i] = rdist(re);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = a[e]+b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i] + b[e][i]);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    a[e] += b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i]);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = c[e] - b[e];
    a[e] -= b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i]);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = a[e]*b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i]*b[e][i]);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    a[e] *= b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i]);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = c[e]/b[e];
    a[e] /= b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i]);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = a[e]*type(10.0f);
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i]*10.0f);
    }
  }
}

TEST_CASE("vtype<double, 4> arithmetic operations") {
  std::random_device rd;
  std::default_random_engine re(rd());
  std::uniform_real_distribution<double> rdist(-10.0, 10.0);
  
  typedef vtype<double, 4> type;
  
  unsigned int const N = 100;
  
  std::vector<type> a(N), b(N), c(N);
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      a[e][i] = rdist(re);
      b[e][i] = rdist(re);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = a[e]+b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i] + b[e][i]);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    a[e] += b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i]);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = c[e] - b[e];
    a[e] -= b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i]);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = a[e]*b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i]*b[e][i]);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    a[e] *= b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i]);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = c[e]/b[e];
    a[e] /= b[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i]);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = a[e]*type(10.0);
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == a[e][i]*10.0);
    }
  }
  // Unary negation
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = -a[e];
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == -a[e][i]);
    }
  }
  // fabs
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = fabs(a[e]);
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == fabs(a[e][i]));
    }
  }
  // abs
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = abs(a[e]);
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == abs(a[e][i]));
    }
  }
  // max
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = max(a[e],b[e]);
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == max(a[e][i],b[e][i]));
    }
  }
  // min
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = min(a[e],b[e]);
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == min(a[e][i],b[e][i]));
    }
  }

  // select
  for(unsigned int e = 0; e < N; ++e) {
    c[e] = select(a[e]<b[e],a[e],b[e]);
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(c[e][i] == (a[e][i]<b[e][i]?a[e][i]:b[e][i]));
    }
  }
}


TEST_CASE("vtype<float, 4> floor, ceil, round operations") {
  std::random_device rd;
  std::default_random_engine re(rd());
  std::uniform_real_distribution<float> rdist(-10.0, 10.0);
  
  typedef vtype<float, 4> type;
  
  unsigned int const N = 100;
  
  std::vector<type> a(N), b(N);
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      a[e][i] = rdist(re);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    b[e] = ceil(a[e]);
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(b[e][i] == std::ceil(a[e][i]));
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    b[e] = floor(a[e]);
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(b[e][i] == std::floor(a[e][i]));
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    b[e] = round(a[e]);
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(b[e][i] == std::round(a[e][i]));
    }
  }
}

TEST_CASE("vtype<float, 8> floor, ceil, round operations") {
  std::random_device rd;
  std::default_random_engine re(rd());
  std::uniform_real_distribution<float> rdist(-10.0, 10.0);
  
  typedef vtype<float, 8> type;
  
  unsigned int const N = 100;
  
  std::vector<type> a(N), b(N);
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      a[e][i] = rdist(re);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    b[e] = ceil(a[e]);
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(b[e][i] == std::ceil(a[e][i]));
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    b[e] = floor(a[e]);
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(b[e][i] == std::floor(a[e][i]));
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    b[e] = round(a[e]);
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(b[e][i] == std::round(a[e][i]));
    }
  }
}

TEST_CASE("vtype<double, 4> floor, ceil, round operations") {
  std::random_device rd;
  std::default_random_engine re(rd());
  std::uniform_real_distribution<double> rdist(-10.0, 10.0);
  
  typedef vtype<double, 4> type;
  
  unsigned int const N = 100;
  
  std::vector<type> a(N), b(N);
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      a[e][i] = rdist(re);
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    b[e] = ceil(a[e]);
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(b[e][i] == std::ceil(a[e][i]));
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    b[e] = floor(a[e]);
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(b[e][i] == std::floor(a[e][i]));
    }
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    b[e] = round(a[e]);
  }
  
  for(unsigned int e = 0; e < N; ++e) {
    for(unsigned int i = 0; i < type::W; ++i) {
      CHECK(b[e][i] == std::round(a[e][i]));
    }
  }
}

