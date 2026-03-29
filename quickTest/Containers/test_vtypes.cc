#include <Tools/vtypes.h>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <random>
#include <array>

using namespace Loci;
using namespace std;

template<typename M>
struct MaskTestFixture {
  size_t size;
  vector<typename M::bitmask_t> lhs_bitmask;
  vector<typename M::bitmask_t> rhs_bitmask;

  MaskTestFixture();
};

template<typename M>
MaskTestFixture<M>::MaskTestFixture() {
  mt19937 rng(23764);
  uniform_int_distribution<unsigned int> dist(0, (1u << M::W)-1);
  size = 100;
  lhs_bitmask.resize(size);
  rhs_bitmask.resize(size);
  for(size_t i = 0; i < size; ++i) {
    lhs_bitmask[i] = dist(rng);
    rhs_bitmask[i] = dist(rng);
  }
}

#define MASK_TYPES vmask<4, 4>, vmask<4, 8>, vmask<4, 16>, vmask<8, 4>, vmask<8, 8>, vmask<8, 16>

TEST_CASE_TEMPLATE("Check constructor vmask(bitmask)", M, MASK_TYPES) {
  MaskTestFixture<M> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    M m(fix.lhs_bitmask[i]);
    CHECK_EQ((typename M::bitmask_t)m, fix.lhs_bitmask[i]);
  }
}

TEST_CASE_TEMPLATE("Check m1 &= m2", M, MASK_TYPES) {
  MaskTestFixture<M> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    M m1(fix.lhs_bitmask[i]);
    M m2(fix.rhs_bitmask[i]);

    typename M::bitmask_t r_bitmask = fix.lhs_bitmask[i] & fix.rhs_bitmask[i];

    m1 &= m2;

    for(unsigned int i = 0; i < M::W; ++i) {
      CHECK_EQ(m1[i], (bool)(r_bitmask & (1u << i)));
    }

    CHECK_EQ((typename M::bitmask_t)m1, r_bitmask);
  }
}

TEST_CASE_TEMPLATE("Check m1 & m2 and m1 && m2", M, MASK_TYPES) {
  MaskTestFixture<M> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    M m1(fix.lhs_bitmask[i]);
    M m2(fix.rhs_bitmask[i]);

    typename M::bitmask_t r_bitmask = fix.lhs_bitmask[i] & fix.rhs_bitmask[i];

    M r1 = m1 & m2;
    M r2 = m1 && m2;

    for(unsigned int i = 0; i < M::W; ++i) {
      CHECK_EQ(r1[i], (bool)(r_bitmask & (1u << i)));
    }

    CHECK_EQ((typename M::bitmask_t)r1, r_bitmask);

    for(unsigned int i = 0; i < M::W; ++i) {
      CHECK_EQ(r2[i], (bool)(r_bitmask & (1u << i)));
    }

    CHECK_EQ((typename M::bitmask_t)r2, r_bitmask);
  }
}

TEST_CASE_TEMPLATE("Check m1 |= m2", M, MASK_TYPES) {
  MaskTestFixture<M> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    M m1(fix.lhs_bitmask[i]);
    M m2(fix.rhs_bitmask[i]);

    typename M::bitmask_t r_bitmask = fix.lhs_bitmask[i] | fix.rhs_bitmask[i];

    m1 |= m2;

    for(unsigned int i = 0; i < M::W; ++i) {
      CHECK_EQ(m1[i], (bool)(r_bitmask & (1u << i)));
    }

    CHECK_EQ((typename M::bitmask_t)m1, r_bitmask);
  }
}

TEST_CASE_TEMPLATE("Check m1 | m2 and m1 || m2", M, MASK_TYPES) {
  MaskTestFixture<M> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    M m1(fix.lhs_bitmask[i]);
    M m2(fix.rhs_bitmask[i]);

    typename M::bitmask_t r_bitmask = fix.lhs_bitmask[i] | fix.rhs_bitmask[i];

    M r1 = m1 | m2;
    M r2 = m1 || m2;

    for(unsigned int i = 0; i < M::W; ++i) {
      CHECK_EQ(r1[i], (bool)(r_bitmask & (1u << i)));
    }

    CHECK_EQ((typename M::bitmask_t)r1, r_bitmask);

    for(unsigned int i = 0; i < M::W; ++i) {
      CHECK_EQ(r2[i], (bool)(r_bitmask & (1u << i)));
    }

    CHECK_EQ((typename M::bitmask_t)r2, r_bitmask);
  }
}

TEST_CASE_TEMPLATE("Check m1 ^= m2", M, MASK_TYPES) {
  MaskTestFixture<M> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    M m1(fix.lhs_bitmask[i]);
    M m2(fix.rhs_bitmask[i]);

    typename M::bitmask_t r_bitmask = fix.lhs_bitmask[i] ^ fix.rhs_bitmask[i];

    m1 ^= m2;

    for(unsigned int i = 0; i < M::W; ++i) {
      CHECK_EQ(m1[i], (bool)(r_bitmask & (1u << i)));
    }

    CHECK_EQ((typename M::bitmask_t)m1, r_bitmask);
  }
}

TEST_CASE_TEMPLATE("Check m1 ^ m2", M, MASK_TYPES) {
  MaskTestFixture<M> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    M m1(fix.lhs_bitmask[i]);
    M m2(fix.rhs_bitmask[i]);

    typename M::bitmask_t r_bitmask = fix.lhs_bitmask[i] ^ fix.rhs_bitmask[i];

    M r = m1 ^ m2;

    for(unsigned int i = 0; i < M::W; ++i) {
      CHECK_EQ(r[i], (bool)(r_bitmask & (1u << i)));
    }

    CHECK_EQ((typename M::bitmask_t)r, r_bitmask);
  }
}

TEST_SUITE_BEGIN("vmask<4, 4>");
TEST_CASE("ctor") {
  vmask<4, 4> mask_ctor_bitmask(0b1101);
  CHECK_EQ(mask_ctor_bitmask[0], true);
  CHECK_EQ(mask_ctor_bitmask[1], false);
  CHECK_EQ(mask_ctor_bitmask[2], true);
  CHECK_EQ(mask_ctor_bitmask[3], true);

  vmask<4, 4> mask_ctor_bools(true, false, false, true);
  CHECK_EQ(mask_ctor_bools[0], true);
  CHECK_EQ(mask_ctor_bools[1], false);
  CHECK_EQ(mask_ctor_bools[2], false);
  CHECK_EQ(mask_ctor_bools[3], true);

  vmask<4, 4> mask_copy_ctor(mask_ctor_bitmask);
  CHECK_EQ(mask_copy_ctor[0], true);
  CHECK_EQ(mask_copy_ctor[1], false);
  CHECK_EQ(mask_copy_ctor[2], true);
  CHECK_EQ(mask_copy_ctor[3], true);

  vmask<4, 4> mask_copy_assignment;
  mask_copy_assignment = mask_ctor_bools;
  CHECK_EQ(mask_copy_assignment[0], true);
  CHECK_EQ(mask_copy_assignment[1], false);
  CHECK_EQ(mask_copy_assignment[2], false);
  CHECK_EQ(mask_copy_assignment[3], true);
}
TEST_SUITE_END();

TEST_SUITE_BEGIN("vmask<4, 8>");
TEST_CASE("ctor") {
  vmask<4, 8> mask_ctor_bitmask(0b10010101);
  CHECK_EQ(mask_ctor_bitmask[0], true);
  CHECK_EQ(mask_ctor_bitmask[1], false);
  CHECK_EQ(mask_ctor_bitmask[2], true);
  CHECK_EQ(mask_ctor_bitmask[3], false);
  CHECK_EQ(mask_ctor_bitmask[4], true);
  CHECK_EQ(mask_ctor_bitmask[5], false);
  CHECK_EQ(mask_ctor_bitmask[6], false);
  CHECK_EQ(mask_ctor_bitmask[7], true);

  vmask<4, 8> mask_ctor_bools(false, false, true, false, true, true, true, false);
  CHECK_EQ(mask_ctor_bools[0], false);
  CHECK_EQ(mask_ctor_bools[1], false);
  CHECK_EQ(mask_ctor_bools[2], true);
  CHECK_EQ(mask_ctor_bools[3], false);
  CHECK_EQ(mask_ctor_bools[4], true);
  CHECK_EQ(mask_ctor_bools[5], true);
  CHECK_EQ(mask_ctor_bools[6], true);
  CHECK_EQ(mask_ctor_bools[7], false);

  vmask<4, 8> mask_copy_ctor(mask_ctor_bitmask);
  CHECK_EQ(mask_copy_ctor[0], true);
  CHECK_EQ(mask_copy_ctor[1], false);
  CHECK_EQ(mask_copy_ctor[2], true);
  CHECK_EQ(mask_copy_ctor[3], false);
  CHECK_EQ(mask_copy_ctor[4], true);
  CHECK_EQ(mask_copy_ctor[5], false);
  CHECK_EQ(mask_copy_ctor[6], false);
  CHECK_EQ(mask_copy_ctor[7], true);

  vmask<4, 8> mask_copy_assignment;
  mask_copy_assignment = mask_ctor_bools;
  CHECK_EQ(mask_copy_assignment[0], false);
  CHECK_EQ(mask_copy_assignment[1], false);
  CHECK_EQ(mask_copy_assignment[2], true);
  CHECK_EQ(mask_copy_assignment[3], false);
  CHECK_EQ(mask_copy_assignment[4], true);
  CHECK_EQ(mask_copy_assignment[5], true);
  CHECK_EQ(mask_copy_assignment[6], true);
  CHECK_EQ(mask_copy_assignment[7], false);
}
TEST_SUITE_END();

TEST_SUITE_BEGIN("vmask<4, 16>");
TEST_CASE("ctor") {
  vmask<4, 16> mask_ctor_bitmask(0b0100100111001010);
  CHECK_EQ(mask_ctor_bitmask[0], false);
  CHECK_EQ(mask_ctor_bitmask[1], true);
  CHECK_EQ(mask_ctor_bitmask[2], false);
  CHECK_EQ(mask_ctor_bitmask[3], true);
  CHECK_EQ(mask_ctor_bitmask[4], false);
  CHECK_EQ(mask_ctor_bitmask[5], false);
  CHECK_EQ(mask_ctor_bitmask[6], true);
  CHECK_EQ(mask_ctor_bitmask[7], true);
  CHECK_EQ(mask_ctor_bitmask[8], true);
  CHECK_EQ(mask_ctor_bitmask[9], false);
  CHECK_EQ(mask_ctor_bitmask[10], false);
  CHECK_EQ(mask_ctor_bitmask[11], true);
  CHECK_EQ(mask_ctor_bitmask[12], false);
  CHECK_EQ(mask_ctor_bitmask[13], false);
  CHECK_EQ(mask_ctor_bitmask[14], true);
  CHECK_EQ(mask_ctor_bitmask[15], false);

  vmask<4, 16> mask_ctor_bools(
    true, true, false, false, false, true, true, false,
    false, false, true, true, true, false, true, true
  );
  CHECK_EQ(mask_ctor_bools[0], true);
  CHECK_EQ(mask_ctor_bools[1], true);
  CHECK_EQ(mask_ctor_bools[2], false);
  CHECK_EQ(mask_ctor_bools[3], false);
  CHECK_EQ(mask_ctor_bools[4], false);
  CHECK_EQ(mask_ctor_bools[5], true);
  CHECK_EQ(mask_ctor_bools[6], true);
  CHECK_EQ(mask_ctor_bools[7], false);
  CHECK_EQ(mask_ctor_bools[8], false);
  CHECK_EQ(mask_ctor_bools[9], false);
  CHECK_EQ(mask_ctor_bools[10], true);
  CHECK_EQ(mask_ctor_bools[11], true);
  CHECK_EQ(mask_ctor_bools[12], true);
  CHECK_EQ(mask_ctor_bools[13], false);
  CHECK_EQ(mask_ctor_bools[14], true);
  CHECK_EQ(mask_ctor_bools[15], true);

  vmask<4, 16> mask_copy_ctor(mask_ctor_bitmask);
  CHECK_EQ(mask_copy_ctor[0], false);
  CHECK_EQ(mask_copy_ctor[1], true);
  CHECK_EQ(mask_copy_ctor[2], false);
  CHECK_EQ(mask_copy_ctor[3], true);
  CHECK_EQ(mask_copy_ctor[4], false);
  CHECK_EQ(mask_copy_ctor[5], false);
  CHECK_EQ(mask_copy_ctor[6], true);
  CHECK_EQ(mask_copy_ctor[7], true);
  CHECK_EQ(mask_copy_ctor[8], true);
  CHECK_EQ(mask_copy_ctor[9], false);
  CHECK_EQ(mask_copy_ctor[10], false);
  CHECK_EQ(mask_copy_ctor[11], true);
  CHECK_EQ(mask_copy_ctor[12], false);
  CHECK_EQ(mask_copy_ctor[13], false);
  CHECK_EQ(mask_copy_ctor[14], true);
  CHECK_EQ(mask_copy_ctor[15], false);

  vmask<4, 16> mask_copy_assignment;
  mask_copy_assignment = mask_ctor_bools;
  CHECK_EQ(mask_copy_assignment[0], true);
  CHECK_EQ(mask_copy_assignment[1], true);
  CHECK_EQ(mask_copy_assignment[2], false);
  CHECK_EQ(mask_copy_assignment[3], false);
  CHECK_EQ(mask_copy_assignment[4], false);
  CHECK_EQ(mask_copy_assignment[5], true);
  CHECK_EQ(mask_copy_assignment[6], true);
  CHECK_EQ(mask_copy_assignment[7], false);
  CHECK_EQ(mask_copy_assignment[8], false);
  CHECK_EQ(mask_copy_assignment[9], false);
  CHECK_EQ(mask_copy_assignment[10], true);
  CHECK_EQ(mask_copy_assignment[11], true);
  CHECK_EQ(mask_copy_assignment[12], true);
  CHECK_EQ(mask_copy_assignment[13], false);
  CHECK_EQ(mask_copy_assignment[14], true);
  CHECK_EQ(mask_copy_assignment[15], true);
}
TEST_SUITE_END();

TEST_SUITE_BEGIN("vmask<8, 4>");
TEST_CASE("ctor") {
  vmask<8, 4> mask_ctor_bitmask(0b1101);
  CHECK_EQ(mask_ctor_bitmask[0], true);
  CHECK_EQ(mask_ctor_bitmask[1], false);
  CHECK_EQ(mask_ctor_bitmask[2], true);
  CHECK_EQ(mask_ctor_bitmask[3], true);

  vmask<8, 4> mask_ctor_bools(true, false, false, true);
  CHECK_EQ(mask_ctor_bools[0], true);
  CHECK_EQ(mask_ctor_bools[1], false);
  CHECK_EQ(mask_ctor_bools[2], false);
  CHECK_EQ(mask_ctor_bools[3], true);

  vmask<8, 4> mask_copy_ctor(mask_ctor_bitmask);
  CHECK_EQ(mask_copy_ctor[0], true);
  CHECK_EQ(mask_copy_ctor[1], false);
  CHECK_EQ(mask_copy_ctor[2], true);
  CHECK_EQ(mask_copy_ctor[3], true);

  vmask<8, 4> mask_copy_assignment;
  mask_copy_assignment = mask_ctor_bools;
  CHECK_EQ(mask_copy_assignment[0], true);
  CHECK_EQ(mask_copy_assignment[1], false);
  CHECK_EQ(mask_copy_assignment[2], false);
  CHECK_EQ(mask_copy_assignment[3], true);
}
TEST_SUITE_END();

TEST_SUITE_BEGIN("vmask<8, 8>");
TEST_CASE("ctor") {
  vmask<8, 8> mask_ctor_bitmask(0b10010101);
  CHECK_EQ(mask_ctor_bitmask[0], true);
  CHECK_EQ(mask_ctor_bitmask[1], false);
  CHECK_EQ(mask_ctor_bitmask[2], true);
  CHECK_EQ(mask_ctor_bitmask[3], false);
  CHECK_EQ(mask_ctor_bitmask[4], true);
  CHECK_EQ(mask_ctor_bitmask[5], false);
  CHECK_EQ(mask_ctor_bitmask[6], false);
  CHECK_EQ(mask_ctor_bitmask[7], true);

  vmask<8, 8> mask_ctor_bools(false, false, true, false, true, true, true, false);
  CHECK_EQ(mask_ctor_bools[0], false);
  CHECK_EQ(mask_ctor_bools[1], false);
  CHECK_EQ(mask_ctor_bools[2], true);
  CHECK_EQ(mask_ctor_bools[3], false);
  CHECK_EQ(mask_ctor_bools[4], true);
  CHECK_EQ(mask_ctor_bools[5], true);
  CHECK_EQ(mask_ctor_bools[6], true);
  CHECK_EQ(mask_ctor_bools[7], false);

  vmask<8, 8> mask_copy_ctor(mask_ctor_bitmask);
  CHECK_EQ(mask_copy_ctor[0], true);
  CHECK_EQ(mask_copy_ctor[1], false);
  CHECK_EQ(mask_copy_ctor[2], true);
  CHECK_EQ(mask_copy_ctor[3], false);
  CHECK_EQ(mask_copy_ctor[4], true);
  CHECK_EQ(mask_copy_ctor[5], false);
  CHECK_EQ(mask_copy_ctor[6], false);
  CHECK_EQ(mask_copy_ctor[7], true);

  vmask<8, 8> mask_copy_assignment;
  mask_copy_assignment = mask_ctor_bools;
  CHECK_EQ(mask_copy_assignment[0], false);
  CHECK_EQ(mask_copy_assignment[1], false);
  CHECK_EQ(mask_copy_assignment[2], true);
  CHECK_EQ(mask_copy_assignment[3], false);
  CHECK_EQ(mask_copy_assignment[4], true);
  CHECK_EQ(mask_copy_assignment[5], true);
  CHECK_EQ(mask_copy_assignment[6], true);
  CHECK_EQ(mask_copy_assignment[7], false);
}
TEST_SUITE_END();

TEST_SUITE_BEGIN("vmask<8, 16>");
TEST_CASE("ctor") {
  vmask<8, 16> mask_ctor_bitmask(0b0100100111001010);
  CHECK_EQ(mask_ctor_bitmask[0], false);
  CHECK_EQ(mask_ctor_bitmask[1], true);
  CHECK_EQ(mask_ctor_bitmask[2], false);
  CHECK_EQ(mask_ctor_bitmask[3], true);
  CHECK_EQ(mask_ctor_bitmask[4], false);
  CHECK_EQ(mask_ctor_bitmask[5], false);
  CHECK_EQ(mask_ctor_bitmask[6], true);
  CHECK_EQ(mask_ctor_bitmask[7], true);
  CHECK_EQ(mask_ctor_bitmask[8], true);
  CHECK_EQ(mask_ctor_bitmask[9], false);
  CHECK_EQ(mask_ctor_bitmask[10], false);
  CHECK_EQ(mask_ctor_bitmask[11], true);
  CHECK_EQ(mask_ctor_bitmask[12], false);
  CHECK_EQ(mask_ctor_bitmask[13], false);
  CHECK_EQ(mask_ctor_bitmask[14], true);
  CHECK_EQ(mask_ctor_bitmask[15], false);

  vmask<8, 16> mask_ctor_bools(
    true, true, false, false, false, true, true, false,
    false, false, true, true, true, false, true, true
  );
  CHECK_EQ(mask_ctor_bools[0], true);
  CHECK_EQ(mask_ctor_bools[1], true);
  CHECK_EQ(mask_ctor_bools[2], false);
  CHECK_EQ(mask_ctor_bools[3], false);
  CHECK_EQ(mask_ctor_bools[4], false);
  CHECK_EQ(mask_ctor_bools[5], true);
  CHECK_EQ(mask_ctor_bools[6], true);
  CHECK_EQ(mask_ctor_bools[7], false);
  CHECK_EQ(mask_ctor_bools[8], false);
  CHECK_EQ(mask_ctor_bools[9], false);
  CHECK_EQ(mask_ctor_bools[10], true);
  CHECK_EQ(mask_ctor_bools[11], true);
  CHECK_EQ(mask_ctor_bools[12], true);
  CHECK_EQ(mask_ctor_bools[13], false);
  CHECK_EQ(mask_ctor_bools[14], true);
  CHECK_EQ(mask_ctor_bools[15], true);

  vmask<8, 16> mask_copy_ctor(mask_ctor_bitmask);
  CHECK_EQ(mask_copy_ctor[0], false);
  CHECK_EQ(mask_copy_ctor[1], true);
  CHECK_EQ(mask_copy_ctor[2], false);
  CHECK_EQ(mask_copy_ctor[3], true);
  CHECK_EQ(mask_copy_ctor[4], false);
  CHECK_EQ(mask_copy_ctor[5], false);
  CHECK_EQ(mask_copy_ctor[6], true);
  CHECK_EQ(mask_copy_ctor[7], true);
  CHECK_EQ(mask_copy_ctor[8], true);
  CHECK_EQ(mask_copy_ctor[9], false);
  CHECK_EQ(mask_copy_ctor[10], false);
  CHECK_EQ(mask_copy_ctor[11], true);
  CHECK_EQ(mask_copy_ctor[12], false);
  CHECK_EQ(mask_copy_ctor[13], false);
  CHECK_EQ(mask_copy_ctor[14], true);
  CHECK_EQ(mask_copy_ctor[15], false);

  vmask<8, 16> mask_copy_assignment;
  mask_copy_assignment = mask_ctor_bools;
  CHECK_EQ(mask_copy_assignment[0], true);
  CHECK_EQ(mask_copy_assignment[1], true);
  CHECK_EQ(mask_copy_assignment[2], false);
  CHECK_EQ(mask_copy_assignment[3], false);
  CHECK_EQ(mask_copy_assignment[4], false);
  CHECK_EQ(mask_copy_assignment[5], true);
  CHECK_EQ(mask_copy_assignment[6], true);
  CHECK_EQ(mask_copy_assignment[7], false);
  CHECK_EQ(mask_copy_assignment[8], false);
  CHECK_EQ(mask_copy_assignment[9], false);
  CHECK_EQ(mask_copy_assignment[10], true);
  CHECK_EQ(mask_copy_assignment[11], true);
  CHECK_EQ(mask_copy_assignment[12], true);
  CHECK_EQ(mask_copy_assignment[13], false);
  CHECK_EQ(mask_copy_assignment[14], true);
  CHECK_EQ(mask_copy_assignment[15], true);
}
TEST_SUITE_END();

template<typename V>
struct VTypeRealTestFixture {
  size_t size;
  vector<std::array<typename V::scalar_t,V::W>> lhs;
  vector<std::array<typename V::scalar_t,V::W>> rhs;

  VTypeRealTestFixture();
};

template<typename V>
VTypeRealTestFixture<V>::VTypeRealTestFixture() {
  mt19937 rng(23764);
  uniform_real_distribution<typename V::scalar_t> dist(-100, 100);
  size = 100;
  lhs.resize(size);
  rhs.resize(size);
  for(size_t i = 0; i < size; ++i) {
    for(unsigned int j = 0; j < V::W; ++j) {
      lhs[i][j] = dist(rng);
      rhs[i][j] = dist(rng);
    }
  }
}

#define VECTOR_REAL_TYPES vtype<float, 4>, vtype<float, 8>, vtype<float, 16>, vtype<double, 4>, vtype<double, 8>, vtype<double, 16>

TEST_CASE_TEMPLATE("Check constructor vtype(scalar_t const * mem)", V, VECTOR_REAL_TYPES) {
  VTypeRealTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v(&fix.lhs[i][0]);

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v[j], fix.lhs[i][j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check v1 += v2", V, VECTOR_REAL_TYPES) {
  VTypeRealTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);
    V v2(&fix.rhs[i][0]);

    v1 += v2;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v1[j], fix.lhs[i][j] + fix.rhs[i][j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check v1 + v2", V, VECTOR_REAL_TYPES) {
  VTypeRealTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);
    V v2(&fix.rhs[i][0]);

    V v3 = v1 + v2;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v3[j], fix.lhs[i][j] + fix.rhs[i][j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check v1 -= v2", V, VECTOR_REAL_TYPES) {
  VTypeRealTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);
    V v2(&fix.rhs[i][0]);

    v1 -= v2;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v1[j], fix.lhs[i][j] - fix.rhs[i][j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check v1 - v2", V, VECTOR_REAL_TYPES) {
  VTypeRealTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);
    V v2(&fix.rhs[i][0]);

    V v3 = v1 - v2;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v3[j], fix.lhs[i][j] - fix.rhs[i][j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check v1 *= v2", V, VECTOR_REAL_TYPES) {
  VTypeRealTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);
    V v2(&fix.rhs[i][0]);

    v1 *= v2;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v1[j], fix.lhs[i][j] * fix.rhs[i][j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check v1 * v2", V, VECTOR_REAL_TYPES) {
  VTypeRealTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);
    V v2(&fix.rhs[i][0]);

    V v3 = v1 * v2;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v3[j], fix.lhs[i][j] * fix.rhs[i][j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check v1 /= v2", V, VECTOR_REAL_TYPES) {
  VTypeRealTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);
    V v2(&fix.rhs[i][0]);

    v1 /= v2;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v1[j], fix.lhs[i][j] / fix.rhs[i][j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check v1 / v2", V, VECTOR_REAL_TYPES) {
  VTypeRealTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);
    V v2(&fix.rhs[i][0]);

    V v3 = v1 / v2;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v3[j], fix.lhs[i][j] / fix.rhs[i][j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check -v", V, VECTOR_REAL_TYPES) {
  VTypeRealTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);

    V v2 = -v1;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v2[j], -fix.lhs[i][j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check v1 < v2", V, VECTOR_REAL_TYPES) {
  VTypeRealTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);
    V v2(&fix.rhs[i][0]);

    typename V::mask_t m = v1 < v2;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(m[j], fix.lhs[i][j] < fix.rhs[i][j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check v1 <= v2", V, VECTOR_REAL_TYPES) {
  VTypeRealTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);
    V v2(&fix.rhs[i][0]);

    typename V::mask_t m = v1 <= v2;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(m[j], fix.lhs[i][j] <= fix.rhs[i][j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check v1 > v2", V, VECTOR_REAL_TYPES) {
  VTypeRealTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);
    V v2(&fix.rhs[i][0]);

    typename V::mask_t m = v1 > v2;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(m[j], fix.lhs[i][j] > fix.rhs[i][j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check v1 >= v2", V, VECTOR_REAL_TYPES) {
  VTypeRealTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);
    V v2(&fix.rhs[i][0]);

    typename V::mask_t m = v1 >= v2;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(m[j], fix.lhs[i][j] >= fix.rhs[i][j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check v1 == v2", V, VECTOR_REAL_TYPES) {
  VTypeRealTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);
    V v2(&fix.rhs[i][0]);

    typename V::mask_t m = v1 == v2;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(m[j], fix.lhs[i][j] == fix.rhs[i][j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check v1 != v2", V, VECTOR_REAL_TYPES) {
  VTypeRealTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);
    V v2(&fix.rhs[i][0]);

    typename V::mask_t m = v1 != v2;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(m[j], fix.lhs[i][j] != fix.rhs[i][j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check select(c, t, f)", V, VECTOR_REAL_TYPES) {
  VTypeRealTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);
    V v2(&fix.rhs[i][0]);

    V v3 = select(v1 < v2, v1, v2);

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v3[j], fix.lhs[i][j] < fix.rhs[i][j] ? fix.lhs[i][j] : fix.rhs[i][j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check store(mem, v)", V, VECTOR_REAL_TYPES) {
  VTypeRealTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v(&fix.lhs[i][0]);
    typename V::scalar_t mem[V::W];

    mem_store(mem, v);

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(mem[j], fix.lhs[i][j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check fabs(v)", V, VECTOR_REAL_TYPES) {
  VTypeRealTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);

    V v2 = fabs(v1);

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v2[j], fabs(fix.lhs[i][j]));
    }
  }
}

TEST_CASE_TEMPLATE("Check abs(v)", V, VECTOR_REAL_TYPES) {
  VTypeRealTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);

    V v2 = abs(v1);

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v2[j], abs(fix.lhs[i][j]));
    }
  }
}

TEST_CASE_TEMPLATE("Check max(v1, v2)", V, VECTOR_REAL_TYPES) {
  VTypeRealTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);
    V v2(&fix.rhs[i][0]);

    V v3 = max(v1, v2);

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v3[j], max(fix.lhs[i][j], fix.rhs[i][j]));
    }
  }
}

TEST_CASE_TEMPLATE("Check min(v1, v2)", V, VECTOR_REAL_TYPES) {
  VTypeRealTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);
    V v2(&fix.rhs[i][0]);

    V v3 = min(v1, v2);

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v3[j], min(fix.lhs[i][j], fix.rhs[i][j]));
    }
  }
}

TEST_CASE_TEMPLATE("Check ceil(v)", V, VECTOR_REAL_TYPES) {
  VTypeRealTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);

    V v2 = ceil(v1);

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v2[j], ceil(fix.lhs[i][j]));
    }
  }
}

TEST_CASE_TEMPLATE("Check floor(v)", V, VECTOR_REAL_TYPES) {
  VTypeRealTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);

    V v2 = floor(v1);

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v2[j], floor(fix.lhs[i][j]));
    }
  }
}

TEST_CASE_TEMPLATE("Check round(v)", V, VECTOR_REAL_TYPES) {
  VTypeRealTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);

    V v2 = round(v1);

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v2[j], round(fix.lhs[i][j]));
    }
  }
}

TEST_CASE_TEMPLATE("Check sqrt(v)", V, VECTOR_REAL_TYPES) {
  VTypeRealTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);

    V v2 = sqrt(fabs(v1));

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v2[j], sqrt(fabs(fix.lhs[i][j])));
    }
  }
}

// TODO: tolerance is too large, but checks for completeness of implementation
// of vtype. Possible remedy: use a couple of Newton iterations to improve
// accuracy of the approximate reciprocal square root calculated using the SIMD
// intrinsics.
TEST_CASE_TEMPLATE("Check rsqrt(v)", V, VECTOR_REAL_TYPES) {
  VTypeRealTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);

    V v2 = rsqrt(fabs(v1));

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v2[j], doctest::Approx(((typename V::scalar_t)1.0)/sqrt(fabs(fix.lhs[i][j]))).epsilon(1e-3));
    }
  }
}

template<typename V>
struct VTypeIntTestFixture {
  size_t size;
  vector<std::array<typename V::scalar_t,V::W>> lhs;
  vector<std::array<typename V::scalar_t,V::W>> rhs;
  vector<typename V::scalar_t> scalar;

  VTypeIntTestFixture();
};

template<typename V>
VTypeIntTestFixture<V>::VTypeIntTestFixture() {
  mt19937 rng(23764);
  uniform_int_distribution<typename V::scalar_t> dist(-100, 100);
  size = 100;
  lhs.resize(size);
  rhs.resize(size);
  scalar.resize(size);
  for(size_t i = 0; i < size; ++i) {
    for(unsigned int j = 0; j < V::W; ++j) {
      lhs[i][j] = dist(rng);
      rhs[i][j] = dist(rng);
      scalar[i] = dist(rng);
    }
  }
}

#define VECTOR_INT_TYPES vtype<int32_t, 8>

TEST_CASE_TEMPLATE("Check constructor vtype(scalar_t)", V, VECTOR_INT_TYPES) {
  VTypeIntTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v(fix.scalar[i]);

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v[j], fix.scalar[i]);
    }
  }
}

TEST_CASE_TEMPLATE("Check constructor vtype(scalar_t...)", V, VECTOR_INT_TYPES) {
  VTypeIntTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    if constexpr(is_same_v<V, vtype<int32_t, 8>>) {
      V v(fix.lhs[i][0], fix.lhs[i][1], fix.lhs[i][2], fix.lhs[i][3],
          fix.lhs[i][4], fix.lhs[i][5], fix.lhs[i][6], fix.lhs[i][7]);

      for(unsigned int j = 0; j < V::W; ++j) {
        CHECK_EQ(v[j], fix.lhs[i][j]);
      }
    }
  }
}

TEST_CASE_TEMPLATE("Check copy constructor and assignment operator", V, VECTOR_INT_TYPES) {
  VTypeIntTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);

    V v2(v1);

    V v3;
    v3 = v2;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v2[j], fix.lhs[i][j]);
    }

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v3[j], fix.lhs[i][j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check constructor vtype(scalar_t const * mem)", V, VECTOR_INT_TYPES) {
  VTypeIntTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v(&fix.lhs[i][0]);

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v[j], fix.lhs[i][j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check v1 += v2", V, VECTOR_INT_TYPES) {
  VTypeIntTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);
    V v2(&fix.rhs[i][0]);

    v1 += v2;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v1[j], fix.lhs[i][j] + fix.rhs[i][j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check v1 + v2", V, VECTOR_INT_TYPES) {
  VTypeIntTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);
    V v2(&fix.rhs[i][0]);

    V v3 = v1 + v2;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v3[j], fix.lhs[i][j] + fix.rhs[i][j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check v1 -= v2", V, VECTOR_INT_TYPES) {
  VTypeIntTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);
    V v2(&fix.rhs[i][0]);

    v1 -= v2;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v1[j], fix.lhs[i][j] - fix.rhs[i][j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check v1 - v2", V, VECTOR_INT_TYPES) {
  VTypeIntTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);
    V v2(&fix.rhs[i][0]);

    V v3 = v1 - v2;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v3[j], fix.lhs[i][j] - fix.rhs[i][j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check v1 *= v2", V, VECTOR_INT_TYPES) {
  VTypeIntTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);
    V v2(&fix.rhs[i][0]);

    v1 *= v2;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v1[j], fix.lhs[i][j] * fix.rhs[i][j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check v1 * v2", V, VECTOR_INT_TYPES) {
  VTypeIntTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    V v1(&fix.lhs[i][0]);
    V v2(&fix.rhs[i][0]);

    V v3 = v1 * v2;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v3[j], fix.lhs[i][j] * fix.rhs[i][j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check v1 /= v2", V, VECTOR_INT_TYPES) {
  VTypeIntTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    std::array<typename V::scalar_t, V::W> denom;
    for(unsigned int j = 0; j < V::W; ++j) {
      if(fix.rhs[i][j] == 0) {
        denom[j] = 1;
      } else {
        denom[j] = fix.rhs[i][j];
      }
    }

    V v1(&fix.lhs[i][0]);
    V v2(&denom[0]);

    v1 /= v2;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v1[j], fix.lhs[i][j] / denom[j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check v1 / v2", V, VECTOR_INT_TYPES) {
  VTypeIntTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    std::array<typename V::scalar_t, V::W> denom;
    for(unsigned int j = 0; j < V::W; ++j) {
      if(fix.rhs[i][j] == 0) {
        denom[j] = 1;
      } else {
        denom[j] = fix.rhs[i][j];
      }
    }

    V v1(&fix.lhs[i][0]);
    V v2(&denom[0]);

    V v3 = v1 / v2;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v3[j], fix.lhs[i][j] / denom[j]);
    }
  }
}

// TODO: SIMD vectorized code for remainder has sign error
//TEST_CASE_TEMPLATE("Check v1 %= v2", V, VECTOR_INT_TYPES) {
//  VTypeIntTestFixture<V> fix;
//  for(size_t i = 0; i < fix.size; ++i) {
//    std::array<typename V::scalar_t, V::W> num;
//    for(unsigned int j = 0; j < V::W; ++j) {
//      num[j] = fix.lhs[i][j];
//    }
//
//    std::array<typename V::scalar_t, V::W> denom;
//    for(unsigned int j = 0; j < V::W; ++j) {
//      if(fix.rhs[i][j] == 0) {
//        denom[j] = 1;
//      } else {
//        denom[j] = fix.rhs[i][j];
//      }
//    }
//
//    V v1(&num[0]);
//    V v2(&denom[0]);
//
//    v1 %= v2;
//
//    for(unsigned int j = 0; j < V::W; ++j) {
//      CHECK_EQ(v1[j], num[j] % denom[j]);
//    }
//  }
//}

// TODO: SIMD vectorized code for remainder has sign error
//TEST_CASE_TEMPLATE("Check v1 %= v2", V, VECTOR_INT_TYPES) {
//  VTypeIntTestFixture<V> fix;
//  for(size_t i = 0; i < fix.size; ++i) {
//    std::array<typename V::scalar_t, V::W> num;
//    for(unsigned int j = 0; j < V::W; ++j) {
//      num[j] = fix.lhs[i][j];
//    }
//
//    std::array<typename V::scalar_t, V::W> denom;
//    for(unsigned int j = 0; j < V::W; ++j) {
//      if(fix.rhs[i][j] == 0) {
//        denom[j] = 1;
//      } else {
//        denom[j] = fix.rhs[i][j];
//      }
//    }
//
//    V v1(&num[0]);
//    V v2(&denom[0]);
//
//    V v3 = v1 % v2;
//
//    for(unsigned int j = 0; j < V::W; ++j) {
//      CHECK_EQ(v3[j], num[j] % denom[j]);
//    }
//  }
//}

TEST_CASE_TEMPLATE("Check v1 <<= v2", V, VECTOR_INT_TYPES) {
  VTypeIntTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    std::array<typename V::scalar_t, V::W> shifts;
    for(unsigned int j = 0; j < V::W; ++j) {
      shifts[j] = fix.rhs[i][j] < 0 ? -fix.rhs[i][j] : fix.rhs[i][j];
      shifts[j] = shifts[j]%(sizeof(typename V::scalar_t)*8);
    }

    V v1(&fix.lhs[i][0]);
    V v2(&shifts[0]);

    v1 <<= v2;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v1[j], fix.lhs[i][j] << shifts[j]);
    }
  }
}

// TODO: The issue with this test is that the right shift operation on signed
// integers is implementation defined. In arithmetic right shift, sign bit is
// filled, in logical right shift, zero is filled. For now, vtype implements
// arithmetic shift. So this test may fail on implementations which use logical
// shift.
TEST_CASE_TEMPLATE("Check v1 >>= v2", V, VECTOR_INT_TYPES) {
  VTypeIntTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    std::array<typename V::scalar_t, V::W> shifts;
    for(unsigned int j = 0; j < V::W; ++j) {
      shifts[j] = fix.rhs[i][j] < 0 ? -fix.rhs[i][j] : fix.rhs[i][j];
      shifts[j] = shifts[j]%(sizeof(typename V::scalar_t)*8);
    }

    V v1(&fix.lhs[i][0]);
    V v2(&shifts[0]);

    v1 >>= v2;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v1[j], fix.lhs[i][j] >> shifts[j]);
    }
  }
}

TEST_CASE_TEMPLATE("Check v1 <<= s", V, VECTOR_INT_TYPES) {
  VTypeIntTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    typename V::scalar_t shift = fix.scalar[i] < 0 ? -fix.scalar[i] : fix.scalar[i];
    shift = shift%(sizeof(typename V::scalar_t)*8);

    V v1(&fix.lhs[i][0]);

    v1 <<= shift;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v1[j], fix.lhs[i][j] << shift);
    }
  }
}

// TODO: The issue with this test is that the right shift operation on signed
// integers is implementation defined. In arithmetic right shift, sign bit is
// filled, in logical right shift, zero is filled. For now, vtype implements
// arithmetic shift. So this test may fail on implementations which use logical
// shift.
TEST_CASE_TEMPLATE("Check v1 >>= s", V, VECTOR_INT_TYPES) {
  VTypeIntTestFixture<V> fix;
  for(size_t i = 0; i < fix.size; ++i) {
    typename V::scalar_t shift = fix.scalar[i] < 0 ? -fix.scalar[i] : fix.scalar[i];
    shift = shift%(sizeof(typename V::scalar_t)*8);

    V v1(&fix.lhs[i][0]);

    v1 >>= shift;

    for(unsigned int j = 0; j < V::W; ++j) {
      CHECK_EQ(v1[j], fix.lhs[i][j] >> shift);
    }
  }
}
