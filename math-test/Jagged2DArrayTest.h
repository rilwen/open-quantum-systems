#pragma once

#include <gtest/gtest.h>
#include "math/Jagged2DArray.h"
#include <list>
#include <vector> 

using namespace rql;

template <class T>
class Jagged2DArrayTest: public testing::Test
{
};

TYPED_TEST_CASE_P(Jagged2DArrayTest);

TYPED_TEST_P(Jagged2DArrayTest, ConstructAndAssign)
{
    math::Jagged2DArray<TypeParam> a(2ul, 3ul);
    EXPECT_EQ(2u, a.size());
    EXPECT_EQ(6u, a.nbr_elements());
    EXPECT_EQ(3u, a.row_size(0));
    EXPECT_EQ(3u, a.row_size(1));
    TypeParam v(1);
    math::Jagged2DArray<TypeParam> b(static_cast<size_t>(2), static_cast<size_t>(3), static_cast<TypeParam>(v));
    EXPECT_EQ(2u, b.size());
    EXPECT_EQ(6u, b.nbr_elements());
    EXPECT_EQ(3u, b.row_size(0));
    EXPECT_EQ(3u, b.row_size(1));
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            EXPECT_EQ(v, b(i,j));
        }
    }
    b(0, 0) = 0;
    EXPECT_EQ(static_cast<TypeParam>(0), b(0,0));
    math::Jagged2DArray<TypeParam> c(b);
    EXPECT_EQ(b.size(), c.size());
    EXPECT_EQ(b.nbr_elements(), c.nbr_elements());
    for (size_t i = 0; i < c.size(); ++i) {
        EXPECT_EQ(b.row_size(i), c.row_size(i));
        for (size_t j = 0; j < c.row_size(i); ++j) {
            EXPECT_EQ(b(i,j), c(i,j));
        }
    }
    math::Jagged2DArray<TypeParam> d;
    EXPECT_EQ(0u, d.size());
    EXPECT_EQ(0u, d.nbr_elements());
    d = a;
    EXPECT_EQ(a.size(), d.size());
    EXPECT_EQ(a.nbr_elements(), d.nbr_elements());
    for (size_t i = 0; i < d.size(); ++i) {
        EXPECT_EQ(a.row_size(i), d.row_size(i));
        for (size_t j = 0; j < d.row_size(i); ++j) {
            EXPECT_EQ(a(i,j), d(i,j));
        }
    }

    std::list<std::vector<TypeParam> > list(3);
    typename std::list<std::vector<TypeParam> >::iterator list_iter = list.begin();
    ++list_iter;
    list_iter->push_back(1);
    ++list_iter;
    list_iter->push_back(2);
    list_iter->push_back(10);
    math::Jagged2DArray<TypeParam> e(list);
    EXPECT_EQ(list.size(), e.size());
    size_t idx = 0;
    size_t total_nbr_elems = 0;
    for (typename std::list<std::vector<TypeParam> >::const_iterator i = list.begin(); i != list.end(); ++i) {
        EXPECT_EQ(i->size(), e.row_size(idx));
        total_nbr_elems += i->size();
        for (size_t j = 0; j < i->size(); ++j) {
            EXPECT_EQ(i->at(j), e(idx, j));
        }
        ++idx;
    }
    EXPECT_EQ(total_nbr_elems, e.nbr_elements());

    std::list<int> sizes(3);
    sizes.push_back(10);
    sizes.push_back(0);
    sizes.push_back(2);
    math::Jagged2DArray<TypeParam> f(sizes.begin(), sizes.end());
    EXPECT_EQ(sizes.size(), f.size());
    idx = 0;
    for (std::list<int>::const_iterator i = sizes.begin(); i != sizes.end(); ++i) {
        EXPECT_EQ(static_cast<size_t>(*i), f.row_size(idx));
        ++idx;
    }
}

TYPED_TEST_P(Jagged2DArrayTest, RowAccess)
{
    TypeParam v(1);
    math::Jagged2DArray<TypeParam> a(static_cast<size_t>(2), static_cast<size_t>(3), static_cast<TypeParam>(2));
    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < a.row_size(i); ++j) {
            EXPECT_EQ(a(i, j), a[i][j]);
            a[i][j] = v;
            math::Jagged2DArrayRowRef<TypeParam> row(a[i]);
            std::fill(row.begin(), row.end(), i);
            EXPECT_EQ(a.row_size(i), row.size());
        }
    }

    const math::Jagged2DArray<TypeParam> b(a);
    for (size_t i = 0; i < b.size(); ++i) {
        for (size_t j = 0; j < b.row_size(i); ++j) {
            EXPECT_EQ(static_cast<TypeParam>(i), b[i][j]);
        }
    }

    a(0u,0u) = static_cast<TypeParam>(1);
    a(1u,0u) = static_cast<TypeParam>(2);
    math::Jagged2DArrayRowRef<TypeParam> r(a[0]);
    EXPECT_EQ(a(0,0), r[0]);
    r.reset(a[1]);
    EXPECT_EQ(a(1,0), r[0]);
}

TYPED_TEST_P(Jagged2DArrayTest, Iterators)
{
  math::Jagged2DArray<TypeParam> a(static_cast<size_t>(2), 3ul);
    TypeParam v(1);
    std::fill(a.flat_begin(), a.flat_end(), v);
    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < a.row_size(i); ++j) {
            EXPECT_EQ(v, a[i][j]);
        }
    }
}

TYPED_TEST_P(Jagged2DArrayTest, Swap)
{
    const size_t ra = 2;
    const size_t ca = 3;
    const size_t rb = 2;
    const size_t cb = 1;
    math::Jagged2DArray<TypeParam> a(ra, ca);
    math::Jagged2DArray<TypeParam> b(rb, cb);
    TypeParam va(1);
    TypeParam vb(2);
    std::fill(a.flat_begin(), a.flat_end(), va);
    std::fill(b.flat_begin(), b.flat_end(), vb);
    a.swap(b);
    EXPECT_EQ(rb, a.size());
    EXPECT_EQ(ra, b.size());
    for (size_t i = 0; i < a.size(); ++i) {
        EXPECT_EQ(cb, a.row_size(i));
        for (size_t j = 0; j < a.row_size(i); ++j) {
            EXPECT_EQ(vb, a(i, j));
        }
    }
    for (size_t i = 0; i < b.size(); ++i) {
        EXPECT_EQ(ca, b.row_size(i));
        for (size_t j = 0; j < b.row_size(i); ++j) {
            EXPECT_EQ(va, b(i, j));
        }
    }

    // go back
    a.swap(b);

    math::swap(a, b);
    EXPECT_EQ(rb, a.size());
    EXPECT_EQ(ra, b.size());
    for (size_t i = 0; i < a.size(); ++i) {
        EXPECT_EQ(cb, a.row_size(i));
        for (size_t j = 0; j < a.row_size(i); ++j) {
            EXPECT_EQ(vb, a(i, j));
        }
    }
    for (size_t i = 0; i < b.size(); ++i) {
        EXPECT_EQ(ca, b.row_size(i));
        for (size_t j = 0; j < b.row_size(i); ++j) {
            EXPECT_EQ(va, b(i, j));
        }
    }
}

TYPED_TEST_P(Jagged2DArrayTest, Triangular)
{
	typedef typename math::Jagged2DArray<TypeParam> my_type;
	TypeParam v(1);
	const size_t N = 3;
	my_type a1(N, false);
	my_type a2(N, true);
	my_type a3(N, false, v);
	my_type a4(N, true, v);
	ASSERT_EQ(N, a1.size());
	ASSERT_EQ(N, a2.size());
	ASSERT_EQ(N, a3.size());
	ASSERT_EQ(N, a4.size());
	const size_t nbr_elems_true = 6;
	const size_t nbr_elems_false = 3;
	ASSERT_EQ(nbr_elems_false, a1.nbr_elements());
	ASSERT_EQ(nbr_elems_true, a2.nbr_elements());
	ASSERT_EQ(nbr_elems_false, a3.nbr_elements());
	ASSERT_EQ(nbr_elems_true, a4.nbr_elements());
	for (size_t i = 0; i < N; ++i) {
		ASSERT_EQ(i, a1.row_size(i));
		ASSERT_EQ(i + 1, a2.row_size(i));
		ASSERT_EQ(i, a3.row_size(i));
		ASSERT_EQ(i + 1, a4.row_size(i));
	}
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < i; ++j) {
			ASSERT_EQ(v, a3(i, j));
		}
		for (size_t j = 0; j <= i; ++j) {
			ASSERT_EQ(v, a4(i, j));
		}
	}
}

//TYPED_TEST_P(Jagged2DArrayTest, Arithmetic)
//{
//	typedef typename math::Jagged2DArray<TypeParam> my_type;
//	my_type a(2u, 2u);
//	a(0,0)=10;
//	a(0,1)=1;
//	a(1,0)=2;
//	a(1,1)=3;
//	my_type b(2u, 2u);
//	b(0,0)=b(1,0)=b(0,1)=b(1,1)=1;
//	a -= b;
//	EXPECT_EQ(static_cast<TypeParam>(9), a(0,0));
//	EXPECT_EQ(static_cast<TypeParam>(0), a(0,1));
//	EXPECT_EQ(static_cast<TypeParam>(1), a(1,0));
//	EXPECT_EQ(static_cast<TypeParam>(2), a(1,1));
//	a += b;
//	EXPECT_EQ(static_cast<TypeParam>(10), a(0,0));
//	EXPECT_EQ(static_cast<TypeParam>(1), a(0,1));
//	EXPECT_EQ(static_cast<TypeParam>(2), a(1,0));
//	EXPECT_EQ(static_cast<TypeParam>(3), a(1,1));
//	a *= static_cast<TypeParam>(2);
//	EXPECT_EQ(static_cast<TypeParam>(20), a(0,0));
//	EXPECT_EQ(static_cast<TypeParam>(2), a(0,1));
//	EXPECT_EQ(static_cast<TypeParam>(4), a(1,0));
//	EXPECT_EQ(static_cast<TypeParam>(6), a(1,1));
//	a /= static_cast<TypeParam>(2);
//	EXPECT_EQ(static_cast<TypeParam>(10), a(0,0));
//	EXPECT_EQ(static_cast<TypeParam>(1), a(0,1));
//	EXPECT_EQ(static_cast<TypeParam>(2), a(1,0));
//	EXPECT_EQ(static_cast<TypeParam>(3), a(1,1));
//}

REGISTER_TYPED_TEST_CASE_P(Jagged2DArrayTest, ConstructAndAssign, RowAccess, Iterators, Swap, Triangular);
