//Link to Boost
#define BOOST_TEST_DYN_LINK

//Define our Module name (prints at testing)
#define BOOST_TEST_MODULE BetheSolver UnitTests

//VERY IMPORTANT - include this last
#include <boost/test/included/unit_test.hpp>
//#include <boost/test/unit_test.hpp>

#include "test.h"
#include "../PowerSeries.h"

// test suite
BOOST_FIXTURE_TEST_SUITE(ExpBetaSeries_suite, SimpleTestFixture, * utf::label("ExpBetaSeries"))

BOOST_DATA_TEST_CASE(constructor_test, bdata::random(-10, 10) ^ bdata::xrange(20), val, index)
{
    ExpBetaSeries ebs(val);
    
    BOOST_TEST(ebs == val);
    if (val != 0) {
        BOOST_TEST(ebs.TermNumber() == 1);
    } else {
        BOOST_TEST(ebs.IsZero());
        BOOST_TEST(ebs.TermNumber() == 0);
    }
    
    BOOST_TEST(ebs.FirstTermPositive() == (val >0));
    BOOST_TEST(ebs.ToNumerical(0.2) == complex<f_type>((f_type)val, 0.0L));
}

BOOST_DATA_TEST_CASE(constructor_test2, bdata::random(-5, 5) ^ bdata::random(-10, 10) ^ bdata::xrange(20), order, val, index)
{
    ExpBetaSeries ebs(order, val);
    
    BOOST_TEST((ebs == val) == (val == 0 || order == 0));
    if (val != 0) {
        BOOST_TEST(ebs.TermNumber() == 1);
    } else {
        BOOST_TEST(ebs.IsZero());
        BOOST_TEST(ebs.TermNumber() == 0);
    }
    
    BOOST_TEST(ebs.FirstTermPositive() == (val >0));
    f_type beta = 0.4;
    BOOST_TEST(ebs.ToNumerical(beta) == ExpI(beta * order) * (f_type)val);
}

BOOST_DATA_TEST_CASE(shiftOrder_test, bdata::random(-5, 5) ^ bdata::random(-10, 10) ^ bdata::random(-5, 5) ^ bdata::xrange(20), order, val, n, index)
{
    f_type beta = 0.6L;
    ExpBetaSeries ebs(order, val);
    
    complex<f_type> c = ebs.ToNumerical(beta);
    ebs.ShiftOrder(n);
    BOOST_TEST(std::abs(ebs.ToNumerical(beta) - c * ExpI(beta * n)) < eps);
}

BOOST_DATA_TEST_CASE(plus_test, bdata::random(-5, 5) ^ bdata::random(-10, 10) ^ bdata::random(-5, 5) ^ bdata::random(-5, 5) ^ bdata::xrange(20), order1, val1, order2, val2, index)
{
    f_type beta = 0.9L;
    ExpBetaSeries ebs1(order1, val1);
    ExpBetaSeries ebs2(order2, val2);
    ExpBetaSeries ebs3 = ebs1 + ebs2;
    
    BOOST_TEST(std::abs(ebs3.ToNumerical(beta) - (ebs1.ToNumerical(beta) + ebs2.ToNumerical(beta))) < eps);
}

BOOST_DATA_TEST_CASE(plus_test2, bdata::random(-5, 5) ^ bdata::random(-10, 10) ^ bdata::random(-5, 5) ^ bdata::random(-10, 10) ^ bdata::xrange(20), order1, val1, order2, val2, index)
{
    f_type beta = 0.9L;
    ExpBetaSeries ebs1(order1, val1);
    ExpBetaSeries ebs2(order2, val2);
    std::complex<f_type> c1 = ebs1.ToNumerical(beta);
    std::complex<f_type> c2 = ebs2.ToNumerical(beta);
    ebs2 += ebs1;
    
    BOOST_TEST(std::abs(ebs2.ToNumerical(beta) - (c1 + c2)) < eps);
}

BOOST_DATA_TEST_CASE(multiply_integer_test, bdata::random(-5, 5) ^ bdata::random(-10, 10) ^ bdata::random(-10, 10) ^ bdata::xrange(20), order1, val1, n, index)
{
    f_type beta = 0.9L;
    ExpBetaSeries ebs1(order1, val1);
    std::complex<f_type> c1 = ebs1.ToNumerical(beta);
    ExpBetaSeries ebs2 = ebs1 * n;
    
    BOOST_TEST(std::abs(ebs2.ToNumerical(beta) - (c1 * (f_type)n)) < eps);
}

BOOST_DATA_TEST_CASE(multiply_integer_test2, bdata::random(-5, 5) ^ bdata::random(-10, 10) ^ bdata::random(-10, 10) ^ bdata::xrange(20), order1, val1, n, index)
{
    f_type beta = 0.9L;
    ExpBetaSeries ebs1(order1, val1);
    std::complex<f_type> c1 = ebs1.ToNumerical(beta);
    ebs1 *= n;
    
    BOOST_TEST(std::abs(ebs1.ToNumerical(beta) - (c1 * (f_type)n)) < eps);
}

BOOST_DATA_TEST_CASE(multiply_test, bdata::random(-5, 5) ^ bdata::random(-10, 10) ^ bdata::random(-5, 5) ^ bdata::random(-10, 10) ^ bdata::xrange(20), order1, val1, order2, val2, index)
{
    f_type beta = 0.9L;
    ExpBetaSeries ebs1(order1, val1);
    ExpBetaSeries ebs2(order2, val2);
    ExpBetaSeries ebs3 = ebs1 * ebs2;
    
    BOOST_TEST(std::abs(ebs3.ToNumerical(beta) - (ebs2.ToNumerical(beta) * ebs1.ToNumerical(beta))) < eps);
}

BOOST_DATA_TEST_CASE(multiply_test2, bdata::random(-5, 5) ^ bdata::random(-10, 10) ^ bdata::random(-5, 5) ^ bdata::random(-10, 10) ^ bdata::xrange(20), order1, val1, order2, val2, index)
{
    f_type beta = 0.9L;
    ExpBetaSeries ebs1(order1, val1);
    ExpBetaSeries ebs2(order2, val2);
    complex<f_type> c1 = ebs1.ToNumerical(beta);
    complex<f_type> c2 = ebs2.ToNumerical(beta);
    ebs2 *= ebs1;
    
    BOOST_TEST(std::abs(ebs2.ToNumerical(beta) - (c1 * c2)) < eps);
}

BOOST_AUTO_TEST_SUITE_END()

// test suite
BOOST_FIXTURE_TEST_SUITE(NExpansionSeries_suite, SimpleTestFixture, * utf::label("NExpansionSeries"))

BOOST_DATA_TEST_CASE(constructor_test, bdata::random(-5, 5) ^ bdata::random(-5, 5) ^ bdata::random(-10, 10) ^ bdata::xrange(20), nOrder, eOrder, val, index)
{
    int N = 11;
    f_type beta = 1.2;
    NExpansionSeries ne(nOrder, ExpBetaSeries(eOrder, val));
    complex<f_type> expected = ExpI(beta * eOrder) * (f_type)val * pow((f_type)N, nOrder);
    
    BOOST_TEST(std::abs(ne.ToNumerical(N, beta) - expected) < eps);
}

BOOST_DATA_TEST_CASE(plus_test, bdata::random(-5, 5) ^ bdata::random(-5, 5) ^ bdata::random(-10, 10) ^ bdata::random(-5, 5) ^ bdata::random(-5, 5) ^ bdata::random(-10, 10) ^ bdata::xrange(20), nOrder1, eOrder1, val1, nOrder2, eOrder2, val2, index)
{
    int N = 11;
    f_type beta = 1.2;
    NExpansionSeries ne1(nOrder1, ExpBetaSeries(eOrder1, val1));
    NExpansionSeries ne2(nOrder2, ExpBetaSeries(eOrder2, val2));
    NExpansionSeries ne3 = ne1 + ne2;
    complex<f_type> c1 = ne1.ToNumerical(N, beta);
    complex<f_type> c2 = ne2.ToNumerical(N, beta);
    
    complex<f_type> expected = c1 + c2;
    
    BOOST_TEST(std::abs(ne3.ToNumerical(N, beta) - expected) < eps);
}

BOOST_DATA_TEST_CASE(plus_test2, bdata::random(-5, 5) ^ bdata::random(-5, 5) ^ bdata::random(-10, 10) ^ bdata::random(-5, 5) ^ bdata::random(-5, 5) ^ bdata::random(-10, 10) ^ bdata::xrange(20), nOrder1, eOrder1, val1, nOrder2, eOrder2, val2, index)
{
    int N = 11;
    f_type beta = 1.2;
    NExpansionSeries ne1(nOrder1, ExpBetaSeries(eOrder1, val1));
    NExpansionSeries ne2(nOrder2, ExpBetaSeries(eOrder2, val2));
    complex<f_type> c1 = ne1.ToNumerical(N, beta);
    complex<f_type> c2 = ne2.ToNumerical(N, beta);
    ne2 += ne1;
    
    complex<f_type> expected = c1 + c2;
    
    BOOST_TEST(std::abs(ne2.ToNumerical(N, beta) - expected) < eps);
}

BOOST_DATA_TEST_CASE(multiply_ExpBetaSeries_test, bdata::random(-5, 5) ^ bdata::random(-5, 5) ^ bdata::random(-10, 10) ^  bdata::random(-5, 5) ^ bdata::random(-10, 10) ^ bdata::xrange(20), nOrder1, eOrder1, val1, eOrder2, val2, index)
{
    int N = 11;
    f_type beta = 1.2;
    NExpansionSeries ne1(nOrder1, ExpBetaSeries(eOrder1, val1));
    ExpBetaSeries ebs(eOrder2, val2);
    NExpansionSeries ne3 = ne1 * ebs;
    complex<f_type> c1 = ne1.ToNumerical(N, beta);
    complex<f_type> c2 = ebs.ToNumerical(beta);
    
    complex<f_type> expected = c1 * c2;
    
    BOOST_TEST(std::abs(ne3.ToNumerical(N, beta) - expected) < eps);
}

BOOST_DATA_TEST_CASE(multiply_ExpBetaSeries_test2, bdata::random(-5, 5) ^ bdata::random(-5, 5) ^ bdata::random(-10, 10) ^  bdata::random(-5, 5) ^ bdata::random(-10, 10) ^ bdata::xrange(20), nOrder1, eOrder1, val1, eOrder2, val2, index)
{
    int N = 11;
    f_type beta = 1.2;
    NExpansionSeries ne1(nOrder1, ExpBetaSeries(eOrder1, val1));
    ExpBetaSeries ebs(eOrder2, val2);
    complex<f_type> c1 = ne1.ToNumerical(N, beta);
    complex<f_type> c2 = ebs.ToNumerical(beta);

    ne1 *= ebs;
    
    complex<f_type> expected = c1 * c2;
    
    BOOST_TEST(std::abs(ne1.ToNumerical(N, beta) - expected) < eps);
}

BOOST_DATA_TEST_CASE(multiply_test, bdata::random(-5, 5) ^ bdata::random(-5, 5) ^ bdata::random(-10, 10) ^ bdata::random(-5, 5) ^ bdata::random(-5, 5) ^ bdata::random(-10, 10) ^ bdata::xrange(20), nOrder1, eOrder1, val1, nOrder2, eOrder2, val2, index)
{
    int N = 7;
    f_type beta = 1.2;
    NExpansionSeries ne1(nOrder1, ExpBetaSeries(eOrder1, val1));
    NExpansionSeries ne2(nOrder2, ExpBetaSeries(eOrder2, val2));
    NExpansionSeries ne3 = ne1 * ne2;
    complex<f_type> c1 = ne1.ToNumerical(N, beta);
    complex<f_type> c2 = ne2.ToNumerical(N, beta);
    
    complex<f_type> expected = c1 * c2;
    
    BOOST_TEST(std::abs(ne3.ToNumerical(N, beta) - expected) < eps);
}

BOOST_DATA_TEST_CASE(multiply_test2, bdata::random(-5, 5) ^ bdata::random(-5, 5) ^ bdata::random(-10, 10) ^ bdata::random(-5, 5) ^ bdata::random(-5, 5) ^ bdata::random(-10, 10) ^ bdata::xrange(20), nOrder1, eOrder1, val1, nOrder2, eOrder2, val2, index)
{
    int N = 7;
    f_type beta = 1.2;
    NExpansionSeries ne1(nOrder1, ExpBetaSeries(eOrder1, val1));
    NExpansionSeries ne2(nOrder2, ExpBetaSeries(eOrder2, val2));
    complex<f_type> c1 = ne1.ToNumerical(N, beta);
    complex<f_type> c2 = ne2.ToNumerical(N, beta);
    ne2 *= ne1;
    
    complex<f_type> expected = c1 * c2;
    
    BOOST_TEST(std::abs(ne2.ToNumerical(N, beta) - expected) < eps);
}


BOOST_AUTO_TEST_SUITE_END()

