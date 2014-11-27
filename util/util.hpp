/*
 *  Copyright (c) 2014 Alp Kucukelbir, alp@cs.columbia.edu
 *
 *  - Adapted from code by Cheng Zhang, Xavi Gratal, Chong Wang,
 *    Matthew D. Hoffman, and David M. Blei
 *  - Licensed under terms of GPLv2 license (see LICENSE)
 *
 */

#ifndef util_hpp
#define util_hpp

#include <iostream>
#include <Eigen/Dense>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>

typedef boost::iostreams::tee_device<std::ostream, std::ofstream> Tee;
typedef boost::iostreams::stream<Tee> TeeStream;

namespace util
{

static const std::string VERSION_STRING = "lda-bump-cpp v1.0";

static const int noBOOTSTRAP = 0;

// wobsRATIO = w_obs/(w_obs+w_ho)
static const double wobsRATIO = 0.5;

static const int FULL   = 10;
static const int VIEW   = 20;

static const int CONVERGED = 1000;

/**
 * define a mersenne twister
 */
std::mt19937_64 mersenne_twister()
{
  std::random_device rd;
  std::mt19937_64 gen(rd());
  return gen;
}

/**
 * returns a bootstrap index of length `length'
 */
std::vector<int> bootstrap_index(int length)
{

  std::vector<int> sample_indices;

  // sample with replacement (bootstrap)
  std::mt19937_64 gen = util::mersenne_twister();
  std::uniform_int_distribution<> unif(0, length - 1 );
  int index;
  for (int i = 0; i < length; ++i)
  {
    index = unif(gen);
    sample_indices.push_back(index);
  }

  return sample_indices;
}

/**
 * returns a range index of length `length'
 */
std::vector<int> range_index(int length)
{

  std::vector<int> range_indices;

  // return range indices
  for (int i = 0; i < length; ++i)
  {
    range_indices.push_back(i);
  }

  return range_indices;
}


double digammad(double x) { return boost::math::digamma<double>(x); }
double lgammad(double x) { return boost::math::lgamma<double>(x); }

/**
 *  Log sum exp trick
 */
double log_sum_exp(Eigen::VectorXd& x)
{
  double result(0.0);
  double xmax(x.maxCoeff());
  x.array() -= xmax;
  result = log((x.array().exp()).sum()) + xmax;
  x.array() += xmax;
  return result;
}

/**
 *  For a vector ~ Dir(alpha), computes E[log(theta)] given alpha.
 */
Eigen::VectorXd dir_exp_log(Eigen::VectorXd const& alpha)
{
  Eigen::VectorXd result = Eigen::VectorXd::Zero(alpha.size());

  // Element-wise evaluation of psi(alpha) [vector <- f(vector)]
  result = alpha.unaryExpr(std::ptr_fun(digammad));

  // Subtract scalar psi(sum(alpha)) from [vector <- vector - scalar]
  result = result.array() - digammad(alpha.sum());

  return result;
}

} // namespace util

#endif
