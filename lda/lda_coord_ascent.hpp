/*
 *  Copyright (c) 2014 Alp Kucukelbir, alp@cs.columbia.edu
 *
 *  - Adapted from code by Cheng Zhang, Xavi Gratal, Chong Wang,
 *    Matthew D. Hoffman, and David M. Blei
 *  - Licensed under terms of GPLv2 license (see LICENSE)
 *
 */

#ifndef lda_coord_ascent_hpp
#define lda_coord_ascent_hpp

#include <chrono>
#include <iostream>
#include <Eigen/Dense>

#include "lda.hpp"

using namespace Eigen;

namespace lda
{

class lda_coord_ascent : public lda
{
public:

  util::corpus train_corpus_;
  util::corpus test_corpus_;

  int const max_itr_;
  double const step_size_;

  std::string const results_dir_;

  lda_coord_ascent(
      util::corpus& train_corpus,
      util::corpus& test_corpus,
      int const max_itr, double const step_size,
      std::string const& results_dir,
      int num_topics, util::vocabulary const& vocab,
      double alpha, double eta,
      double tau0,  double kappa
      ):
  lda(num_topics, vocab, alpha, eta, tau0, kappa),
  train_corpus_(train_corpus),
  test_corpus_(test_corpus),
  max_itr_(max_itr),
  step_size_(step_size),
  results_dir_(results_dir)
  { }


  int run_lda_coord_ascent(bool compute_elbo, MatrixXd& lambda_init)
  {

    // If every element of lambda_init is zero, initialize from scratch
    if (lambda_init.isZero(0))
      init_lda();
    else
      init_lda(lambda_init);

    double elbo     = 0.0;
    double log_pred = 0.0;

    std::ofstream my_file;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> elapsed_seconds;

    start = std::chrono::system_clock::now();
    for (int i = 0; i < max_itr_; ++i)
    {

      update_lambda(train_corpus_, util::noBOOTSTRAP, step_size_);
      if (compute_elbo == true)
        elbo = calculate_ELBO(train_corpus_, util::VIEW);
      log_pred = calculate_log_pred(test_corpus_,
                                    util::FULL,
                                    util::wobsRATIO);

      std::cout
      << "lda_coord_ascent | itr " << i
      << "\t| elbo = " << elbo
      << "\t| lp = " << log_pred
      << "\t| rho_t = " << rho_t_
      << std::endl;

      my_file.open(results_dir_ + "/runtime_coord_ascent.dat",
                   std::ios_base::app);
      if (my_file.is_open())
      {
        my_file
        << i << "\t,\t"
        << elbo << "\t,\t"
        << log_pred << "\t,\t"
        << rho_t_ << "\t,\t"
        << std::endl;
      }
      my_file.close();

      my_file.open(results_dir_ + "/lambda_coord_ascent.dat");
      if (my_file.is_open())
      {
        my_file << lambda_;
      }
      my_file.close();

    }
    end = std::chrono::system_clock::now();

    std::cout << "----------------------------------------------" << std::endl;
    elapsed_seconds = end-start;
    std::cout << "total seconds: \t" << elapsed_seconds.count()   << std::endl;
    std::cout << "----------------------------------------------" << std::endl;

    return util::CONVERGED;

  }

}; // class lda_coord_ascent

} // namespace lda

#endif


















