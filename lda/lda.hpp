/*
 *  Copyright (c) 2014 Alp Kucukelbir, alp@cs.columbia.edu
 *
 *  - Adapted from code by Cheng Zhang, Xavi Gratal, Chong Wang,
 *    Matthew D. Hoffman, and David M. Blei
 *  - Licensed under terms of GPLv2 license (see LICENSE)
 *
 */

#ifndef lda_hpp
#define lda_hpp

#include <random>
#include <iostream>
#include <Eigen/Dense>

#include "corpus.hpp"

using namespace Eigen;

namespace lda
{

class lda
{
public:

  int const num_topics_;          // Number of topics
  util::vocabulary const vocab_;  // Vocabulary
  int const vocab_length_;        // Length of vobaculary

  double const alpha_;    // Hyperparameter prior on topic proportions (\theta)
  double const eta_;      // Hyperparameter prior on topics (\beta)

  double const tau0_;     // Parameter that downweights early iterations
  double const kappa_;    // Learning rate: exponential decay rate. Should be
                          // between (0.5, 1.0] for asymptotic convergence.

  double rho_t_;          // Step-size for stochastic gradient optimization
  int current_iter_;      // Current iteration, used for step size calculation

  MatrixXd lambda_;       // has shape (num_topics_, vocab_length_)
  MatrixXd Elogbeta_;     // has shape (num_topics_, vocab_length_)
  MatrixXd expElogbeta_;  // has shape (num_topics_, vocab_length_)

  std::gamma_distribution<double> gamma_dist_;  // Gamma distribution

  double const e_step_thresh_;   // Threshold for coordinate ascent in E-step
  int    const e_step_iter_lim_; // Iteration limit for E-step

  lda(int num_topics, util::vocabulary const &vocab,
      double alpha, double eta,
      double tau0,  double kappa):
    vocab_(vocab),
    num_topics_(num_topics),
    vocab_length_(vocab_.length_),
    alpha_(alpha),
    eta_(eta),
    tau0_(tau0),
    kappa_(kappa),
    rho_t_(1.0),
    current_iter_(0),
    gamma_dist_(1e2, 1e-2),
    e_step_thresh_(0.01),
    e_step_iter_lim_(5)
  { }


  /**
   * initialize LDA class with random Gamma draws
   */
  void init_lda()
  {
    rho_t_        = 1.0;
    current_iter_ = 0;

    // Initialize the variational distribution q(\beta | \lambda)
    std::mt19937_64 gen = util::mersenne_twister();
    lambda_ = MatrixXd::Zero(num_topics_, vocab_length_);
    for (int j = 0, nRows = lambda_.rows(),
                    nCols = lambda_.cols(); j < nCols; ++j)
    {
      for (int i = 0; i < nRows; ++i)
      {
        lambda_(i, j) = gamma_dist_(gen);
      }
    }

    // Initialize E[log(\beta)]
    Elogbeta_ = MatrixXd::Zero(num_topics_, vocab_length_);
    for (int i = 0, nRows = Elogbeta_.rows(); i < nRows; ++i)
    {
      Elogbeta_.row(i) = util::dir_exp_log(lambda_.row(i));
    }

    // Initialize exp( E[log(\beta)] )
    expElogbeta_ = Elogbeta_.array().exp();
  }


  /**
   * initialize LDA class with MatrixXd parameter
   * @param lambda_init
   */
  void init_lda(MatrixXd& lambda_init)
  {
    rho_t_        = 1.0;
    current_iter_ = 0;

    // Initialize lambda from parameter
    lambda_ = lambda_init;

    // Initialize E[log(\beta)]
    Elogbeta_ = MatrixXd::Zero(num_topics_, vocab_length_);
    for (int i = 0, nRows = Elogbeta_.rows(); i < nRows; ++i)
    {
      Elogbeta_.row(i) = util::dir_exp_log(lambda_.row(i));
    }

    // Initialize exp( E[log(\beta)] )
    expElogbeta_ = Elogbeta_.array().exp();
  }


  /**
   * updates local variables (gamma)
   * runs coordinate ascent over documents in a dataset
   *
   * @param dataset
   * @param whereCODE
   * @param gamma       WRITES gamma updates to this MatrixXd
   * @param wobsRATIO   wobs/(wobs+who) ratio. default is 1.0 during training.
   *                    can be around 0.5 for testing.
   */
  void update_gamma(util::corpus& dataset,
                    int whereCODE,
                    MatrixXd& gamma,
                    double wobsRATIO = 1.0)
  {

    int number_of_documents;
    std::vector<util::document*> const* dataset_to_process;

    switch(whereCODE)
    {
      case util::FULL:
        number_of_documents = dataset.num_docs_full_;
        dataset_to_process  = &(dataset.docs_full_);
        break;
      case util::VIEW:
        number_of_documents = dataset.num_docs_view_;
        dataset_to_process  = &(dataset.docs_view_);
        break;
      default:
        number_of_documents = 0;
        break;
    }

    // Initialize the variational distribution q(\theta | \gamma)
    std::mt19937_64 gen = util::mersenne_twister();
    gamma = MatrixXd::Zero(number_of_documents, num_topics_);
    for (int j = 0, nRows = gamma.rows(), nCols = gamma.cols(); j < nCols; ++j)
    {
      for (int i = 0; i < nRows; ++i)
      {
        gamma(i, j) = gamma_dist_(gen);
      }
    }

    // Initialize E[log(\theta)]
    MatrixXd Elogtheta = MatrixXd::Zero(number_of_documents, num_topics_);
    for (int i = 0, nRows = Elogtheta.rows(); i < nRows; ++i)
    {
      Elogtheta.row(i) = util::dir_exp_log(gamma.row(i));
    }

    // Initialize exp( E[log(\theta)] )
    MatrixXd expElogtheta = Elogtheta.array().exp();

    // For each document, update its gamma and phi variational parameters
    int doc_length;
    for (int d = 0; d < number_of_documents; ++d)
    {
      util::document const* doc = dataset_to_process->at(d);

      // Compute the number of words to consider in estimating the local
      // latent variables for each document.
      // TRAINING: wobsRATIO should be 1.0
      // TESTING: wobsRATIO is typically 0.5 (used for predictive calculations)
      doc_length = static_cast<int>( std::floor(wobsRATIO * doc->length_) );

      VectorXd gammad        = gamma.row(d);
      VectorXd Elogthetad    = Elogtheta.row(d);
      VectorXd expElogthetad = expElogtheta.row(d);

      // phi_{dwk} is proportional to expElogthetad_k * expElogbetad_w
      // phinorm is the normalizer
      MatrixXd expElogbetad = MatrixXd::Zero(num_topics_, doc_length);
      for (int i = 0; i < doc_length; ++i)
      {
        expElogbetad.col(i) = expElogbeta_.col(doc->words_.at(i));
      }
      VectorXd phinorm = expElogthetad.transpose() * expElogbetad;
      phinorm = phinorm.array() + DBL_MIN;

      // Iterate between gamma and phi until convergence
      double meanchange = INFINITY;
      for (int it = 0; it < e_step_iter_lim_; ++it)
      {
        VectorXd new_gamma;

        // We represent phi implicitly to save memory and time.
        // Substituting the value of the optimal phi back into
        // the update for gamma gives this update. Cf. Lee&Seung 2001.
        VectorXd cts_over_phinorm = VectorXd::Zero(doc_length);
        for (int i = 0; i < doc_length; ++i)
        {
          cts_over_phinorm(i) = doc->counts_.at(i) / phinorm(i);
        }
        VectorXd tmp  = cts_over_phinorm.transpose() * expElogbetad.transpose();
        new_gamma     = alpha_ + expElogthetad.array() * tmp.array();

        Elogthetad    = util::dir_exp_log(new_gamma);
        expElogthetad = Elogthetad.array().exp();
        phinorm       = expElogthetad.transpose() * expElogbetad;
        phinorm       = phinorm.array() + DBL_MIN;

        meanchange = (new_gamma - gammad).norm();
        if (meanchange < e_step_thresh_)
          break;

        gammad = new_gamma;
      }
      gamma.row(d) = gammad;
    }
  }

  /**
   * calculate sufficient statistics
   *
   * call after updating local latent variables (gamma)
   *
   * @param dataset
   * @param gamma
   * @param view_index
   * @param sstats      WRITES sufficient statistics to this MatrixXd
   */
  void calculate_suff_stats(util::corpus& dataset,
                            MatrixXd const& gamma,
                            std::vector<int> const& view_index,
                            MatrixXd& sstats)
   {

    // compute sufficient statistics on util::VIEW
    int number_of_documents = dataset.num_docs_view_;
    std::vector<util::document*>
      const* dataset_to_process = &(dataset.docs_view_);

    // Calculate E[log(\theta)]
    MatrixXd Elogtheta = MatrixXd::Zero(number_of_documents, num_topics_);
    for (int i = 0, nRows = Elogtheta.rows(); i < nRows; ++i)
    {
      Elogtheta.row(i) = util::dir_exp_log(gamma.row(i));
    }

    // Calculate exp( E[log(\theta)] )
    MatrixXd expElogtheta = Elogtheta.array().exp();

    // Initialize sufficient statistics for \lambda
    sstats = MatrixXd::Zero(num_topics_, vocab_length_);

    // Iterate over documents indeced by index into view
    for (auto d : view_index)
    {
      util::document const* doc = dataset_to_process->at(d);

      VectorXd Elogthetad    = Elogtheta.row(d);
      VectorXd expElogthetad = expElogtheta.row(d);

      // phi_{dwk} is proportional to expElogthetad_k * expElogbetad_w
      // phinorm is the normalizer
      MatrixXd expElogbetad = MatrixXd::Zero(num_topics_, doc->length_);
      for (int i = 0; i < doc->length_; ++i)
      {
        expElogbetad.col(i) = expElogbeta_.col(doc->words_.at(i));
      }
      VectorXd phinorm = expElogthetad.transpose() * expElogbetad;
      phinorm = phinorm.array() + DBL_MIN;

      // Contribution of document d to the expected sufficient
      // statistics for the M step.
      VectorXd cts_over_phinorm = VectorXd::Zero(doc->length_);
      for (int i = 0; i < doc->length_; ++i)
      {
        cts_over_phinorm(i) = doc->counts_.at(i) / phinorm(i);
      }

      MatrixXd outerprod = expElogthetad * cts_over_phinorm.transpose();

      for (int i = 0; i < doc->length_; ++i)
      {
        sstats.col(doc->words_.at(i)) += outerprod.col(i);
      }

    }
    sstats = sstats.cwiseProduct(expElogbeta_);

   }

  /**
   * update global latent variables (lambda)
   * @param dataset       input dataset
   * @param num_bootstrap number of bootstraps to use for bumping
   * @param step_size     if 0.0, will use Robbins-Monro
   */
  void update_lambda(util::corpus& dataset,
                     int num_bootstrap,
                     double step_size = 0.0)
  {

    // lambda update is on util::VIEW
    std::vector<util::document*>
      const* dataset_to_process = &(dataset.docs_view_);
    double doc_ratio = (static_cast<double>(dataset.num_docs_full_) /
                       dataset.num_docs_view_);

    // Increment iteration counter (used for step_size calculation)
    ++ current_iter_;

    // Calculate step size
    if (step_size == 0.0)
    {
      rho_t_ = pow( tau0_ + current_iter_, -1*kappa_);
    }
    else
    {
      rho_t_ = step_size;
    }

    // Declare gamma and sufficient statistics matrices
    Eigen::MatrixXd gamma;
    Eigen::MatrixXd sstats;

    // NO BUMPING
    if (num_bootstrap == util::noBOOTSTRAP)
    {
      // Update local latent variables
      update_gamma(dataset, util::VIEW, gamma);

      // Calculate sufficient statistics on all documents in util::VIEW
      std::vector<int> view_index = util::range_index(dataset.num_docs_view_);
      calculate_suff_stats(dataset, gamma, view_index, sstats);

      // Take global latent variable gradient step
      if (rho_t_ == 1.0)
      {
        lambda_ = eta_ + doc_ratio * sstats.array();
      }
      else
      {
        lambda_ = (1.0-rho_t_) * lambda_.array()
                + rho_t_ * (eta_ + doc_ratio * sstats.array());
      }

      // Update internal variables that depend on lambda
      for (int i = 0, nRows = Elogbeta_.rows(); i < nRows; ++i)
      {
        Elogbeta_.row(i) = util::dir_exp_log(lambda_.row(i));
      }
      expElogbeta_ = Elogbeta_.array().exp();

      return;
    }
    // USE BUMPING
    else
    {
      // Initialize bumping related variables
      Eigen::VectorXd bumping_metric(num_bootstrap);
      std::vector<Eigen::MatrixXd> lambda_bumping;
      std::vector<Eigen::MatrixXd> Elogbeta_bumping;
      std::vector<Eigen::MatrixXd> expElogbeta_bumping;

      // Store global parameters from previous iteration
      Eigen::MatrixXd lambda_orig(lambda_);
      Eigen::MatrixXd Elogbeta_orig(Elogbeta_);
      Eigen::MatrixXd expElogbeta_orig(expElogbeta_);

      // Update gamma on documents in VIEW once!
      update_gamma(dataset, util::VIEW, gamma);

      for (int b = 0; b < num_bootstrap; ++b)
      {

        // Calculate bootstrap index. Include all documents in util::VIEW
        // by convention (Tibshirani & Knight 1999)
        std::vector<int> bt_index;
        if (b == 0)
          bt_index = util::range_index(dataset.num_docs_view_);
        else
          bt_index = util::bootstrap_index(dataset.num_docs_view_);

        // Calculate sufficient statistics using bootstrap index
        calculate_suff_stats(dataset, gamma, bt_index, sstats);

        // STORE gradient step in global latent variable lambda
        lambda_ = (1.0-rho_t_) * lambda_orig.array() + rho_t_ * ( eta_ +
                  doc_ratio * sstats.array() );
        lambda_bumping.push_back(lambda_);

        // STORE variables that dpeend on lambda
        for (int i = 0, nRows = Elogbeta_.rows(); i < nRows; ++i)
        {
          Elogbeta_.row(i) = util::dir_exp_log(lambda_.row(i));
        }
        Elogbeta_bumping.push_back(Elogbeta_);
        expElogbeta_ = Elogbeta_.array().exp();
        expElogbeta_bumping.push_back(expElogbeta_);

        // Calculate bumping metric on all documents in util::VIEW.
        // This will use the recently updated global variables above.
        bumping_metric(b) = calculate_log_pred(dataset,
                                               util::VIEW,
                                               util::wobsRATIO);

        // Revert to the global variable settings before taking bumping step
        lambda_      = lambda_orig;
        Elogbeta_    = Elogbeta_orig;
        expElogbeta_ = expElogbeta_orig;

      }

      // Find index that maximizes the bumping  metric
      VectorXd::Index maxIndex;
      bumping_metric.maxCoeff(&maxIndex);
      int maxIndexint = static_cast<double>(maxIndex);

      // Update lambda based on bumped gradient
      lambda_      = lambda_bumping.at(maxIndexint);
      Elogbeta_    = Elogbeta_bumping.at(maxIndexint);
      expElogbeta_ = expElogbeta_bumping.at(maxIndexint);

      return;

    }

  }

  /**
   * calculate evidence lower bound (always re-calculates local variables)
   *
   * @param  dataset
   * @param  whereCODE
   * @return
   */
  double calculate_ELBO(util::corpus& dataset, int whereCODE)
  {

    int number_of_documents;
    std::vector<util::document*> const* dataset_to_process;
    switch(whereCODE)
    {
      case util::FULL:
        number_of_documents = dataset.num_docs_full_;
        dataset_to_process  = &(dataset.docs_full_);
        break;
      case util::VIEW:
        number_of_documents = dataset.num_docs_view_;
        dataset_to_process  = &(dataset.docs_view_);
        break;
      default:
        number_of_documents = 0;
        break;
    }

    double elbo(0.0);

    // Do an E step to update gamma
    Eigen::MatrixXd gamma;
    update_gamma(dataset, whereCODE, gamma);

    // Initialize E[log(\theta)]
    MatrixXd Elogtheta = MatrixXd::Zero(gamma.rows(), gamma.cols());
    for (int i = 0, nRows = Elogtheta.rows(); i < nRows; ++i)
    {
      Elogtheta.row(i) = util::dir_exp_log(gamma.row(i));
    }

    // Initialize exp( E[log(\theta)] )
    MatrixXd expElogtheta = Elogtheta.array().exp();

    // Compute E[log p(dataset | \theta, \beta)]
    for (int d = 0; d < number_of_documents; ++d)
    {
      VectorXd phinorm, tmp, tmp_doc_counts;

      util::document const* doc = dataset_to_process->at(d);

      phinorm        = VectorXd::Zero(doc->length_);
      tmp_doc_counts = VectorXd::Zero(doc->length_);
      for (int i = 0; i < doc->length_; ++i)
      {
        tmp = Elogtheta.row(d).transpose() + Elogbeta_.col(doc->words_.at(i));
        phinorm(i) = util::log_sum_exp(tmp);
        tmp_doc_counts(i) = doc->counts_.at(i);
      }
      elbo += phinorm.dot(tmp_doc_counts);
    }

    // Compute E[log p(\theta | \alpha) - log q(\theta | \gamma)]
    elbo += ((alpha_ - gamma.array()).cwiseProduct(Elogtheta.array())).sum();

    MatrixXd lgammagamma = gamma.unaryExpr(std::ptr_fun(util::lgammad));
    elbo += (lgammagamma.array() - util::lgammad(alpha_)).sum();

    VectorXd gammarowsum = gamma.rowwise().sum();
    gammarowsum = gammarowsum.unaryExpr(std::ptr_fun(util::lgammad));

    elbo += (util::lgammad(alpha_*num_topics_) - gammarowsum.array()).sum();

    // E[log p(\beta | \eta) - log q (\beta | \lambda)]
    elbo += ((eta_ - lambda_.array()).cwiseProduct(Elogbeta_.array())).sum();

    MatrixXd lgammalambda = lambda_.unaryExpr(std::ptr_fun(util::lgammad));
    elbo += (lgammalambda.array() - util::lgammad(eta_)).sum();

    VectorXd lambdarowsum = lambda_.rowwise().sum();
    lambdarowsum = lambdarowsum.unaryExpr(std::ptr_fun(util::lgammad));

    elbo += (util::lgammad(eta_*vocab_length_) - lambdarowsum.array()).sum();

    return elbo;
  }

  /**
   * average per-word log predictive likelihood
   *
   * @param  dataset
   * @param  whereCODE
   * @param  wobsRATIO
   * @return double
   */
  double calculate_log_pred(util::corpus& dataset, int whereCODE,
                              double wobsRATIO)
  {
    double log_pred = 0.0;

    int number_of_documents;
    std::vector<util::document*> const* dataset_to_process;
    switch(whereCODE)
    {
      case util::FULL:
        number_of_documents = dataset.num_docs_full_;
        dataset_to_process  = &(dataset.docs_full_);
        break;
      case util::VIEW:
        number_of_documents = dataset.num_docs_view_;
        dataset_to_process  = &(dataset.docs_view_);
        break;
      default:
        number_of_documents = 0;
        break;
    }

    // Do an E step to update gamma
    Eigen::MatrixXd gamma;
    update_gamma(dataset, whereCODE, gamma, wobsRATIO);

    // Normalize gamma to get E_q[\theta]
    Eigen::MatrixXd Etheta = MatrixXd::Zero(number_of_documents, num_topics_);
    for (int i = 0, nRows = Etheta.rows(); i < nRows; ++i)
    {
      Etheta.row(i) = gamma.row(i).array() / gamma.row(i).sum() ;
    }

    // Normalize lambda to get E_q[\beta]
    Eigen::MatrixXd Ebeta = MatrixXd::Zero(num_topics_, vocab_length_);
    for (int i = 0, nRows = Ebeta.rows(); i < nRows; ++i)
    {
      Ebeta.row(i) = lambda_.row(i).array() / lambda_.row(i).sum() ;
    }

    // Calculate inner product for each w_new in all test documents
    int test_doc_length;
    int normalizer = 0;
    for (int d = 0; d < number_of_documents; ++d)
    {
      util::document const* doc = dataset_to_process->at(d);

      VectorXd Etheta_d = Etheta.row(d);

      test_doc_length = 1 + static_cast<int>(
                              std::floor(wobsRATIO * doc->length_) );

      // NOTE: for loop goes from test_doc_length to doc->length_ !!!
      for (int i = test_doc_length; i < doc->length_; ++i)
      {
        log_pred += log( Etheta_d.dot( Ebeta.col(doc->words_.at(i)) ) );
        ++ normalizer;
      }

    }
    return log_pred / static_cast<double>(normalizer);

  }

}; // class lda


std::ostream& operator<<(std::ostream& os, lda const& lda)
{
  os << "lda.lambda_ = " << std::endl;
  os << lda.lambda_ << std::endl;

  os << "lda.expElogbeta_ = " << std::endl;
  os << lda.expElogbeta_ << std::endl;

  os << "lda.expElogbeta_.rowwise().sum() = " << std::endl;
  os << lda.expElogbeta_.rowwise().sum() << std::endl;

  return os;
}

} // namespace lda

#endif


















