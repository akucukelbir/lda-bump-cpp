/*
 *  Copyright (c) 2014 Alp Kucukelbir, alp@cs.columbia.edu
 *
 *  - Adapted from code by Cheng Zhang, Xavi Gratal, Chong Wang,
 *    Matthew D. Hoffman, and David M. Blei
 *  - Licensed under terms of GPLv2 license (see LICENSE)
 *
 */

#ifndef corpus_hpp
#define corpus_hpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <string>
#include <cassert>
#include <random>

#include "util.hpp"

namespace util
{

/**
 * class vocabulary
 */
class vocabulary
{
public:
  std::vector<std::string> vocab_;
  int length_;

  vocabulary(): length_(0) {};

  void read_data(std::string filename)
  {
    vocab_.clear();
    std::ifstream input(filename.c_str());
    std::copy(std::istream_iterator<std::string>(input),
          std::istream_iterator<std::string>(),
          std::back_inserter(vocab_));
    length_ = vocab_.size();
  }
}; // vocabulary

std::ostream& operator<<(std::ostream& os, vocabulary const& obj)
{
  // write obj to stream
  os << "vocabulary.vocab_ =" << std::endl;
  int i = 0;
  for (auto vv : obj.vocab_)
  {
    os << i << " : " << vv << std::endl;
    ++i;
  }
  return os;
};

/**
 * class document
 */
class document
{
public:
  std::vector<int> words_;      // vector of words
  std::vector<int> counts_;     // vector of counts
  int length_;                  // number of words
  int total_;                   // sum of word counts

  std::vector<int> words_bt_;   // bootstrapped vector of words
  std::vector<int> counts_bt_;  // bootstrapped vector of counts
  int length_bt_;               // number of bootstrapped words
  int total_bt_;                // sum of bootstrapped word counts

  document(): length_(0),   total_(0),
              length_bt_(0), total_bt_(0) {};


  void shuffle_words_and_counts()
  {
    auto seed = unsigned(std::time(0));

    std::srand(seed);
    std::random_shuffle(words_.begin(), words_.end());

    std::srand(seed);
    std::random_shuffle(counts_.begin(), counts_.end());
  }

}; // document

std::ostream& operator<<(std::ostream& os, document const& doc)
{
  // write doc to stream
  os << "document.words_ = " << std::endl;
  for (auto ww : doc.words_)
  {
    os << ww << std::endl;
  }
  os << "document.counts_ = " << std::endl;
  for (auto cc : doc.counts_)
  {
    os << cc << std::endl;
  }
  os << "document.length_ = " << doc.length_ << std::endl;
  os << "document.total_ = " << doc.total_ << std::endl;

  os << "document.words_bt_ = " << std::endl;
  for (auto ww : doc.words_bt_)
  {
    os << ww << std::endl;
  }
  os << "document.counts_bt_ = " << std::endl;
  for (auto cc : doc.counts_bt_)
  {
    os << cc << std::endl;
  }
  os << "document.length_bt_ = " << doc.length_ << std::endl;
  os << "document.total_bt_ = " << doc.total_ << std::endl;
  return os;
};

/**
 * class corpus
 */
class corpus
{
public:
  std::vector<document> docs_;         // all the documents in the corpus
  int num_docs_;                       // # of documents in the corpus

  std::vector<document*> docs_full_;   // full corpus (convenience object)
  int num_docs_full_;                  // # of documents in the corpus

  std::vector<document*> docs_view_;   // arbitrary view into corpus
  int num_docs_view_;                  // # of docs in view

  corpus():
    num_docs_full_(0),
    num_docs_view_(0) { };

  corpus(std::vector<document> docs):
    docs_(docs),
    num_docs_full_(0),
    num_docs_view_(0) { init (); };

  corpus(util::corpus const& corpus):
    docs_(corpus.docs_),
    num_docs_(corpus.num_docs_),
    num_docs_full_(0),
    num_docs_view_(0) { init(); };

  void init()
  {
    // copy (pointers of) full corpus into docs_full_
    for (int i = 0; i < num_docs_; ++i)
    {
      docs_full_.push_back( &docs_[i] );
    }
    num_docs_full_ = docs_.size();

    // copy full into view
    copy_docs_into_view();
  }

  void corpus_clear()
  {
    num_docs_full_ = 0;
    docs_full_.clear();
  }

  void view_clear()
  {
    num_docs_view_ = 0;
    docs_view_.clear();
  }


  /**
   * samples documents into view with replacement (bootstrap)
   * @param number_of_samples  if 0, use size of full
   */
  std::vector<int> sample_from_docs_into_view(int number_of_samples = 0)
  {
    view_clear();

    if (number_of_samples == 0)
    {
      number_of_samples = num_docs_full_;
    }

    std::vector<int> sample_indices;

    // sample with replacement (bootstrap)
    std::mt19937_64 gen = util::mersenne_twister();
    std::uniform_int_distribution<> unif(0, num_docs_full_ - 1 );
    int index;
    for (int i = 0; i < number_of_samples; ++i)
    {
      index = unif(gen);
      sample_indices.push_back(index);
      docs_view_.push_back( &docs_[index] );
    }
    num_docs_view_ = number_of_samples;

    return sample_indices;
  }

  /**
   * copy documents into view, no resampling
   */
  void copy_docs_into_view()
  {
    view_clear();

    for (int i = 0; i < num_docs_full_; ++i)
    {
      docs_view_.push_back( docs_full_[i] );
    }
    num_docs_view_ = num_docs_full_;
  }


  /**
   * read in data file stored in LDA-C format
   * @param filename  string containing filename
   */
  void read_data(std::string filename)
  {
    std::cout << "[ corpus::read_data ]" << std::endl;
    std::cout << "  processing file " << "\"" << filename << "\"" << std::endl;

    corpus_clear();

    std::ifstream input(filename.c_str());
    std::string line;

    while (getline(input, line))
    {
      document doc_now;
      int max_word = 0; // the max word id of this document
      std::vector<std::string> tokens;
      std::istringstream iss(line);
      int tt = 0;
      int num_terms;
      do
      {
        std::string sub;
        iss >> sub; // get a token with space as delim
        if (tt == 0)
        {
          num_terms = atoi(sub.c_str());
          // std::cout << "num_terms: " << num_terms << std::endl;
        }
        else
        {
          std::string token; //get every number
          std::stringstream stream(sub);
          int ii = 0;
          while (getline(stream, token, ':'))
          {
            if (ii == 0)
            {
              doc_now.words_.push_back(atoi(token.c_str()));
            }
            else
            {
              doc_now.counts_.push_back(atof(token.c_str()));
              doc_now.total_ += atoi(token.c_str());
              doc_now.length_ ++;
            }
            ii++;
          }
        }
        tt++;
      }
      while (iss);

      assert( num_terms == tt-1 );
      assert( doc_now.length_ == num_terms );
      assert( doc_now.words_.size() == doc_now.counts_.size() );

      // shuffle the order of words so that log pred calculation isn't biased
      // (this should've been done during generation of LDA-C format file)
      doc_now.shuffle_words_and_counts();

      // push back the doc to the corpus
      docs_.push_back(doc_now);

    }

    num_docs_ = docs_.size();
    std::cout << "  number of docs in corpus  : " << num_docs_ << std::endl;
    std::cout << std::endl;
  }

}; // corpus

std::ostream& operator<<(std::ostream& os, corpus const& corp)
{
  // write corpus to stream
  os << "-=-=-=-=-=-=- corpus operator<< -=-=-=-=-=-=-" << std::endl;
  os << "corpus.docs_full_ = " << std::endl;
  for (auto ww : corp.docs_full_)
  {
    os << *ww << std::endl;
  }
  os << "corpus.num_docs_full_ = " << corp.num_docs_full_ << std::endl;
  os << "---------------------------------------------" << std::endl;
  os << "corpus.docs_view_ = " << std::endl;
  for (auto ww : corp.docs_view_)
  {
    os << *ww << std::endl;
  }
  os << "corpus.num_docs_view_ = " << corp.num_docs_view_ << std::endl;
  os << "-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-" << std::endl;
  return os;
}

} // namespace util

#endif
