/****************************************************************************************
 * File: test_common.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 26/11/2014

 ***************************************************************************************/
#ifndef LOUVAIN_TEST_COMMON_HPP
#define LOUVAIN_TEST_COMMON_HPP

#include "distance/comparable.hpp"
#include <vector>

struct Simi: public samogwas::SimilarityMatrix {
  Simi(std::vector< std::vector<double> > d): sim(d) {}
  virtual double compute( const size_t varA, const size_t varB ) { return sim[varA][varB]; }

  virtual size_t nbrVariables() const {
    return sim.size();
  }

  virtual void invalidate() {}

  std::vector< std::vector<double> > sim;
};

/****************************************************************************************/
#endif // _TEST_COMMON_HPP
