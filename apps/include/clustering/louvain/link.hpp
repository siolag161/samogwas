/****************************************************************************************
 * File: link.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 26/11/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_LOUVAIN_LINK_HPP
#define SAMOGWAS_LOUVAIN_LINK_HPP


#include <memory>
#include <vector>
#include "distance/comparable.hpp"

namespace samogwas
{

namespace louvain {

class WeightedLink: public SimilarityMatrix {
  typedef std::vector<double> Weights;
  typedef std::shared_ptr< std::vector<Weights> > WeightPtr;

 public:
  // WeightedLink( WeightPtr w );
  WeightedLink( const size_t sz );

  virtual double compute( const size_t varA, const size_t varB );
  virtual size_t nbrVariables() const;
  virtual void invalidate() {}

  void setWeight( const size_t varA, const size_t varB, const double val );
 protected:  
  WeightPtr weights;

};

}

} // namespace samogwasends here. samogwas

/****************************************************************************************/
#endif // SAMOGWAS_LINK_HPP
