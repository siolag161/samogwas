/****************************************************************************************
 * File: naive_bayes_em.hpp
 * Description: Specific class of EM Algorithm which computes under the Naive-Bayes
 * -----------  multinomial model
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 31/07/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_NAIVE_BAYES_EM_HPP
#define SAMOGWAS_NAIVE_BAYES_EM_HPP

#include "core_em.hpp"
namespace samogwas
{


struct MultiEM: public EMFunc {
  typedef std::vector< std::vector<int> > Matrix;
  
 public:
  MultiEM( int nbrRts, int imMode ): nbrRestarts(nbrRts), imputMode(imMode) {}
  ~MultiEM() {}

  virtual void run( ResultEM& result,
                    const Variable& latentVar,
                    const Variables& variables,
                    const Matrix& dataTable,
                    const std::vector< std::vector<bool> > & defTable,
                    const double threshold );


  virtual void impute( ResultEM& result,                 
                       const plSymbol& latentVar,
                       const Matrix& dataTable,
                       EMLearner& bestModel,
                       plMatrixDataDescriptor<int>& dataDesc );

  int nbrRestarts;
  int imputMode; // methods for imputing missing values  (ARGMAX or DRAW)
  
};
  
} // namespace samogwasends here. samogwas


/****************************************************************************************/
#endif // SAMOGWAS_NAIVE_BAYES_EM_HPP
