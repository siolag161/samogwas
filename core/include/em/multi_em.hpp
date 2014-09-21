/****************************************************************************************
 * File: multi_em.hpp
 * Description: 
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 03/08/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_MULTI_EM_HPP
#define SAMOGWAS_MULTI_EM_HPP

#include "core_em.hpp"

namespace samogwas
{


struct NV_EM: public EMFunc {
  typedef std::vector< std::vector<int> > Matrix;  
  

 public:
  NV_EM( int nbrRts, int imMode ): nbrRestarts(nbrRts), imputMode(imMode) {
    clear_values();
  }
  
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
  // protected:
  virtual void clear_values();
  virtual void init_values( int N, int K, const Variables& variables );

  virtual void step_E( std::vector< std::vector<double> >& theta,
                       const std::vector<double>& pY,
                       const std::vector< std::vector< std::vector<double> > >& pYX,
                       const NV_EM::Matrix& data );

  virtual void step_M( std::vector<double>& pY,
                       std::vector< std::vector< std::vector<double> > >& pYX,
                       const std::vector< std::vector<double> >& theta,
                       const NV_EM::Matrix& data );
  
  // private:
  int nbrRestarts;
  int imputMode; // methods for imputing missing values  (ARGMAX or DRAW)

  std::vector< std::vector<double> > theta;
  std::vector<double> pY;
  std::vector< std::vector< std::vector<double> > > pYX;
  
};

} // namespace samogwasends here. samogwas

/****************************************************************************************/
#endif // SAMOGWAS_MULTI_EM_HPP

