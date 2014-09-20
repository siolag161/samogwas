/****************************************************************************************
 * File: EMCardFunc.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 11/07/2014

 ***************************************************************************************/
#ifndef FLTM_EMCARDFUNC_HPP
#define FLTM_EMCARDFUNC_HPP

namespace samogwas
{

struct CardFunc {
  virtual int operator()(const std::vector<int>& cluster) = 0;
};

/** Computes the cardinality according to the following formula:
 * card = min(alpha * nbrChildren + beta, maxCardinality)
 * where alpha, beta and maxCard are given.
 */
struct LinearCardinality: public CardFunc
{
  LinearCardinality(const double a, const double b, const int card):
      alpha(a), beta(b), maxCard(card) {}
    virtual int operator()(const std::vector<int>& cluster)  {
    int nbrVariables = cluster.size();
    return std::min( int(alpha * nbrVariables + beta), maxCard);
  }

 protected:
  double alpha;
  double beta;
  int maxCard;
};

} // namespace samogwasends here. fltm

/****************************************************************************************/
#endif // FLTM_EMCARDFUNC_HPP