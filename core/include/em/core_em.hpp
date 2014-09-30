/****************************************************************************************
 * File: core_em.hpp
 * Description: All EM algorithms must implement this interface.
 *
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 01/08/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_CORE_EM_HPP
#define SAMOGWAS_CORE_EM_HPP

#include <pl.h> // ProBT library
#include <vector>

namespace samogwas
{

typedef plSymbol Variable;
typedef plVariablesConjunction Variables;
typedef plComputableObject DistTab; // Table of Distributions
typedef plComputableObjectList DistTabList;
typedef plJointDistribution JointDistribution;

/** This structure represents a set of attributes that should be returned by
  *  any EM algorithm dedicated to FLTM construction.
  *  Attributes:
  *  - jointDistribution of type ProBT jointDistribution, containing parameters of model
  *    to be learnt by EM (i.e. conditional distribution of each observed variable given
  *    the latent variable and the distribution of the latent variable itself).
  *  - imputedData: a vector of data imputed for the latent variable: this vector contains a value for each observation.
 */
struct ResultEM
{
  ResultEM() { } 
  JointDistribution jointDistribution;
  std::vector<int> imputedData;
};

/** Interface common to by any type of EM underlying algorithm.
 *
 */
struct EMInterface {
  enum ImputationType { ARGMAX = 0, DRAW };
  typedef size_t Cardinality;
  typedef size_t Size; 
  typedef plEMLearner EMLearner;
  typedef std::vector<plLearnObject*> LearnObjectPtrs;
  typedef std::vector<EMLearner> CandidateModels;
  
  typedef std::vector< std::vector<int> > Matrix;

  /** This functor takes as parameters a latent variable, a set of observed variables,
   *  a data matrix and performs an EM algorithm using a stopping threshold.
   *  This threshold provides a stopping criterion: when two consecutive likelihoods
   *  are close enough, the algorithm is said to have converged.
   */
  virtual void operator()( ResultEM& result,
                           const Variable& latentVar,
                           const Variables& variables,
                           const Matrix& dataTable,                   
                           const double threshold );
  
  /** Underlying method for operator()
   *
   */
  virtual void run( ResultEM& result,
                    const Variable& latentVar,
                    const Variables& variables,
                    const Matrix& dataTable,
                    const double threshold );

  virtual ~EMInterface() {}

 protected:

  /** Imputes values for the latent variable based on the estimation made by the EM algorithm.
   *  There are two mode for doing this:
   *    - Drawing randomly from the estimated joint distribution,
   *    - Taking the value that has the maximum likelihood.
   */
  virtual void impute( ResultEM& result,                 
                       const plSymbol& latentVar,
                       const Matrix& dataTable,
                       EMLearner& bestModel,
                       plMatrixDataDescriptor<int>& dataDesc ) = 0;

  /** This run() method that does most of the job. It will be called internally by the run() function above.
   *
   */
  virtual void run( ResultEM& result,
                    const Variable& latentVar,
                    const Variables& variables,
                    const Matrix& dataTable,
                    const double threshold,
                    const std::vector< std::vector<bool> > & defTable ) = 0;

  /** Creates the "learn objects" needed by the EM algorithm method provided by ProBT.
   *
   */
  virtual LearnObjectPtrs createLearnObjects( const Variable& latentVar, const Variables& variables );

  /** Creates the "computable objects" needed by the EM algorithm method provided by ProbBT.
   * 
   */
  virtual plComputableObjectList createComputableObjects( const Variable& latentVar, const Variables& variables );

  /** Takes a reference to the "EM learner" and the data set and returns the logLikelihood.
   *
   */
  virtual double logLikelihood( EMLearner& learner, plMatrixDataDescriptor<int>& dataDesc );

  /** Takes a set of models and returns one that has the best likelihood.
   *
   */
  virtual EMLearner getBestModel( CandidateModels& learners,
                                  plMatrixDataDescriptor<int>& dataDesc );
  /** 
   * The ProBT EM learner requires a parameter called "definition table" which indicates
   * whether the value is missing in the initial data matrix ( missing: 0; non-missing: 1).
   */
  virtual std::vector< std::vector<bool> > createDefinitionTable( const Matrix& dataMat );

};

} // namespace samogwas ends here.

/****************************************************************************************/
#endif // SAMOGWAS_CORE_EM_HPP
