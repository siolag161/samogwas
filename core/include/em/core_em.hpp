/****************************************************************************************
 * File: core_em.hpp
 * Description: This provides several common methods and constants used by classes in this
 * -----------  EM module.
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 01/08/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_CORE_EM_HPP
#define SAMOGWAS_CORE_EM_HPP

// #include <boost/filesystem.hpp> // to obtain the program's name
// #include <boost/program_options.hpp>
#include <pl.h>
#include <vector>

namespace samogwas
{

typedef plVariablesConjunction Variables;
typedef plComputableObjectList DistTabList; 
typedef plComputableObject DistTab;
typedef plSymbol Variable;
typedef plJointDistribution JointDistribution;

// /** This structure represents a set of attributes that should be returned by 
//  *  any EM algorithm dedicated to FLTM construction.
//  *  Attributes:
//  *  - varConj Of type proPT variable conjnction, containing observed variables first
//  *       and the latent variable at the end
//  *  - jointDistribution of type proBT jointDistribution, containing parameters of model
//  *       to be learnt by EM (i.e. conditional distribution of observed given
//  the latent variable and the distribution of the latent variable itself.
// */
struct ResultEM
{
  ResultEM() { } 
  JointDistribution jointDistribution;
  std::vector<int> imputedData; 
};

/** Common interface by any type of EM underlying algorithm.
 *
 */
struct EMFunc {
  enum ImputType { ARGMAX = 0, DRAW }; // imputation type, an algorithm
  typedef plVariablesConjunction Variables;
  typedef plSymbol Variable ;
  typedef size_t Cardinality;
  typedef size_t Size; 
  typedef plEMLearner EMLearner;
  typedef std::vector<plLearnObject*> LearnObjectPtrs;
  typedef std::vector<EMLearner> CandidateModels;
  
  typedef std::vector< std::vector<int> > Matrix;

  /** This is the main method, which takes a latent variable, a set of observed variables,
   *  a data sample for those observed variales and performs EM algorithms using the threshold 
   *  provided. The theshold provides the stoppage criteria: whenever two consecutive likelihood
   *  are considered close enough.
   */
  virtual void operator()( ResultEM& result,
                           const Variable& latentVar,
                           const Variables& variables,
                           const Matrix& dataTable,                   
                           const double threshold );
  
  /** Alternative method for operator()
   *
   */
  virtual void run( ResultEM& result,
                    const Variable& latentVar,
                    const Variables& variables,
                    const Matrix& dataTable,
                    const double threshold );

  virtual ~EMFunc() {}

 protected:

  /** Imputes values for the latent variable based on the estimation made by the EM.
   *  There are two mode for doing this:
   *    + Drawing randomly from the estimated joint distribution
   *    + Taking the value that has the maximum likelihood
   */
  virtual void impute( ResultEM& result,                 
               const plSymbol& latentVar,
               const Matrix& dataTable,
               EMLearner& bestModel,
               plMatrixDataDescriptor<int>& dataDesc ) = 0;

  /** The method that does most of the job. Will be called internally by the run above
   *
   */
  virtual void run( ResultEM& result,
                    const Variable& latentVar,
                    const Variables& variables,
                    const Matrix& dataTable,
                    const std::vector< std::vector<bool> > & defTable,
                    const double threshold ) = 0;  

  /** Creates learning_objects needed by em_algorithm method provided by ProBT
   *
   */
  virtual LearnObjectPtrs createLearnObjects( const Variable& latentVar, const Variables& variables );

  /** Sets up Computatble Objects needed by em_algorithm by ProbBT
   * 
   */
  virtual plComputableObjectList createComputableObjects( const Variable& latentVar, const Variables& variables );

  /** Takes a refenece to the current learner and the actual dataset and returns the current logLikelihood
   *
   */
  virtual double logLikelihood( EMLearner& learner, plMatrixDataDescriptor<int>& dataDesc );

  /** Take a set of models and returns one that having the best likelihood
   *
   */
  virtual EMLearner getBestModel( CandidateModels& learners,
                                  plMatrixDataDescriptor<int>& dataDesc );
  /** 
   *
   */
  virtual std::vector< std::vector<bool> > createDefinitionTable( const Matrix& dataMat );

};



} // namespace fltm ends here. fltm

/****************************************************************************************/
#endif // SAMOGWAS_CORE_EM_HPP
