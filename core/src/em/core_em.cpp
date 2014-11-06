/****************************************************************************************
* File: core_em.cpp
* Description: Implementation of the core_em.hpp, the common interface for all EM algorithms.
*
* @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
* @date: 31/07/2014

***************************************************************************************/

#define NOT_MISSING 1
#define MISSING 0
#define DEFAULT_VALUE -1
#define MISSING_VALUE -1

#include "em/core_em.hpp"
#include "utils/matrix_utils.hpp"

namespace samogwas
{

/////////////////////////////////////////////////////////////////////////////////////
void EMInterface::run( ResultEM& result,
                       const Variable& latentVar,
                       const Variables& variables,
                       const Matrix& dataTable,
                       const double threshold ) {
  std::vector< std::vector<bool> >* defTable = createDefinitionTable( dataTable ); 
  run( result, latentVar, variables, dataTable, threshold, defTable );
}

/////////////////////////////////////////////////////////////////////////////////////
EMInterface::LearnObjectPtrs EMInterface::createLearnObjects( const Variable& latentVar,
                                                              const Variables& variables )
{
  LearnObjectPtrs learnObjects;
  plLearnObject* learnLatent = new plLearnHistogram(latentVar);
  learnObjects.push_back(learnLatent);
  
  for (size_t var = 0; var < variables.size(); ++var) {
    learnObjects.push_back(new plCndLearnObject <plLearnHistogram> (variables[var], latentVar));
  }

  return learnObjects;
}

/////////////////////////////////////////////////////////////////////////////////////
plComputableObjectList EMInterface::createComputableObjects( const Variable& latentVar,
                                                             const Variables& variables )
{
  plComputableObjectList jointTable; // joint distribution P(X,Z) = P(Z)*P(X_1 | Z)*...*P(X_n | Z)

  const plProbTable latentProbInit(latentVar, true); // Z denotes the latent variable.
  jointTable *= latentProbInit;
  
  for (size_t x = 0; x < variables.size(); ++x) {
    plDistributionTable distTab_Xi_Z(variables[x], latentVar); // Conditional distribution P(Xi | Z)
    for (size_t h = 0; h < latentVar.cardinality(); ++h) {
      distTab_Xi_Z.push( plProbTable(variables[x], true), (int)h );
    }
    jointTable *= distTab_Xi_Z; // adds the conditional table
  }  
  
  return jointTable;
}
////////////////////////////////////////////////////////////////////

double EMInterface::logLikelihood( EMLearner& learner, plMatrixDataDescriptor<int>& dataDesc )
{
  double result = 0.0;
  try {
    result = learner.get_last_computed_loglikelihood(); // if stored value available
  } catch (plError& e) {
    result = (double)learner.compute_loglikelihood(dataDesc); // otherwise, computes the value
  }
  return result;
}

/////////////////////////////////////////////////////////////////////
EMInterface::EMLearner EMInterface::getBestModel( CandidateModels& learners,
                                                  plMatrixDataDescriptor<int>& dataDesc )
{
  EMLearner bestModel = learners[0];  
  for (size_t i = 1; i < learners.size(); ++i) {
    if ( logLikelihood( learners[i], dataDesc) > logLikelihood(bestModel, dataDesc)) {
      bestModel = learners[i];
    }
  }  
  return bestModel;
}

////////////////////////////////////////////////////////////////////
/** ProBT requires that the latent variable appears as the first variable (i.e. the first column)
  * in the definition table. The column corresponding to the latent variable is entirely missing.
  *
  */
std::vector< std::vector<bool> >* EMInterface::createDefinitionTable( const Matrix& dataMat ) {
    
  const Size nbrInds = utility::nrows(dataMat);
  const Size nbrVars = utility::ncols(dataMat);
  auto defTable = new std::vector<std::vector<bool>>(nbrInds, std::vector<bool>(nbrVars, NOT_MISSING));
  for (size_t ind = 0; ind < nbrInds; ++ind) { // Each row will contain (FALSE TRUE TRUE ... TRUE).
    for ( size_t v = 0; v < nbrVars; ++v ) {
      if (dataMat[ind][v]== MISSING_VALUE)
        (*defTable)[ind][v] = MISSING;
    }
    
  }
   return defTable;
} 

}
