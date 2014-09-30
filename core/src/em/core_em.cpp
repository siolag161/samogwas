#define NOT_MISSING 1
#define MISSING 0
#define DEFAULT_VALUE = -1

#include "em/core_em.hpp"
#include "utils/matrix_utils.hpp"


namespace samogwas
{
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EMInterface::operator()( ResultEM& result,
                         const Variable& latentVar,
                         const Variables& variables,
                         const Matrix& dataTable,                   
                         const double threshold ) {
  run( result, latentVar, variables, dataTable, threshold );
}

void EMInterface::run( ResultEM& result,
                         const Variable& latentVar,
                         const Variables& variables,
                         const Matrix& dataTable,                   
                         const double threshold ) {
  std::vector< std::vector<bool> > defTable = createDefinitionTable( dataTable ); 
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

////////////////////////////////////////////////////////////////////
plComputableObjectList EMInterface::createComputableObjects( const Variable& latentVar,
                                                        const Variables& variables )
{
  plComputableObjectList dataCndTable; 

  const plProbTable latentProbInit(latentVar, true);
  dataCndTable *= latentProbInit;
  
  for (size_t x = 0; x < variables.size(); ++x) {
    plDistributionTable distTab_Xi_Z(variables[x], latentVar); 
    for (size_t h = 0; h < latentVar.cardinality(); ++h) {
      distTab_Xi_Z.push( plProbTable(variables[x], true), (int)h );
    }
    dataCndTable *= distTab_Xi_Z; // adds the conditional table to result     

  }  
  
  return dataCndTable;
}
////////////////////////////////////////////////////////////////////

double EMInterface::logLikelihood( EMLearner& learner, plMatrixDataDescriptor<int>& dataDesc )
{
  double result = 0.0;
  try {
    result = learner.get_last_computed_loglikelihood();
  } catch (plError& e) {
    result = (double)learner.compute_loglikelihood(dataDesc);
  }
  return result;
}

/////////////////////////////////////////////////////////////////////
EMInterface::EMLearner EMInterface::getBestModel( CandidateModels& learners,
                        plMatrixDataDescriptor<int>& dataDesc )
{
  EMLearner bestModel = learners[0];  
  for (size_t i = 0; i < learners.size(); ++i) {
    if ( logLikelihood( learners[i], dataDesc) > logLikelihood(bestModel, dataDesc)) {
      bestModel = learners[i];
    }
  }  
  return bestModel;
}

std::vector< std::vector<bool> > EMInterface::createDefinitionTable( const Matrix& dataMat ) {
    
  const Size nbrInds = utility::nrows(dataMat);
  const Size nbrVars = utility::ncols(dataMat);
  
  std::vector< std::vector<bool> > defTable;
  defTable.reserve( nbrInds );
  for (size_t ind = 0; ind < nbrInds; ++ind) {
    defTable.push_back(std::vector<bool>(nbrVars, NOT_MISSING));
    defTable[ind][0] = false; // the first columns = FALSE
  }

  return defTable;
} 

}
