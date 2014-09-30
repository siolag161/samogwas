/****************************** IMLEMENTATION BELOW THIS POINT **************************/
#define NOT_MISSING 1
#define MISSING 0
#define DEFAULT_VALUE = -1

#include "em/naive_bayes_em.hpp"
#include "utils/matrix_utils.hpp"

namespace samogwas
{

/**
 *
 */
void MultiEM::run( ResultEM& result,
                   const Variable& latentVar,
                   const Variables& variables,
                   const Matrix& dataTable,
                   const double threshold,
                   const std::vector< std::vector<bool> > & defTable ) {

  LearnObjectPtrs learnObjs = createLearnObjects(latentVar, variables);
  plMatrixDataDescriptor<int> dataDesc(latentVar ^ variables, &dataTable, &defTable);

  CandidateModels candidateModels;  // vector to store different learners

  for (size_t run = 0; run < nbrRestarts; ++run) {
    // std::cout << "run: " << run << std::endl;
    plComputableObjectList compObjList = createComputableObjects(latentVar, variables); // randomly initialized
    EMLearner learner(compObjList, learnObjs);
    learner.run(dataDesc, threshold); // executes the EM learning process    
    candidateModels.push_back(learner); // puts the learner into the collection
  }

  plEMLearner bestModel = getBestModel( candidateModels, dataDesc);
  impute( result, latentVar, dataTable, bestModel, dataDesc );
  result.jointDistribution = bestModel.get_joint_distribution();

}

void MultiEM::impute( ResultEM& result,                 
                      const plSymbol& latentVar,
                      const Matrix& dataTable,
                      EMLearner& bestModel,
                      plMatrixDataDescriptor<int>& dataDesc )
{
  std::vector<plValues> missingMostProbVals;
  std::vector<std::vector<plProbValue> > missingProbTable;
  bestModel.compute_missing_values_infos(dataDesc, missingMostProbVals, missingProbTable);

  const size_t nbrInds = missingMostProbVals.size(); 
  result.imputedData = std::vector<int> ( missingMostProbVals.size(), 0);

  if (imputMode == ARGMAX) { // argmax imputation
    for (size_t ind = 0; ind < nbrInds; ++ind) {
      result.imputedData[ind] = missingMostProbVals[ind][0];
    }
  } else { // random drawing from the distribution
    for (size_t ind = 0; ind < nbrInds; ++ind) {
      plProbTable probTab(latentVar, missingProbTable[ind], true);
      result.imputedData[ind] = probTab.draw()[0];
    } 
  }
}

} // namespace samogwasends here. samogwas

