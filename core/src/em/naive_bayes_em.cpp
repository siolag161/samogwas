/****************************************************************************************
* File: naive_bayes_em.cpp
* Description: Implementation of the EM algorithm dedicated to the parameter learning of
* -----------  the Naive Bayes multinomial model.
*
* @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
* @date: 31/07/2014

***************************************************************************************/
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
void NaiveBayesEM::run( ResultEM& result,
                        const Variable& latentVar,
                        const Variables& variables,
                        const Matrix& dataTable,
                        const double threshold,
                        const std::vector< std::vector<bool> > & defTable ) {

  LearnObjectPtrs learnObjs = createLearnObjects(latentVar, variables);
  plMatrixDataDescriptor<int> dataDesc(latentVar ^ variables, &dataTable, &defTable);

  CandidateModels candidateModels;  // vector to store different learners

  for (size_t it = 0; it < nbrRestarts; ++it) {
    //  joint distribution P(X,Z) = P(Z)*P(X_1 | Z)*...*P(X_n | Z), Z: latent variable
    plComputableObjectList jointDist = createComputableObjects(latentVar, variables); // all distributions
                                                                                      // randomly initialized
    EMLearner learner(jointDist, learnObjs);
    learner.run(dataDesc, threshold); // executes the EM learning process    
    candidateModels.push_back(learner); // puts the learnt model into the collection
  }

  plEMLearner bestModel = getBestModel(candidateModels, dataDesc);
  impute( result, latentVar, dataTable, bestModel, dataDesc );
  result.jointDistribution = bestModel.get_joint_distribution();

}

///////////////////////////////////////////////////////////////////////////////////
void NaiveBayesEM::impute( ResultEM& result,
                           const plSymbol& latentVar,
                           const Matrix& dataTable,
                           EMLearner& bestModel,
                           plMatrixDataDescriptor<int>& dataDesc )
{
  std::vector<plValues> missingValues; /** For each data row X^{i}, stores the value v_i with the highest probability
                                           for the missing variable: v_i = argmax_z P(Z=z|X^{i}). */
  std::vector<std::vector<plProbValue> > latentVarDistri; /** For each data row X^{i}, stores the current probability
                                                        table for the missing variable: P(Z|X^{i}) */

  bestModel.compute_missing_values_infos(dataDesc, missingValues, latentVarDistri);

  const size_t nbrInds = missingValues.size();
  result.imputedData = std::vector<int> ( missingValues.size(), 0);

  if (imputMode == ARGMAX) { // argmax imputation
    for (size_t ind = 0; ind < nbrInds; ++ind) {
      result.imputedData[ind] = missingValues[ind][0];
    }
  } else { // random drawing from the distribution
    for (size_t ind = 0; ind < nbrInds; ++ind) {
      plProbTable probTab(latentVar, latentVarDistri[ind], true);
      result.imputedData[ind] = probTab.draw()[0];
    } 
  }
}

} // namespace samogwas ends here.

