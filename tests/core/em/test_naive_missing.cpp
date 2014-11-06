#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
#   define BOOST_TEST_MODULE
#endif  
#include <boost/test/unit_test.hpp>
#include <pl.h>
#include <fstream>
#include <iostream>

#include "em/core_em.hpp"

class Data 
{ 
};


static const int CARDI = 2;
static const plProbValue PR[] = {0.3, 0.7}; // actual

std::vector< std::vector< std::vector<float> > > PROBTABS = {
  { {0.5,0.5}, {0.2,0.8} },
  { {0.7,0.3}, {0.4,0.6} }
};

static const int NOT_MISSING = 1;
static const int MISSING = 0;

typedef std::vector<int> Row;
typedef std::vector<Row> Matrix;
typedef plMatrixDataDescriptor<int> MatrixDesc;

using namespace samogwas;
/////////////////////////////////////////////////////////////////////

BOOST_FIXTURE_TEST_SUITE( Test_Naive_Missing_Data, Data )

Matrix generate_data(unsigned n);
void run_em( Matrix& data,
             unsigned int &nparams, plFloat &llk,
             plFloat &bic, plJointDistribution &model );

// MatrixDesc create_matrix_desc( Matrix& data );



typedef plSymbol Variable;
typedef plVariablesConjunction Variables;
//////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE( Test_EM_bool ) {
  auto data = generate_data(1000); 

  
  unsigned int nparams; // The number of parameters
  plFloat llk; // The log-likelihood
  plFloat bic; // The bic score
  plJointDistribution model; // the learnt model
  
  run_em( data, nparams, llk, bic, model );

  // std::cout << model.get_computable_object_list()[0] << std::endl;
  // std::cout << model.get_computable_object_list()[1] << std::endl;
}

BOOST_AUTO_TEST_CASE( Test_EM_bunbun ) {
  // auto dataTable = generate_data(1000);
  // auto defTable = EMInterface::createDefinitionTable(dataTable);

  // for (int i = 0; i < dataTable.size();++i) {
  //   if (i%20==0) {
  //     (*defTable)[i][i] =MISSING;
  //     dataTable[i][i] = -1;
  //   }
  // }   

  
  // const plSymbol latentVar("C", plIntegerType(0, CARDI-1));
  // plVariablesConjunction variables;
  // variables ^= plSymbol("X", plIntegerType(0, CARDI-1));
  // variables ^= plSymbol("Y", plIntegerType(0, CARDI-1));
  
  // EMInterface::LearnObjectPtrs learnObjs = EMInterface::createLearnObjects(latentVar, variables);
  // plMatrixDataDescriptor<int> dataDesc(latentVar ^ variables, &dataTable, defTable);

  // // CandidateModels candidateModels;  // vector to store different learners

  // // for (size_t it = 0; it < nbrRestarts; ++it) {
  //   //  joint distribution P(X,Z) = P(Z)*P(X_1 | Z)*...*P(X_n | Z), Z: latent variable
  //   plComputableObjectList jointDist = EMInterface::createComputableObjects(latentVar, variables); // all distributions
  //   // randomly initialized
  //   EMInterface::EMLearner learner(jointDist, learnObjs);
  //   learner.set_same_missing_variables(false);
  //   learner.run(dataDesc, 0.0000000001); 
    // executes the EM learning process    
    // candidateModels.push_back(learner); // puts the learnt model into the collection
  // }

  // plEMLearner bestModel = getBestModel(candidateModels, dataDesc);
  // imputeLatent( result, latentVar, dataTable, bestModel, dataDesc );
  // result.jointDistribution = bestModel.get_joint_distribution();

}


///////////////////////////////////////////////////////////////////
Matrix generate_data(unsigned nrows) {
  
  const plSymbol C("C", plIntegerType(0, CARDI-1));
  const plProbTable PC(C, PR);

  plVariablesConjunction vars;
  vars ^= plSymbol("X", plIntegerType(0, CARDI-1));
  vars ^= plSymbol("Y", plIntegerType(0, CARDI-1));

  plComputableObjectList jointTable;
   jointTable *= PC;
  for (int v = 0; v < vars.size(); ++v) {
    plDistributionTable pv(vars[v], C);
    for(unsigned int i = 0; i < CARDI; ++i) { 
      pv.push( plProbTable( vars[v], PROBTABS[v][i]), (int)i );
    }
    jointTable *= pv;
  }
  
  
  // Constructing the model
  //
  plVariablesConjunction varConj = C^vars;
  const plJointDistribution mixture(varConj, jointTable);

  plValues vals(varConj);

  Matrix data(nrows, Row(varConj.size(), 0.0));

  for(unsigned int i = 0; i < nrows; ++i) {
    mixture.draw(vals);
    data[i][0] = -1;
    
    for (int v = 0; v < vars.size(); ++v) {
      data[i][v+1] = vals[vars[v]];
    }
  }

  return data;
}


void run_em( Matrix& data,
             unsigned int &nparams, plFloat &llk,
             plFloat &bic, plJointDistribution &model )
{
  const plSymbol C("C", plIntegerType(0, CARDI-1));
  
  plVariablesConjunction vars;
  vars ^= plSymbol("X", plIntegerType(0, CARDI-1));
  vars ^= plSymbol("Y", plIntegerType(0, CARDI-1));

  // EM Initial distribution on the class (kernel) variable: P(C)
  const bool random_prob = true;
  const plProbTable pc_init(C, random_prob);
  
  plVariablesConjunction varConj = C^vars;
  std::vector< std::vector<bool> >* def =
      new std::vector< std::vector<bool> >( data.size(),
                                            std::vector<bool>(varConj.size(), NOT_MISSING) );
  for (int i = 0; i < data.size();++i) {
    (*def)[i][0] = MISSING;
    if (i%10==0) {
      data[i][2] = -1;
      (*def)[i][2] = MISSING;
    }
  }   
    
  plComputableObjectList jointTable = EMInterface::createComputableObjects(C,vars); 

  // plComputableObjectList jointTable;
  //  jointTable *= pc_init;
  // for (int v = 0; v < vars.size(); ++v) {
  //   plDistributionTable pv(vars[v], C);
  //   for(unsigned int i = 0; i < CARDI; ++i) {
  //     pv.push( plProbTable( vars[v], random_prob), (int)i );
  //   }
  //   jointTable *= pv;
  // }

  // std::vector<plLearnObject*> learn_objs(varConj.size());
  // learn_objs[0] = new plLearnHistogram(C);
  // for (int i=0; i<vars.size();++i)
  //   learn_objs[i+1] = new plCndLearnObject<plLearnHistogram>(vars[i],C);;

  std::vector<plLearnObject*> learn_objs = EMInterface::createLearnObjects(C,vars);

  MatrixDesc dataDesc(varConj, &data, def); // error
  plEMLearner myEM( jointTable, learn_objs);
  // Run untill convergence
  myEM.run( dataDesc, 0.0001);

  // Fill the output parameters
  nparams = myEM.get_n_parameters();
  llk = myEM.get_last_computed_loglikelihood();
  bic = llk - 0.5*nparams*std::log( dataDesc.get_n_records());
  model = myEM.get_joint_distribution();
}

////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END() 
