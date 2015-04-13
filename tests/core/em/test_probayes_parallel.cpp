// #define BOOST_TEST_DYN_LINK
// #ifdef STAND_ALONE
// #   define BOOST_TEST_MODULE
// #endif  
// #include <boost/test/unit_test.hpp>
// #include <fstream>
// #include <map>
// #include <boost/random.hpp>
// #include <boost/random/normal_distribution.hpp>
// #include <math.h>
// #include <boost/lexical_cast.hpp>

// #include <thread>

// #include <pl.h>
// #define NOT_MISSING 1
// #define MISSING 0
// #define MISSING_VALUE -1
// #define MISSING_OBS 10

// typedef std::vector< std::vector<bool>> DefTab;
// typedef std::shared_ptr<DefTab> DefTabPtr;
// typedef std::vector<std::vector<int>> Matrix;

// typedef plSymbol Variable;
// typedef plVariablesConjunction Variables;
// typedef plComputableObject DistTab; // Table of Distributions
// typedef plComputableObjectList DistTabList;
// typedef plJointDistribution JointDistribution;
// typedef size_t Size; 
// typedef plEMLearner EMLearner;
// typedef std::vector<plLearnObject*> LearnObjectPtrs;
// typedef std::vector<EMLearner> CandidateModels;

// class Data 
// { 
// };

// // Simulate and generate a 2-components (kernels) Gaussian mixture data. This data will be used for learning using EM based on a BIC score


// BOOST_FIXTURE_TEST_SUITE( Test_EM_Parallel, Data )


// ////////////////////////////////////////////

// void run_em( int tid ) {
//   plSymbol latentVar( boost::lexical_cast<std::string>(tid), plIntegerType(0, 1) );
//   plVariablesConjunction vars;
//   vars ^= plSymbol( "X_0", plIntegerType(0, 2) );  vars ^= plSymbol( "X_1", plIntegerType(0, 2) );  vars ^= plSymbol( "X_2", plIntegerType(0, 2) );

//   LearnObjectPtrs learnObjects;
//   plLearnObject* learnLatent = new plLearnHistogram(latentVar);
//   learnObjects.push_back(learnLatent);  
//   for (size_t var = 0; var < vars.size(); ++var) {
//     learnObjects.push_back(new plCndLearnObject <plLearnHistogram> (vars[var], latentVar));
//   }

//   plComputableObjectList jointDist; // joint distribution P(X,Z) = P(Z)*P(X_1 | Z)*...*P(X_n | Z)
//   const plProbTable latentProbInit(latentVar, true); // Z denotes the latent variable.
//   jointDist *= latentProbInit;  
//   for (size_t x = 0; x < vars.size(); ++x) {
//     plDistributionTable distTab_Xi_Z(vars[x], latentVar); // Conditional distribution P(Xi | Z)
//     for (size_t h = 0; h < latentVar.cardinality(); ++h) {
//       distTab_Xi_Z.push( plProbTable(vars[x], true), (int)h );
//     } 
//     jointDist *= distTab_Xi_Z; // adds the conditional table
//   }
    
//   auto learner = new EMLearner(jointDist, learnObjects);
// }

// static const int NBR_RESTARTS = 5;
// BOOST_AUTO_TEST_CASE( Test_EM_Parallel ) {

//   auto dat = std::make_shared<Matrix>(std::initializer_list<std::vector<int>>{
//       std::vector<int>{ MISSING_VALUE,1,1,1},
//       std::vector<int>{ MISSING_VALUE,0,1,2},
//       std::vector<int>{ MISSING_VALUE,0,1,2},
//       std::vector<int>{ MISSING_VALUE,2,2,0} });
//   auto def = std::make_shared<DefTab>(std::initializer_list< std::vector<bool> >{
//     std::vector<bool>{MISSING, NOT_MISSING, NOT_MISSING, NOT_MISSING},
//     std::vector<bool>{MISSING, MISSING, NOT_MISSING, NOT_MISSING},
//     std::vector<bool>{MISSING, NOT_MISSING, NOT_MISSING, NOT_MISSING},
//     std::vector<bool>{MISSING, NOT_MISSING, NOT_MISSING, NOT_MISSING} });


//   std::thread threads[NBR_RESTARTS];
//   for(int i = 0; i < NBR_RESTARTS; ++i)
//   {   
//     threads[i] = std::thread(run_em, i);
//   }

//   for (int i = 0; i < NBR_RESTARTS; ++i) {
//     threads[i].join();
//   }
// }

// BOOST_AUTO_TEST_SUITE_END()  /// Test InfoTheo ends here
