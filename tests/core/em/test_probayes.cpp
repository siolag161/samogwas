#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
#   define BOOST_TEST_MODULE
#endif  
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <map>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <math.h>
#include <boost/lexical_cast.hpp>
#include <thread>

#include <pl.h>

#define NOT_MISSING 1
#define MISSING 0
#define MISSING_VALUE -1
#define MISSING_OBS 10

typedef std::vector< std::vector<bool>> DefTab;
typedef std::shared_ptr<DefTab> DefTabPtr;
typedef std::vector<std::vector<int>> Matrix;

typedef plSymbol Variable;
typedef plVariablesConjunction Variables;
typedef plComputableObject DistTab; // Table of Distributions
typedef plComputableObjectList DistTabList;
typedef plJointDistribution JointDistribution;
typedef size_t Size; 
typedef plEMLearner EMLearner;
typedef std::vector<plLearnObject*> LearnObjectPtrs;
typedef std::vector<EMLearner> CandidateModels;

class Data 
{};

BOOST_FIXTURE_TEST_SUITE( Test_EM_With_Missing_Data, Data )
////////////////////////////////////////////
static const int NBR_RESTARTS = 10;
BOOST_AUTO_TEST_CASE( Test_EM_Parallel ) {

  auto dat = std::make_shared<Matrix>(std::initializer_list<std::vector<int>>{
      std::vector<int>{ MISSING_VALUE,1,1,1},
      std::vector<int>{ MISSING_VALUE,0,1,2},
      std::vector<int>{ MISSING_VALUE,0,1,2},
      std::vector<int>{ MISSING_VALUE,2,2,0} });
  auto def = std::make_shared<DefTab>(std::initializer_list< std::vector<bool> >{
    std::vector<bool>{MISSING, NOT_MISSING, NOT_MISSING, NOT_MISSING},
    std::vector<bool>{MISSING, NOT_MISSING, NOT_MISSING, NOT_MISSING},
    std::vector<bool>{MISSING, NOT_MISSING, NOT_MISSING, NOT_MISSING},
    std::vector<bool>{MISSING, NOT_MISSING, NOT_MISSING, NOT_MISSING} });  

  
  #pragma omp parallel for schedule(static,1) num_threads(7) 
  for( int i = 0; i < NBR_RESTARTS; ++i )
  {
    plSymbol latentVar( boost::lexical_cast<std::string>("Z"), plIntegerType(0, 1) );
    plVariablesConjunction vars;
    vars ^= plSymbol( "X_0", plIntegerType(0, 2) );
    vars ^= plSymbol( "X_1", plIntegerType(0, 2) );
    vars ^= plSymbol( "X_2", plIntegerType(0, 2) );
    //plError::treat_this_warning_as_error(62);
    LearnObjectPtrs learnObjects;
    plLearnObject* learnLatent = new plLearnHistogram(latentVar);
    learnObjects.push_back(learnLatent);  
    for (size_t var = 0; var < vars.size(); ++var) {
      learnObjects.push_back(new plCndLearnObject <plLearnHistogram> (vars[var], latentVar));
    }

    plComputableObjectList jointTab; 
    const plProbTable latentProbInit(latentVar, true); 
    jointTab *= latentProbInit;  
    for (size_t x = 0; x < vars.size(); ++x) {
      plDistributionTable distTab_Xi_Z(vars[x], latentVar); 
      for (size_t h = 0; h < latentVar.cardinality(); ++h) {
        distTab_Xi_Z.push( plProbTable(vars[x], true), (int)h );
      } 
      jointTab *= distTab_Xi_Z; 
    }
    
     EMLearner learner (jointTab, learnObjects);
     plMatrixDataDescriptor<int> dataDesc( latentVar^vars, dat.get(), def.get() );
     learner.run(dataDesc, 0.001); // executes the EM learning process       
  }
}
 








// //========================================================
// void output_mixture(const plJointDistribution &mixture)
// {
//   std::cout << mixture.get_computable_object_list()[0] << std::endl;
//   std::cout << mixture.get_computable_object_list()[1] << std::endl;
// }

// std::vector<std::vector<double>>* generate_data(unsigned int ndata);
// std::vector< std::vector<bool>>* generate_data_def(unsigned int ndata);
// plDataDescriptor* get_data_desc(unsigned int ndata, const plSymbol& C, const plSymbol& X);


// void run_em(  const unsigned int n_data, const std::vector<unsigned int> &n_mixture_candidates,
//              std::vector<unsigned int> &nparams,
//              std::vector<plFloat> &llk,
//              std::vector<plFloat> &bic,
//              std::vector<plJointDistribution> &model);

// void run_em(  const unsigned int n_data, const std::vector<unsigned int> &n_mixture_candidates,
//              std::vector<unsigned int> &nparams,
//              std::vector<plFloat> &llk,
//              std::vector<plFloat> &bic,
//              std::vector<plJointDistribution> &model);

// BOOST_AUTO_TEST_CASE( Test_EM_With_GaussianMixture_Desc ) {
//   const unsigned int ndata = 10000;
//   //  EM-BIC BASED LEARNING
//   // Number of components (kernels) candidates
//   std::vector<unsigned int> n_mixture_candidates;
//   n_mixture_candidates.push_back(1); 
//   n_mixture_candidates.push_back(2); 
//   n_mixture_candidates.push_back(3);
//   n_mixture_candidates.push_back(4);
//   // Output parameters for each number of components (kernels) candidate
//   std::vector<unsigned int> nparams; // The number of parameters
//   std::vector<plFloat> llk; // The log-likelihood
//   std::vector<plFloat> bic; // The bic score
//   std::vector<plJointDistribution> model; // the learnt model
//   //  Run an EM for each number of components (kernels) candidate
//   run_em( ndata, n_mixture_candidates, 
//          nparams, llk, bic, model);

 
//   // Get the model with the best BIC score
//   unsigned int best_candidate_index = std::max_element(bic.begin(), bic.end()) - bic.begin();

//   BOOST_CHECK_EQUAL(n_mixture_candidates[best_candidate_index], 2);

//   output_mixture(model[best_candidate_index]);

// }

// ///////////////////////////////
// std::vector<std::vector<double>>* generate_data(unsigned int ndata) {
//   auto data = new std::vector<std::vector<double>>(ndata, std::vector<double>(2, -1.0));
//   const unsigned int nc = 2;
//   const plSymbol C("C", plIntegerType(0, nc-1));
//   // X observed variable
//   const plSymbol X("X", plRealType(-100.0, 100.0));
//   // Actual (PC) table
//   const plProbValue pr[] = {0.3, 0.7};
//   const plProbTable PC(C, pr);
//   // Actual P(X | C) distributions
//   plDistributionTable PX(X, C);
//   const plFloat mean[] = {-10.0, 10.0};
//   const plFloat sd[] = {1.0, 3.0};
//   for(unsigned int i = 0; i < nc; ++i) {
//     PX.push( plNormal(X, mean[i], sd[i]), int(i));
//   }
//   // Constructing the model
//   const plJointDistribution mixture(C^X, PC*PX);
//   plValues val_CX(C^X);
//   for(unsigned int i = 0; i < ndata; ++i) {
//     mixture.draw(val_CX);
//     if (i!=MISSING_OBS) {
//       (*data)[i][1] = val_CX[1];
//     } 
//   }
//   return data;
// }

// std::vector< std::vector<bool>>* generate_data_def(unsigned int ndata) {
//   auto dataDef = new std::vector<std::vector<bool>>(ndata, std::vector<bool>(2, MISSING));
//   for(unsigned int i = 0; i < ndata; ++i) {
//     if (i!=MISSING_OBS) {
//       (*dataDef)[i][1] = NOT_MISSING;  
//     }
//   }

//   return dataDef;
// }


// plDataDescriptor* get_data_desc(unsigned int ndata, const plSymbol& C, const plSymbol& X) {
//   auto generated_data = generate_data(ndata);
//   auto data_def = generate_data_def(ndata);

//   // for ( int i = 0; i < ndata; ++i ) {
//   //   std::cout << i << ": " << (*data_def)[i][0] << "," << (*data_def)[i][1] << " - "
//   //             << (*generated_data)[i][0] << "," << (*generated_data)[i][1] << std::endl;
//   // }
  
//   return new plMatrixDataDescriptor<double>(C^X, generated_data, data_def);

// }


// void run_em(  const unsigned int n_data, unsigned int nc,
//              unsigned int &nparams, plFloat &llk,
//              plFloat &bic, plJointDistribution &model) {
  
//   const plSymbol C("C", plIntegerType(0, nc-1));
//   const plSymbol X("X", plRealType(-100.0, 100.0));
//   // EM Initial distribution on the class (kernel) variable: P(C)
//   const bool random_prob = true;
//   const plProbTable pc_init(C, random_prob);
//   // EM Initial Gaussians : P(X | C)
//   plDistributionTable px_init(X, C);
//   for(unsigned int i = 0; i < nc; ++i) {
//     px_init.push( plNormal(X, -10 + plRandomFloat(20.), 1.0), 
//                   int(i));
//   }
//   // P(C) is learnt as an histogram
//   plLearnHistogram LC(C);
//   // P(X | C) is learnt as a set of gaussians (a gaussian for each value of C)
//   plCndLearnObject <plLearn1dNormal> LX(X, C);
//   // Creating the EM learner instance
//   auto myCSVdata = dynamic_cast<plMatrixDataDescriptor<double>*>(get_data_desc(n_data, C, X));
//   std::vector <plLearnObject*> learn_objs(2); learn_objs[0] = &LC; learn_objs[1] = &LX;
//   plEMLearner myEM(pc_init*px_init, learn_objs);
//   // Run untill convergence
//   myEM.run(*myCSVdata, 0.0001);
//   // Fill the output parameters
//   nparams = myEM.get_n_parameters();
//   llk = myEM.get_last_computed_loglikelihood();
//   bic = llk - 0.5*nparams*std::log(myCSVdata->get_n_records());
//   model = myEM.get_joint_distribution();


// }

// void run_em( const unsigned int n_data, const std::vector<unsigned int> &n_mixture_candidates,
//              std::vector<unsigned int> &nparams,
//              std::vector<plFloat> &llk,
//              std::vector<plFloat> &bic,
//              std::vector<plJointDistribution> &model) {
//   nparams.resize( n_mixture_candidates.size() );
//   llk.resize( n_mixture_candidates.size() );
//   bic.resize( n_mixture_candidates.size() );
//   model.resize( n_mixture_candidates.size() );
//   for(unsigned int i = 0; i < n_mixture_candidates.size(); ++i) {
//     const unsigned int nc = n_mixture_candidates[i];
//     run_em( n_data, nc, nparams[i], llk[i], bic[i], model[i]);
//   }

  
// }

BOOST_AUTO_TEST_SUITE_END()  /// Test InfoTheo ends here
