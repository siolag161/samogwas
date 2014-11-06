#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
#   define BOOST_TEST_MODULE
#endif  
#include <boost/test/unit_test.hpp>
#include <pl.h>
#include <fstream>
#include <iostream>



class Data 
{ 
};

static const std::string FILE_NAME = "data/em_gaussian_data.csv";
static const int CARDI = 2;
static const plProbValue PR[] = {0.3, 0.7}; // actual
static const int NOT_MISSING = 1;
static const int MISSING = 0;

typedef std::vector<double> Row;
typedef std::vector<Row> Matrix;
typedef plMatrixDataDescriptor<double> MatrixDesc;
/////////////////////////////////////////////////////////////////////

BOOST_FIXTURE_TEST_SUITE( Test_Missing_Data, Data )

Matrix generate_data(unsigned n);
void run_em( Matrix& data,
             unsigned int &nparams, plFloat &llk,
             plFloat &bic, plJointDistribution &model );


// MatrixDesc create_matrix_desc( Matrix& data );


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


///////////////////////////////////////////////////////////////////
Matrix generate_data(unsigned nrows) {
  Matrix data(nrows, std::vector<double>(2, 0.0));
  
  const plSymbol C("C", plIntegerType(0, CARDI-1));
  const plSymbol X("X", plRealType(-100.0, 100.0));
  const plProbTable PC(C, PR);

  // Actual P(X | C) distributions
  plDistributionTable PX(X, C);
  const plFloat mean[] = {-10.0, 10.0};
  const plFloat sd[] = {1.0, 3.0};
  for(unsigned int i = 0; i < CARDI; ++i) {
    PX.push( plNormal(X, mean[i], sd[i]), int(i));
  }
  
  // Constructing the model
  //
  plVariablesConjunction varConj = C^X;
  const plJointDistribution mixture(varConj, PC*PX);

  plValues val_CX(C^X);
  for(unsigned int i = 0; i < nrows; ++i) {
    mixture.draw(val_CX);
    data[i][0] = -1;
    data[i][1] =  val_CX[X];
  }

  return data;
}


void run_em( Matrix& data,
             unsigned int &nparams, plFloat &llk,
             plFloat &bic, plJointDistribution &model )
{
  const plSymbol C("C", plIntegerType(0, CARDI-1));
  const plSymbol X("X", plRealType(-100.0, 100.0));

  // EM Initial distribution on the class (kernel) variable: P(C)
  const bool random_prob = true;
  const plProbTable pc_init(C, random_prob);

  // EM Initial Gaussians : P(X | C)
  plDistributionTable px_init(X, C);
  for(unsigned int i = 0; i < CARDI; ++i) {
    px_init.push( plNormal(X, -10 + plRandomFloat(20.), 1.0), 
                  int(i));
  }

  // P(C) is learnt as an histogram
  plLearnHistogram LC(C);
  // P(X | C) is learnt as a set of gaussians (a gaussian for each value of C)
  plCndLearnObject <plLearn1dNormal> LX(X, C);

  plVariablesConjunction varConj = C^X;
  std::vector< std::vector<bool> >* def =
      new std::vector< std::vector<bool> >( data.size(), std::vector<bool>(2, NOT_MISSING));
  for (int i = 0; i < data.size();++i) {  
    (*def)[i][0] = MISSING;
  }

  MatrixDesc dataDesc(varConj, &data, def); // error

  std::vector <plLearnObject*> learn_objs(2); learn_objs[0] = &LC; learn_objs[1] = &LX;
  plEMLearner myEM(pc_init*px_init, learn_objs);
  // Run untill convergence
  myEM.run( dataDesc, 0.0001);

  // Fill the output parameters
  nparams = myEM.get_n_parameters();
  llk = myEM.get_last_computed_loglikelihood();
  bic = llk - 0.5*nparams*std::log( dataDesc.get_n_records());
  model = myEM.get_joint_distribution();
}



BOOST_AUTO_TEST_SUITE_END() 
