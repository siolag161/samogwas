/****************************************************************************************
* File: gwas_strategy.hpp
* Description: This module provides common interfaces for GWAS strategies.
*
* @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
* @date: 01/08/2014

***************************************************************************************/

#ifndef SAMOGWAS_GWAS_STRATEGY_HPP
#define SAMOGWAS_GWAS_STRATEGY_HPP

#include <vector>

#include "fltm/core_fltm.hpp"
#include "statistics/association_test.hpp"

namespace samogwas
{
    class GWAS_Strategy {
    public:
        typedef std::vector<int> Vector;
        typedef std::vector<Vector> Matrix;
        virtual void execute(FLTM_Result& result, Matrix& genotype, Vector& phenotype, stats::StatTest* statTest) = 0;
    };

} // namespace samogwas ends here.

/**************************************** IMPLEMENTATION BELOW THIS POINT ****************************************/


/****************************************************************************************/
#endif // SAMOGWAS_GWAS_STRATEGY_HPP
