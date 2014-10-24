/****************************************************************************************
 * File: score_node_criteria.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 21/10/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_SCORE_NODE_CRITERIA_HPP
#define SAMOGWAS_SCORE_NODE_CRITERIA_HPP

#include <functional>   // std::less
#include <statistics/association_test.hpp>
#include "node_criteria.hpp"
namespace samogwas
{

struct LevelScoreNodeCriterion: public NodeCriterion<double, std::less<double> > {
  LevelScoreNodeCriterion( std::vector<double>& thresholdsByLevel,
                           const std::vector<double>& actualScores )
      : levelThresholds(thresholdsByLevel), scores(actualScores) {}

  virtual double nodeValue( const Graph& g, const vertex_t& vertex) const {
    return scores[vertex];
  }

  virtual double referenceValue( const Graph& g, const vertex_t& vertex) const {
    int level = g[vertex].level;
    return levelThresholds[level];
  }
  
  
 private:
  std::vector<double>& levelThresholds;
  const std::vector<double>& scores;
};

////////////////////////////////////////////////////////////////////////////////////
// template<Test
struct StatTestNodeCriterion: public NodeCriterion<double, std::less<double> > {
  StatTestNodeCriterion( std::vector<double>& thresholdsByLevel,
                         std::shared_ptr<stats::GWASAssociationTest> test )
      : levelThresholds(thresholdsByLevel), statTest(test) {}

  virtual double nodeValue( const Graph& g, const vertex_t& vertex) const {
    return statTest->execute((size_t)vertex, g[vertex].variable.cardinality());
  }

  virtual double referenceValue( const Graph& g, const vertex_t& vertex) const {
    int level = g[vertex].level;
    return levelThresholds[level];
  }
  
  
 private:
  std::vector<double>& levelThresholds;  
  std::shared_ptr<stats::GWASAssociationTest> statTest;
};

} // namespace samogwasends here. samogwas

/****************************************************************************************/
#endif // SAMOGWAS_SCORE_NODE_CRITERIA_HPP
