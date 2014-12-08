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
////////////////////////////////////////////////////////////////////////////////////////

struct MultiSourceNodeCriterion: public NodeCriterion<double, std::less<double> > {
  // typedef std::binary_function<double,double,bool> Comp;
  MultiSourceNodeCriterion( const std::vector<std::vector<double>>& thresholdsByLevel,
                            const std::vector<std::vector<double>>& actualScores,
                            const double overral_thres = 1.0 )
      : levelThresholds(thresholdsByLevel), scores(actualScores), overall_threshold(overral_thres) {}

  virtual double nodeValue( const Graph& g, const vertex_t& vertex) const {
    assert(scores.size() > 0);
    return scores[0][vertex];
  }

  virtual double referenceValue( const Graph& g, const vertex_t& vertex) const {
    int level = g[vertex].level;
    return levelThresholds[0][level];
  }  

  virtual bool isValid( const Graph& g, const vertex_t& vertex) const {
    
    int level = g[vertex].level;    
    double overall_score = 0.0;
    printf("vertex: %d, level: %d, ", vertex, level);
    for ( size_t c = 0; c < scores.size(); ++c ) {
      double level_val = levelThresholds[c][level];
      double score_val = scores[c][vertex];
      overall_score += func(score_val,level_val);
    }

    printf("overall_s: %f, overall_t: %f, score <= thres: %d\n", overall_score, overall_threshold,
           (overall_score >= overall_threshold));

    return overall_score >= overall_threshold;
  }

  
  void setOverallThreshold( const double thres ) {
    this->overall_threshold = thres;
  }
  
 private:
  const std::vector<std::vector<double>>& levelThresholds; // scores[i][j] = criteria-i-level-j
  const std::vector<std::vector<double>>& scores; // scores[i][j] = criteria-i-level-j
  double overall_threshold;
};

} // namespace samogwasends here. samogwas

/****************************************************************************************/
#endif // SAMOGWAS_SCORE_NODE_CRITERIA_HPP
