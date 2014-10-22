/****************************************************************************************
 * File: score_node_criteria.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 21/10/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_SCORE_NODE_CRITERIA_HPP
#define SAMOGWAS_SCORE_NODE_CRITERIA_HPP

#include "node_criteria.hpp"
namespace samogwas
{

struct LevelScoreNodeCriterion: public NodeCriterion {
  LevelScoreNodeCriterion( std::vector<double>& thresholdsByLevel,
                           const std::vector<double>& actualScores)
      : levelThresholds(thresholdsByLevel), scores(actualScores) {}

  virtual bool isValid( const Graph& graph, const vertex_t& vertex, bool less = true ) const {
    int level = graph[vertex].level;
    double score = scores[vertex];
    double threshold = levelThresholds[level];
    return less ? (score < threshold) : score > threshold;
  }
  
 private:
  std::vector<double>& levelThresholds;
  const std::vector<double>& scores;
};

} // namespace samogwasends here. samogwas

/****************************************************************************************/
#endif // SAMOGWAS_SCORE_NODE_CRITERIA_HPP
