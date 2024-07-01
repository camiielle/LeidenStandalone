#ifndef Leiden_H__
#define Leiden_H__

#include <memory>
#include <array>
#include "TICLGraph.h"

namespace ticl {

  void leidenAlgorithm(TICLGraph &graph, Partition &partition, std::vector<Flat> &flatFinalPartition, int gamma, double theta);

  bool isAlgorithmDone(TICLGraph const &graph, Partition const &partition, bool hasNodeBeenMoved);

  Partition &removeEmptyCommunities(Partition &partition);

  Partition &refinePartition(
      TICLGraph const &graph, Partition const&partition, Partition &singlePartition, int gamma, int nEdges, double theta);

  Partition &moveNodesFast(TICLGraph const &graph, Partition &partition, int nEdges, bool& hasNodeBeenMoved);

  Partition &singletonPartition(TICLGraph const &graph, Partition &singlePartition);

  Partition &mergeNodesSubset(Partition &partition, Community const &subset, int gamma, int nEdges, double theta);

  TICLGraph &aggregateGraph(TICLGraph &graph, Partition const& partition, Partition const &refinedPartition);
}  // namespace ticl

#endif

