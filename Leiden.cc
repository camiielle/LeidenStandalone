#include "TICLGraph.h"
#include "Leiden.h"
#include <cmath>
#include <queue>
#include <cassert>
#include <cmath>
#include <random>
#include <string>
#include <utility>
#include <iostream>

namespace ticl {

  void leidenAlgorithm(TICLGraph &graph, Partition &partition, std::vector<Flat> &flatFinalPartition, int gamma, double theta) {
    moveNodesFast(graph, partition);
    //moveNodesFast(partition, gamma_);
    std::cout << "MOVED NODES FAST\n";

    if (!(isAlgorithmDone(graph, partition))) {
      std::cout << __LINE__ << std::endl;
      Partition refinedPartition = Partition{std::vector<Community>{}};
      std::cout << __LINE__ << std::endl;
      assert(refinedPartition.getCommunities().empty());
      std::cout << __LINE__ << std::endl;
      //PROBLEM SEEMS TO BE HERE
      refinePartition(graph, partition, refinedPartition, gamma, theta);
      std::cout << __LINE__ << std::endl;
      //create aggregate graph based on refined partition P_ref
      aggregateGraph(graph, refinedPartition);
      std::cout << __LINE__ << std::endl;
      //but maintain partition P
      auto &communities = partition.getCommunities();
      std::cout << __LINE__ << std::endl;
      std::vector<Community> aggregatedCommunities{};
      std::cout << __LINE__ << std::endl;
      for (auto const &community : communities) {
        Community aggregatedCommunity{{}, degree(community) + 1};
        std::cout << __LINE__ << std::endl;
        for (auto const &aggregateNode : graph.getNodes()) {
          if (isCommunityContained(std::get<Community>(aggregateNode), community)) {
            aggregatedCommunity.getNodes().push_back(aggregateNode);
          }
        }
        std::cout << __LINE__ << std::endl;
        aggregatedCommunities.push_back(aggregatedCommunity);
        std::cout << __LINE__ << std::endl;
      }
      communities = aggregatedCommunities;
      std::cout << __LINE__ << std::endl;
      leidenAlgorithm(graph, partition, flatFinalPartition, gamma, theta);
    }

    else {
      std::cout << __LINE__ << std::endl;
      partition.flattenPartition(flatFinalPartition);
      std::cout << __LINE__ << std::endl;
    }
  }

  bool isAlgorithmDone(TICLGraph const &graph, Partition const &partition) {
    return partition.getCommunities().size() == graph.getNodes().size();
  }

  //quality function, Constant Potts Model
  //tested on godbolt
  auto CPM(Partition const &partition, int gamma) {
    int CPMResult{};
    for (auto const &community : partition.getCommunities()) {
      int n{communitySize(community, 0)};
      CPMResult += (numberOfEdges(community, community) - gamma * n * (n - 1) / 2);
    }
    return CPMResult;
  }

  //interpreting E(C,C) as non null even if C has a single node, if degree(node)>0
  /*auto CPM_contribution_from_new_community(Node const &node, int gamma) {
  Community newCommunity{std::vector<Node>{node}, degree(node) + 1};
  int n = {communitySize(newCommunity, 0)};
  int result{numberOfEdges(newCommunity, newCommunity) - gamma * n * (n - 1) / 2};
  return result;
}*/

  /*auto CPM_after_move(Partition const &partition,
                    int gamma,
                    Community const &communityFrom,
                    Community const &communityTo,
                    Node const &node) {
  int CPMResult{};
  auto const &communities = partition.getCommunities();
  for (auto const &community : communities) {
    if (community == communityFrom && communityFrom.getNodes().size() > 1) {
      std::vector<Node> vectorWithoutNode{};
      std::copy_if(communityFrom.getNodes().begin(),
                   communityFrom.getNodes().end(),
                   std::back_inserter(vectorWithoutNode),
                   [&](Node const &n) { return !(n == node); });
      Community communityWithoutNode{vectorWithoutNode, communityFrom.getDegree()};
      CPMResult += (numberOfEdges(communityWithoutNode, communityWithoutNode) -
                    gamma * (communitySize(communityWithoutNode, 0) * communitySize(communityWithoutNode, 0) / 2));
    } else if (community == communityTo) {
      Community communityWithNewNode{community};
      communityWithNewNode.getNodes().push_back(node);
      CPMResult += (numberOfEdges(communityWithNewNode, communityWithNewNode) -
                    gamma * (communitySize(communityWithNewNode, 0) * communitySize(communityWithNewNode, 0) / 2));

    } else {
      CPMResult += (numberOfEdges(community, community) -
                    gamma * (communitySize(community, 0) * communitySize(community, 0) / 2));
      //PROBLEM CAUSED AFTER THIS LINE
    }
  }
  return CPMResult;
}*/

  auto delta_CPM_after_move(int gamma, Community const &communityFrom, Community const &communityTo, Node const &node) {
    std::vector<Node> vectorWithoutNode{};
    std::copy_if(communityFrom.getNodes().begin(),
                 communityFrom.getNodes().end(),
                 std::back_inserter(vectorWithoutNode),
                 [&](Node const &n) { return !(n == node); });
    Community communityWithoutNode{vectorWithoutNode, communityFrom.getDegree()};
    int first_term{numberOfEdges(communityWithoutNode, communityWithoutNode) -
                   numberOfEdges(communityFrom, communityFrom) + communitySize(communityFrom, 0) - 1};
    Community communityWithNewNode{communityTo};
    communityWithNewNode.getNodes().push_back(node);
    int second_term{numberOfEdges(communityWithNewNode, communityWithNewNode) -
                    numberOfEdges(communityTo, communityTo) - communitySize(communityTo, 0)};

    return first_term + second_term;
  }

  double delta_modularity_after_move(int totalEdges,
                                     Community const &communityFrom,
                                     Community const &communityTo,
                                     Node const &node) {
    std::vector<Node> vectorWithoutNode{};
    std::copy_if(communityFrom.getNodes().begin(),
                 communityFrom.getNodes().end(),
                 std::back_inserter(vectorWithoutNode),
                 [&](Node const &n) { return !(n == node); });
    Community communityWithoutNode{vectorWithoutNode, communityFrom.getDegree()};
    auto first_term =
        numberOfEdges(communityWithoutNode, communityWithoutNode) - numberOfEdges(communityFrom, communityFrom);

    Community communityWithNewNode{communityTo};
    communityWithNewNode.getNodes().push_back(node);
    auto second_term =
        numberOfEdges(communityWithNewNode, communityWithNewNode) - numberOfEdges(communityTo, communityTo);

    auto sumDegreeFrom = kappa(communityFrom);
    auto sumDegreeTo = kappa(communityTo);
    double third_term = kappa(node) * (sumDegreeFrom - sumDegreeTo) / (2. * totalEdges);

    return (first_term + second_term + third_term) / totalEdges;
  }

  auto delta_CPM_from_empty(Node const &node, int gamma, Community const &communityFrom) {
    std::vector<Node> vectorWithoutNode{};
    std::copy_if(communityFrom.getNodes().begin(),
                 communityFrom.getNodes().end(),
                 std::back_inserter(vectorWithoutNode),
                 [&](Node const &n) { return !(n == node); });
    Community communityWithoutNode{vectorWithoutNode, communityFrom.getDegree()};
    int first_term{numberOfEdges(communityWithoutNode, communityWithoutNode) -
                   numberOfEdges(communityFrom, communityFrom) + communitySize(communityFrom, 0) - 1};

    Community newCommunity{std::vector<Node>{node}, degree(node) + 1};
    return first_term + numberOfEdges(newCommunity, newCommunity);
  }

  double delta_modularity_from_empty(int totalEdges, Community const &communityFrom, Node const &node) {
    std::vector<Node> vectorWithoutNode{};
    std::copy_if(communityFrom.getNodes().begin(),
                 communityFrom.getNodes().end(),
                 std::back_inserter(vectorWithoutNode),
                 [&](Node const &n) { return !(n == node); });
    Community communityWithoutNode{vectorWithoutNode, communityFrom.getDegree()};
    int first_term{numberOfEdges(communityWithoutNode, communityWithoutNode) -
                   numberOfEdges(communityFrom, communityFrom)};

    Community newCommunity{{node}, degree(node) + 1};
    int second_term{numberOfEdges(newCommunity, newCommunity)};

    int sumDegreeFrom{kappa(communityFrom)};
    double third_term = kappa(node) * sumDegreeFrom / (2. * totalEdges);

    return ((first_term + second_term + third_term) / totalEdges);
  }

  void moveNode(Community &communityFrom, Community &communityTo, Node const &node) {
    auto it = std::find(communityFrom.getNodes().begin(), communityFrom.getNodes().end(), node);
    assert(it != communityFrom.getNodes().end());
    assert(&(*it) != &node);
    communityFrom.getNodes().erase(it);
    communityTo.getNodes().push_back(node);
  }

  /*void queueCommunity(Community &community, std::queue<Node> &queue) {
    //elements are added to the queue in random order
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(community.getNodes().begin(), community.getNodes().end(), g);

    for (auto const &node : community.getNodes()) {
      queue.push(node);
    }
  }*/

  /*Partition &removeEmptyCommunities(Partition &partition) {
    auto &communities = partition.getCommunities();
    communities.erase(std::remove_if(communities.begin(), communities.end(), [](Community const &community) {
      return community.getNodes().size() == 0;
    }));

    auto const &communitiesAfterRemoval = partition.getCommunities();
    for (auto const &communityAfterRemoval : communitiesAfterRemoval) {
      assert(communityAfterRemoval.getNodes().size() != 0);
    }

    return partition;
  }*/

  /*int bestCommunityIndex(std::vector<Community> const &communities,
                         int &bestDeltaCPM,
                         Community const &currentCommunity,
                         Node const &currentNode,
                         int gamma) {
    //variation of quality function if I move the node to a different community
    assert(bestDeltaCPM == 0);
    int indexBestCommunity{};
    for (unsigned int i = 0; i < communities.size(); ++i) {
      //PROBLEM IS IN THE LINE BELOW!!!!!
      auto deltaCPM = delta_CPM_after_move(gamma, currentCommunity, communities[i], currentNode);
      if (deltaCPM > bestDeltaCPM) {
        bestDeltaCPM = deltaCPM;
        indexBestCommunity = i;
      }
    }
    return indexBestCommunity;
  }*/

  int bestCommunityIndexModularity(std::vector<Community> const &communities,
                                   double &bestDeltaModularity,
                                   Community const &currentCommunity,
                                   Node const &currentNode,
                                   int totalEdges) {
    //variation of quality function if I move the node to a different community
    assert(bestDeltaModularity == 0.);
    int indexBestCommunity{};
    for (unsigned int i = 0; i < communities.size(); ++i) {
      auto deltaModularity = delta_modularity_after_move(totalEdges, currentCommunity, communities[i], currentNode);
      if ((deltaModularity - bestDeltaModularity) > 0.0001) {
        bestDeltaModularity = deltaModularity;
        indexBestCommunity = i;
      }
    }
    return indexBestCommunity;
  }

  /*template <class ADAPTER>
  struct hack : private ADAPTER {
    static auto &get(ADAPTER &a) { return a.*(&hack::c); }
  };

  template <class ADAPTER>
  const auto &get_container(ADAPTER &a) {
    return hack<ADAPTER>::get(a);
  }*/

  /*Partition &moveNodesFast(Partition &partition, int gamma) {
  //all nodes are added to queue in random order
  std::vector<Node> allNodes{};
  std::for_each(
      partition.getCommunities().begin(), partition.getCommunities().end(), [&allNodes](Community const &community) {
        for (auto const &node : community.getNodes())
          allNodes.push_back(node);
      });
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(allNodes.begin(), allNodes.end(), g);
  std::queue<Node> queue{};
  //nodes are added to the queue in random order
  for (auto const &node : allNodes) {
    queue.push(node);
  }

  //std::cout << "PARTITION SIZE: " << partition.getCommunities().size() << std::endl;
  //std::cout << "TOTAL N OF NODES: " << allNodes.size() << std::endl;

  while (!(queue.empty())) {
    Node currentNode{queue.front()};
    queue.pop();
    auto &queue_members = get_container(queue);
    //std::cout << "NODES IN QUEUE: " << queue_members.size() << std::endl;
    auto &currentCommunity = partition.getCommunities()[partition.findCommunityIndex(currentNode)];
    auto &communities = partition.getCommunities();
    //std::cout << "NODE ID " << getId(currentNode) << "COMMUNITY INDEX " << partition.findCommunityIndex(currentNode)
    // << std::endl;
    //std::cout << "N OF COMMUNITIES IN PARTITION: " << communities.size() << std::endl;

    int bestDeltaCPM{0};
    int indexBestCommunity{bestCommunityIndex(communities, bestDeltaCPM, currentCommunity, currentNode, gamma)};

    //variation of quality function if I move the node to an empty community
    int deltaCPMFromEmpty{0};
    if (currentCommunity.getNodes().size() != 1) {
      deltaCPMFromEmpty = delta_CPM_from_empty(currentNode, gamma, currentCommunity);
    }
    //std::cout << "BEST CPM " << bestDeltaCPM << "CPM FROM NEW " << deltaCPMFromEmpty << std::endl;

    if (deltaCPMFromEmpty > bestDeltaCPM) {
      bestDeltaCPM = -1;
      Community newCommunity{std::vector<Node>{}, degree(currentNode) + 1};
      communities.push_back(newCommunity);
      moveNode(currentCommunity, newCommunity, currentNode);
      // making sure all nbrs of currentNode who are not in bestCommunity will be visited
      std::for_each(communities.begin(), communities.end() - 1, [&](auto const &community) {
        std::for_each(community.getNodes().begin(), community.getNodes().end(), [&](auto const &node) {
          if (areNeighbours(currentNode, node) &&
              std::find(queue_members.begin(), queue_members.end(), node) == queue_members.end()) {
            queue.push(node);
          }
        });
      });
    }

    if (bestDeltaCPM > 0) {
      moveNode(currentCommunity, communities[indexBestCommunity], currentNode);
      // making sure all nbrs of currentNode who are not in bestCommunity will be visited
      for (unsigned int i = 0; i < communities.size(); ++i) {
        if (i != static_cast<unsigned int>(indexBestCommunity)) {
          std::for_each(communities[i].getNodes().begin(), communities[i].getNodes().end(), [&](auto const &node) {
            if (areNeighbours(currentNode, node) &&
                std::find(queue_members.begin(), queue_members.end(), node) == queue_members.end()) {
              queue.push(node);
            }
          });
        }
      }
    }

    if (currentCommunity.getNodes().size() == 0) {
      communities.erase(std::remove(communities.begin(), communities.end(), currentCommunity));
    }
  }

  return partition;
}*/

  Partition &moveNodesFast(TICLGraph const &graph, Partition &partition) {
    auto const nEdges = totalEdges(graph);
    auto const nNodes = graph.getNodes().size();

    //all nodes are added to queue in random order
    std::vector<Node> queue{};
    queue.reserve(nNodes * 4);
    auto &communities = partition.getCommunities();
    std::for_each(communities.begin(), communities.end(), [&](Community const &community) {
      auto &nodes = community.getNodes();
      queue.insert(queue.end(), nodes.begin(), nodes.end());
    });
    std::random_device rd;
    std::default_random_engine g(rd());
    std::shuffle(queue.begin(), queue.end(), g);

    std::size_t front_index{0};
    while (front_index != queue.size()) {
      auto currentNode = queue[front_index];
      ++front_index;
      auto currentCommunityIndex = partition.findCommunityIndex(currentNode);
      auto &communities = partition.getCommunities();
      auto &currentCommunity = communities[currentCommunityIndex];

      double bestDeltaModularity{0};
      int indexBestCommunity{
          bestCommunityIndexModularity(communities, bestDeltaModularity, currentCommunity, currentNode, nEdges)};

      //variation of quality function if I move the node to an empty community
      int deltaModularityFromEmpty{0};
      if (currentCommunity.getNodes().size() != 1) {
        deltaModularityFromEmpty = delta_modularity_from_empty(nEdges, currentCommunity, currentNode);
      }

      if ((deltaModularityFromEmpty - bestDeltaModularity) > 0.0001) {
        Community newCommunity{{}, degree(currentNode) + 1};
        communities.push_back(newCommunity);
        moveNode(currentCommunity, newCommunity, currentNode);
        // making sure all nbrs of currentNode who are not in bestCommunity will be visited
        std::for_each(communities.begin(), communities.end() - 1, [&](auto const &community) {
          std::for_each(community.getNodes().begin(), community.getNodes().end(), [&](auto const &node) {
            if (areNeighbours(currentNode, node) &&
                std::find(queue.begin() + front_index, queue.end(), node) == queue.end()) {
              queue.push_back(node);
            }
          });
        });
      } else if (bestDeltaModularity > 0.0001) {
        moveNode(currentCommunity, communities[indexBestCommunity], currentNode);
        // making sure all nbrs of currentNode who are not in bestCommunity will be visited
        for (unsigned int i = 0; i < communities.size(); ++i) {
          if (i != static_cast<unsigned int>(indexBestCommunity)) {
            std::for_each(communities[i].getNodes().begin(), communities[i].getNodes().end(), [&](auto const &node) {
              if (areNeighbours(currentNode, node) &&
                  std::find(queue.begin() + front_index, queue.end(), node) == queue.end()) {
                queue.push_back(node);
              }
            });
          }
        }
      }

      if (currentCommunity.getNodes().empty()) {
        auto it = communities.begin() + currentCommunityIndex;
        assert(it != communities.end());
        communities.erase(it);
      }
    }

    return partition;
  }

  //fills an empty partition with a singleton partition
  Partition &singletonPartition(TICLGraph const &graph, Partition &singlePartition) {
    assert(singlePartition.getCommunities().empty());
    auto const &nodes = graph.getNodes();
    auto &communities = singlePartition.getCommunities();
    for (auto const &node : nodes) {
      Community singletonCommunity{std::vector<Node>{node}, degree(node) + 1};
      communities.push_back(singletonCommunity);
    }
    assert(!(singlePartition.getCommunities().empty()));
    assert(singlePartition.getCommunities().size() == nodes.size());
    for (auto const &community : singlePartition.getCommunities()) {
      assert(community.getNodes().size() == 1);
    }
    return singlePartition;
  }

  bool isNodeWellConnected(Node const &node, Community const &subset, int gamma) {
    auto nodes = subset.getNodes();
    auto it = std::find(nodes.begin(), nodes.end(), node);
    if (it != nodes.end()) {
      nodes.erase(it);
    }
    Community singletonCommunity{{node}, degree(node) + 1};
    Community subsetWithoutNode{nodes, degree(subset)};
    auto edges{numberOfEdges(singletonCommunity, subsetWithoutNode)};
    assert(edges >= 0);
    auto nodeSize{communitySize(singletonCommunity, 0)};
    auto subsetSize{communitySize(subset, 0)};
    return (edges >= (gamma * nodeSize * (subsetSize - nodeSize)));
  }

  bool isCommunityWellConnected(Community const &community, Community const &subset, int gamma) {
    Community subsetMinuscommunity{{}, degree(subset)};
    for (auto const &node : subset.getNodes()) {
      if (std::find(community.getNodes().begin(), community.getNodes().end(), node) == community.getNodes().end()) {
        subsetMinuscommunity.getNodes().push_back(node);
      }
    }
    auto edges{numberOfEdges(community, subsetMinuscommunity)};
    assert(edges >= 0);
    auto comSize{communitySize(community, 0)};
    auto subsetSize{communitySize(subset, 0)};
    return (edges >= (gamma * comSize * (subsetSize - comSize)));
  }

  int extractRandomCommunityIndex(std::vector<Community> const &communities,
                                  Partition const &partition,
                                  Node const &node,
                                  Community const &nodeCommunity,
                                  Community const &subset,
                                  int gamma,
                                  double theta) {
    std::vector<double> deltaCPMs{};

    //calculating delta_H for all communities
    for (auto const &community : communities) {
      if (isCommunityWellConnected(community, subset, gamma)) {
        deltaCPMs.push_back((delta_CPM_after_move(gamma, nodeCommunity, community, node)));
      } else {
        // communities not well connected are not considered
        deltaCPMs.push_back(-1.);
      }
    }

    //creating the discrete probability function
    std::vector<double> distribution{};
    for (auto const &deltaCPM : deltaCPMs) {
      if (deltaCPM < 0.) {
        distribution.push_back(0.);
      } else {
        assert(theta > 0.);
        distribution.push_back(std::exp(deltaCPM / theta));
      }
    }

    //extracting a random community
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> d(distribution.begin(), distribution.end());
    //extracts a random index
    int resultIndex = d(gen);

    return resultIndex;
  }

  Partition &mergeNodesSubset(Partition &partition, Community const &subset, int gamma, double theta) {
    auto &communities = partition.getCommunities();

    //consider only nodes that are well connected within subset S
    for (auto const &node : subset.getNodes()) {
      if (isNodeWellConnected(node, subset, gamma)) {
        int index = partition.findCommunityIndex(node);
        auto &nodeCommunity = communities[index];
        std::cout << "COMM (non recursive) SIZE " << nodeCommunity.getNodes().size();
        std::cout << " COMM SIZE " << communitySize(nodeCommunity, 0);
        assert(communitySize(nodeCommunity, 0) != 0);
        //consider only nodes that have not yet been merged
        if (communitySize(nodeCommunity, 0) == 1) {
          int communityToIndex{
              extractRandomCommunityIndex(communities, partition, node, nodeCommunity, subset, gamma, theta)};
          auto &communityTo = communities[communityToIndex];
          moveNode(nodeCommunity, communityTo, node);
          //removing empty communities
          if (nodeCommunity.getNodes().empty()) {
            auto it = communities.begin() + index;
            assert(it != communities.end());
            communities.erase(it);
          }
        }
      }
    }
    return partition;
  }

  Partition &refinePartition(
      TICLGraph const &graph, Partition &partition, Partition &singlePartition, int gamma, double theta) {
    //fills an empty partition with a singleton partition
    auto &refinedPartition = singletonPartition(graph, singlePartition);
    auto const &communities = partition.getCommunities();
    for (auto const &community : communities) {
      mergeNodesSubset(refinedPartition, community, gamma, theta);
    }
    return refinedPartition;
  }

  //communities become nodes in aggregate graph
  TICLGraph &aggregateGraph(TICLGraph &graph, Partition const &partition) {
    std::vector<Community> const &communities{partition.getCommunities()};
    std::vector<Node> aggregatedNodes{};
    aggregatedNodes.reserve(communities.size());
    std::for_each(communities.begin(), communities.end(), [&aggregatedNodes](auto const &community) {
      aggregatedNodes.push_back(Node{community});
    });

    assert(aggregatedNodes.size() == communities.size());
    graph.setNodes(aggregatedNodes);
    return graph;
  }
}  // namespace ticl