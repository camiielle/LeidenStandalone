#include <iostream>
#include <vector>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <iostream>
#include <vector>
#include <fstream>
#include "Trackster.h"
#include "TICLLayerTiles.h"
#include "TICLGraph.h"
#include "Leiden.h"

using namespace ticl;
// Function to read the file paths from a text file
std::vector<std::string> readFilePaths(const std::string& txtFileName) {
  std::vector<std::string> fileNames;
  std::ifstream infile(txtFileName);
  std::string line;
  while (std::getline(infile, line)) {
    if (!line.empty() && line.front() != '#') {
      fileNames.push_back(line);
    }
  }
  return fileNames;
}

void TICLGraphProducer(std::vector<Trackster>& trackstersclue3d, TICLGraph& graph) {
  TICLLayerTile tracksterTilePos;
  TICLLayerTile tracksterTileNeg;

  for (size_t id_t = 0; id_t < trackstersclue3d.size(); ++id_t) {
    auto t = trackstersclue3d[id_t];
    if (t.barycenter().eta() > 0.) {
      tracksterTilePos.fill(t.barycenter().eta(), t.barycenter().phi(), id_t);
    } else if (t.barycenter().eta() < 0.) {
      tracksterTileNeg.fill(t.barycenter().eta(), t.barycenter().phi(), id_t);
    }
  }

  std::vector<Elementary> allElemNodes;

  for (size_t id_t = 0; id_t < trackstersclue3d.size(); ++id_t) {
    auto t = trackstersclue3d[id_t];

    Elementary tNode(id_t);

    auto bary = t.barycenter();
    double del = 0.1;

    double eta_min = std::max(abs(bary.eta()) - del, (double)TileConstants::minEta);
    double eta_max = std::min(abs(bary.eta()) + del, (double)TileConstants::maxEta);

    if (bary.eta() > 0.) {
      std::array<int, 4> search_box =
          tracksterTilePos.searchBoxEtaPhi(eta_min, eta_max, bary.phi() - del, bary.phi() + del);
      if (search_box[2] > search_box[3]) {
        search_box[3] += TileConstants::nPhiBins;
      }

      for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
        for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
          auto& neighbours = tracksterTilePos[tracksterTilePos.globalBin(eta_i, (phi_i % TileConstants::nPhiBins))];
          for (auto n : neighbours) {
            if (trackstersclue3d[n].barycenter().z() < bary.z()) {
              tNode.addInnerNeighbour(n);
              tNode.addNeighbour(n);
            } else if (trackstersclue3d[n].barycenter().z() > bary.z()) {
              tNode.addOuterNeighbour(n);
              tNode.addNeighbour(n);
            }
          }
        }
      }
    }

    else if (bary.eta() < 0.) {
      std::array<int, 4> search_box =
          tracksterTileNeg.searchBoxEtaPhi(eta_min, eta_max, bary.phi() - del, bary.phi() + del);
      if (search_box[2] > search_box[3]) {
        search_box[3] += TileConstants::nPhiBins;
      }

      for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
        for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
          auto& neighbours = tracksterTileNeg[tracksterTileNeg.globalBin(eta_i, (phi_i % TileConstants::nPhiBins))];
          for (auto n : neighbours) {
            if (abs(trackstersclue3d[n].barycenter().z()) < abs(bary.z())) {
              tNode.addInnerNeighbour(n);
              tNode.addNeighbour(n);
            } else if (abs(trackstersclue3d[n].barycenter().z()) > abs(bary.z())) {
              tNode.addOuterNeighbour(n);
              tNode.addNeighbour(n);
            }
          }
        }
      }
    }
    allElemNodes.push_back(tNode);
  }
  std::vector<Node> allNodes{};
  allNodes.reserve(allElemNodes.size());
  for (auto const& e : allElemNodes) {
    Node node{e};
    allNodes.push_back(node);
  }
  graph.setNodes(allNodes);
}

void processEvent(TTree* tree, int ev) {
  std::cout << "Began processing event: " << ev << std::endl;
  // Access specific leaves of the "tracksters" TTree
  std::vector<float>*barycenter_x = nullptr, *barycenter_y = nullptr, *barycenter_z = nullptr, *raw_energy = nullptr,
  *barycenter_eta = nullptr, *barycenter_phi = nullptr;

  tree->SetBranchAddress("barycenter_x", &barycenter_x);
  tree->SetBranchAddress("barycenter_y", &barycenter_y);
  tree->SetBranchAddress("barycenter_z", &barycenter_z);
  tree->SetBranchAddress("barycenter_eta", &barycenter_eta);
  tree->SetBranchAddress("barycenter_phi", &barycenter_phi);
  tree->SetBranchAddress("raw_energy", &raw_energy);
  std::vector<Trackster> tracksters;
  tree->GetEntry(ev);
  for (size_t j = 0; j < barycenter_x->size(); ++j) {
    Trackster t;
    t.setBarycenter(
        barycenter_eta->at(j), barycenter_phi->at(j), barycenter_x->at(j), barycenter_y->at(j), barycenter_z->at(j));
    t.setRawEnergy(raw_energy->at(j));
    tracksters.push_back(t);
    //t.Print();
  }
  if (tracksters.size() != 0) {
    std::cout << "NUM OF TRACKSTERS: " << tracksters.size() << std::endl;
    TICLGraph graph;
    TICLGraphProducer(tracksters, graph);
    std::cout << "Running algo on event: " << ev << std::endl;
    //applying Leiden algo
    Partition partition{std::vector<Community>{}};
    std::vector<Flat> flatFinalPartition;
    singletonPartition(graph, partition);
    std::cout << "INIITIAL GRAPH NODE SIZE " << partition.getCommunities().size() << std::endl;
    int gamma{1};
    double theta{0.01};
    leidenAlgorithm(graph, partition, flatFinalPartition, gamma, theta);
  }
}

void processRootFile(const std::string& fileName) {
  TFile* file = TFile::Open(fileName.c_str());
  if (!file || file->IsZombie()) {
    std::cerr << "Error opening file: " << fileName << std::endl;
    return;
  }

  // Navigate to the ticlDumper subdirectory
  TDirectory* dir = file->GetDirectory("ticlDumper");
  if (!dir) {
    std::cerr << "ticlDumper directory not found in file: " << fileName << std::endl;
    file->Close();
    return;
  }

  // Access the "tracksters" TTree
  TTree* tree = (TTree*)dir->Get("tracksters");
  if (!tree) {
    std::cerr << "tracksters TTree not found in file: " << fileName << std::endl;
    file->Close();
    return;
  }

  auto maxNumEvents = tree->GetEntries();
  size_t numOfEvents = 40;
  // Loop over the entries and fill the vector with Tracksters
  for (size_t ev = 0; ev < numOfEvents; ++ev) {
    std::cout << "ev " << ev << std::endl;
    processEvent(tree, ev);
  }

  file->Close();
}

// process list of root files
void readMultipleRootFiles(const std::vector<std::string>& fileNames) {
  for (const auto& fileName : fileNames) {
    processRootFile(fileName);
  }
}

int main(int argc, char** argv) {
  // Check if the path to the text file is provided
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <file_list.txt>" << std::endl;
    return 1;
  }

  // Read the file paths from the text file
  std::vector<std::string> fileNames = readFilePaths(argv[1]);
  if (fileNames.empty()) {
    std::cerr << "No valid file paths found in: " << argv[1] << std::endl;
    return 1;
  }
  readMultipleRootFiles(fileNames);

  return 0;
}