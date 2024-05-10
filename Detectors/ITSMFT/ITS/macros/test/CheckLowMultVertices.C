// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <vector>
#include <set>
#include <unordered_map>
#include <string>
#include <bits/stdc++.h>  // for popcount

#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>

// #include "CommonUtils/RootSerializableKeyValueStore.h"
#include "Framework/Logger.h"
#include "ITSBase/GeometryTGeo.h"
#include "ReconstructionDataFormats/Vertex.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/TrackReference.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITS/TrackITS.h"
#include "ITSMFTSimulation/Hit.h"
#endif
#include "DataFormatsITSMFT/CompCluster.h"

o2::MCCompLabel getMainLabel(std::vector<o2::MCCompLabel>& labs);

struct ParticleInfo {
  int event;
  int pdg;
  float pt;
  float eta;
  float phi;
  int mother;
  int first;
  unsigned short clusters = 0u;
  unsigned char isReco = 0u;
  unsigned char isFake = 0u;
  bool isPrimary = 0u;
  unsigned char storedStatus = 2; /// not stored = 2, fake = 1, good = 0
  bool canContribToVertex = false;
  std::array<int, 7> rofs = {-1, -1, -1, -1, -1, -1, -1}; /// readout frames of corresponding clusters
  o2::its::TrackITS track;
  o2::MCCompLabel lab;
};

struct TrackInfo {
  int nHits;
  double pT;
  double totMomentum;
  int maxLayer;
};

struct SimVertInfo {
  int nContributors;                   // (likely) number of tracklets associated to vertex
  int nCells;                          // number of cells associated to vertex
  bool wasRecod;                       // whether the vertex was reconstructed
  float x, y, z;                       // position of the vertex
  std::map<int, TrackInfo> trackInfo;  // number of hits associated to each trackID
};

struct RecoVertInfo {
  int eventID;        // ID of event (-1 if fake)
  int nContributors;  // (likely) number of cells associated to vertex
  float x, y, z;      // position of the vertex
  int iROF;           // in which ROF this vertex was reco'd
};

struct GlobalRecoInfo {  // this is only really for lowMultBeamDistCut changes
  int nLowMult;      // number of low mult events
  int nRecoLowMult;  // number of reconstructed low mult events
  int nZeros;        // number of zero track events
  double globalEff;  // global (average) efficiency
  double globalPur;  // global (average) purity
};

struct GlobalInfo {
  GlobalRecoInfo globalRecoInfo;           // Info to compare lowMultBeamDistCuts for
  std::vector<RecoVertInfo> recoVertInfo;  // Info on all reco'd vertices
  std::map<int, SimVertInfo> simVertInfo;  // Info on simVerts
};

struct RofInfo {
  void print();
  void uniqeff();
  int id = 0;
  std::vector<int> eventIds;                                                       // ID of events in rof
  std::vector<bool> usedIds;                                                       // EvtID used to calculate actual efficiency
  std::vector<ParticleInfo> parts;                                                 // Particle usable for vertexing
  std::vector<std::vector<o2::MCCompLabel>> vertLabels;                            // Labels associated to contributors to vertex
  std::unordered_map<int, std::array<double, 3>> simVerts;                         // Simulated vertices of events that can be spot in current rof <evtId, pos[3]>
  std::vector<o2::dataformats::Vertex<o2::dataformats::TimeStamp<int>>> recoVerts; // Vertices found in current ROF
  float recoeff = 0.f;                                                             // Vertexing efficiency
  float purity = 0.f;                                                              // Purity of reconstructed vertices
  float duplRate = 0.f;                                                            // Duplication rate of reconstructed vertices in ROF
  GlobalInfo* globalInfo;                                                          // Central Global info struct
};

void RofInfo::print()
{
  std::cout << "\n=================================== ROF " << id << " ============================================ \n";
  // Simulated vertices
  for (auto& sV : simVerts) {
    std::cout << "\tSimulated vertex for event: " << sV.first << " vertex:"
              << " x= " << sV.second[0]
              << " y= " << sV.second[1]
              << " z= " << sV.second[2]
              << std::endl;
    std::cout << "\t\tPotentially contributing tracks:\n";
    for (auto& part : parts) {
      if (part.lab.getEventID() == sV.first && part.canContribToVertex) {
        std::cout << "\t\t\t" << part.lab << "\t" << part.pt << " [GeV]\t" << part.pdg << std::endl;
      }
    }
    std::cout << std::endl;
  }

  // Reconstructed vertices
  for (size_t iV{0}; iV < recoVerts.size(); ++iV) {
    auto l = getMainLabel(vertLabels[iV]);
    auto eventID = l.isSet() ? l.getEventID() : -1;
    std::cout << "\tReconstructed vertex for event: " << eventID << " (-1: fake):"
              << " x= " << recoVerts[iV].getX()
              << " y= " << recoVerts[iV].getY()
              << " z= " << recoVerts[iV].getZ()
              << std::endl;
    std::cout << "\t\tContributor labels:\n";
    for (auto& l : vertLabels[iV]) {
      std::cout << "\t\t\t" << l << std::endl;
    }
  }

  // Efficiency
  if (simVerts.size()) {  // without simulated vertex efficiency means nothing
    std::cout << "\n\tEfficiency: " << recoeff * 100 << " %\n";
  } else {
    std::cout << "\n\tEfficiency: N/A\n";
  }
  // if (!recoeff) {
    // int sv = simVerts.begin()->first;  // usually just one vertex
    // std::cout << "Efficiency zero, ROF: " << id << "; Number of tracks in vertex: "
              // << std::count_if(parts.begin(), parts.end(), [](auto part) {part.lab.getEventID() == simVerts.begin()->first}) << std::endl;
  // Purity
  if (recoVerts.size()) {
    std::cout << "\tPurity: " << purity *  100 << "%\n";
    std::cout << "\tDuplication Rate: " << duplRate *  100 << "%\n";
    if (duplRate) {
      std::cout << "Duplication rate nonzero for rof " << id << std::endl;
    }
  }
}

void RofInfo::uniqeff()
{
  int c{0};
  int nFakes{0};
  int current{-42};
   // sorting at this point should be harmless.
   std::sort(parts.begin(), parts.end(),
             [](ParticleInfo& lp, ParticleInfo& rp) {return lp.lab.getEventID() > rp.lab.getEventID(); });
  for (auto& p : parts) {
    if (p.lab.getEventID() != current) {
      eventIds.push_back(p.lab.getEventID());
      current = p.lab.getEventID();
    }
  }

  usedIds.resize(eventIds.size(), false);
  std::set<int> uniqueRecoVertIDs;
  for (size_t iV{0}; iV < vertLabels.size(); ++iV) {
    auto label = getMainLabel(vertLabels[iV]);
    int eventID = label.getEventID();
    if (!label.isSet()) {
      nFakes++;  // unset label means unmatched
      eventID = -1;
    }
    uniqueRecoVertIDs.insert(eventID);
    for (size_t iEvID{0}; iEvID < eventIds.size(); ++iEvID) {
      if (eventIds[iEvID] == label.getEventID() && !usedIds[iEvID]) {
        usedIds[iEvID] = true;
        globalInfo->simVertInfo[ eventIds[iEvID] ].wasRecod = true;
        ++c;
        break;
      }
    }
  }
  recoeff = (float)c / (float)eventIds.size();
  purity = 1 - (nFakes / (float)vertLabels.size());
  duplRate = 1. - (float)uniqueRecoVertIDs.size() / (float)vertLabels.size();
}

#pragma link C++ class ParticleInfo + ;
#pragma link C++ class RofInfo + ;

o2::MCCompLabel getMainLabel(std::vector<o2::MCCompLabel>& labs)
{
  o2::MCCompLabel lab;
  size_t max_count = 0;
  for (size_t i = 0; i < labs.size(); i++) {
    size_t count = 1;
    for (size_t j = i + 1; j < labs.size(); j++) {
      if (labs[i] == labs[j] && (labs[i].isSet() && labs[j].isSet()))
        count++;
    }
    if (count > max_count)
      max_count = count;
  }

  if (max_count == 1) { // pick first valid label in case of no majority
    for (size_t i = 0; i < labs.size(); i++) {
      if (labs[i].isSet())
        return labs[i];
    }
  }

  for (size_t i = 0; i < labs.size(); i++) {
    size_t count = 1;
    for (size_t j = i + 1; j < labs.size(); j++)
      if (labs[i] == labs[j])
        count++;
    if (count == max_count)
      lab = labs[i];
  }
  return lab;
}

GlobalInfo CheckVerticesSingle(
  const int dumprof = -1, std::string path = "exp300-0-0_1-0.05/lowMultBeamDistCut-0/",
  std::string tracfilePath = "o2trac_its.root", std::string clusfilePath = "o2clus_its.root",
  std::string kinefilePath = "sgn_1_Kine.root", std::string itsHitFilePath = "o2sim_HitsITS.root")
{
  using Vertex = o2::dataformats::Vertex<o2::dataformats::TimeStamp<int>>;
  using namespace o2::dataformats;
  using namespace o2::itsmft;
  using namespace o2::its;

  // Geometry
  o2::base::GeometryManager::loadGeometry(path.data());
  auto gman = o2::its::GeometryTGeo::Instance();

  // MC tracks and event header
  TFile* kineFile = TFile::Open((path + kinefilePath).data());
  TTree* mcTree = (TTree*)kineFile->Get("o2sim");
  mcTree->SetBranchStatus("*", 0); // disable all branches
  mcTree->SetBranchStatus("MCEventHeader*", 1);
  mcTree->SetBranchStatus("MCTrack*", 1);

  std::vector<o2::MCTrack>* mcArr = nullptr;
  mcTree->SetBranchAddress("MCTrack", &mcArr);
  MCEventHeader* eventHeader = nullptr;
  mcTree->SetBranchAddress("MCEventHeader.", &eventHeader);

  // MC Hits
  TFile* itsHitFile = TFile::Open((path + itsHitFilePath).data());
  TTree* hitTree = (TTree*)itsHitFile->Get("o2sim");
  std::vector<o2::itsmft::Hit>* hitArr = nullptr;
  hitTree->SetBranchAddress("ITSHit", &hitArr);

  // Clusters
  TFile* clusFile = TFile::Open((path + clusfilePath).data());
  TTree* clusTree = (TTree*)clusFile->Get("o2sim");
  std::vector<CompClusterExt>* clusArr = nullptr;
  clusTree->SetBranchAddress("ITSClusterComp", &clusArr);
  std::vector<o2::itsmft::ROFRecord>* clusROFRecords = nullptr;
  clusTree->SetBranchAddress("ITSClustersROF", &clusROFRecords);

  // Cluster MC labels
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* clusLabArr = nullptr;
  clusTree->SetBranchAddress("ITSClusterMCTruth", &clusLabArr);

  // Reconstructed vertices
  TFile* tracFile = TFile::Open((path + tracfilePath).data());
  TTree* recTree = (TTree*)tracFile->Get("o2sim");

  std::vector<Vertex>* recVerArr = nullptr;
  recTree->SetBranchAddress("Vertices", &recVerArr);
  std::vector<ROFRecord>* recVerROFArr = nullptr;
  recTree->SetBranchAddress("VerticesROF", &recVerROFArr);
  std::vector<o2::MCCompLabel>* recLabelsArr = nullptr;
  recTree->SetBranchAddress("ITSVertexMCTruth", &recLabelsArr);
  // Process
  // Fill MC info
  auto nev{mcTree->GetEntriesFast()};
  std::vector<std::vector<ParticleInfo>> info(nev);
  std::vector<std::array<double, 3>> simVerts;
  GlobalInfo globalInfo;
  globalInfo.globalRecoInfo = {0, 0, 0, 0., 0.};  // initialise all to zero
  int prev_eventID = 0;
  for (auto n{0}; n < nev; ++n) {
    mcTree->GetEvent(n);
    hitTree->GetEvent(n);
    info[n].resize(mcArr->size());
    
    simVerts.push_back({eventHeader->GetX(), eventHeader->GetY(), eventHeader->GetZ()});
    globalInfo.simVertInfo.insert({(int) eventHeader->GetEventID() - 1,
                                    SimVertInfo{ eventHeader->GetNPrim(), 0, false,
                                                (float) eventHeader->GetX(),
                                                (float) eventHeader->GetY(),
                                                (float) eventHeader->GetZ(), {}}});
    
    for (unsigned int mcI{0}; mcI < mcArr->size(); ++mcI) {
      auto part = mcArr->at(mcI);
      info[n][mcI].event = n;
      info[n][mcI].pdg = part.GetPdgCode();
      info[n][mcI].pt = part.GetPt();
      info[n][mcI].phi = part.GetPhi();
      info[n][mcI].eta = part.GetEta();
      info[n][mcI].isPrimary = part.isPrimary();      
    }
    // now fill track info
    std::map<int, TrackInfo>* trInfoPtr = &(globalInfo.simVertInfo[eventHeader->GetEventID() - 1].trackInfo);
    
    for (int hitI = 0; hitI < hitArr->size(); hitI++) {
      o2::itsmft::Hit hit = hitArr->at(hitI);
      int trackID = hit.GetTrackID();
      int layer = gman->getLayer( hit.GetDetectorID() );
      // find where trackID would fit in the map
      std::map<int, TrackInfo>::iterator lowerBoundTrackID = trInfoPtr->lower_bound(trackID);
      // check if key already exist, then lowerBoundTrackID->first == trackID
      if (lowerBoundTrackID != trInfoPtr->end() && lowerBoundTrackID->first == trackID)
        lowerBoundTrackID->second.nHits++;
        // also increment maxLayer if further out
        if ( layer > lowerBoundTrackID->second.maxLayer)
          lowerBoundTrackID->second.maxLayer = layer;
      else  // key does not exist yet, set it, add a hit and the momentum
        trInfoPtr->insert(
          lowerBoundTrackID,
          std::pair<int, TrackInfo>(
            trackID, TrackInfo{1, mcArr->at(trackID).GetPt(), mcArr->at(trackID).GetP(), layer }));
    }
    if ((eventHeader->GetEventID() - prev_eventID) != 1)
      std::cout << "Event: " << eventHeader->GetEventID() << " not properly configured.\n";
    prev_eventID = eventHeader->GetEventID();
  }

  // Fill ROF info and complement MC info with cluster info
  std::vector<RofInfo> rofinfo;
  for (int frame = 0; frame < clusTree->GetEntriesFast(); frame++) { // Cluster frames
    if (!clusTree->GetEvent(frame))
      continue;
    rofinfo.resize(clusROFRecords->size());
    std::cout << "Size of rofinfo now: " << rofinfo.size() << "\n";
    for (size_t rof{0}; rof < clusROFRecords->size(); ++rof) {
      rofinfo[rof].globalInfo = &globalInfo;  // make sure all point to the same one
      for (int iClus{clusROFRecords->at(rof).getFirstEntry()}; iClus < clusROFRecords->at(rof).getFirstEntry() + clusROFRecords->at(rof).getNEntries(); ++iClus) {
        auto lab = (clusLabArr->getLabels(iClus))[0];
        if (!lab.isValid() || lab.getSourceID() != 0 || !lab.isCorrect())
          continue;

        int trackID, evID, srcID;
        bool fake;
        lab.get(trackID, evID, srcID, fake);
        if (evID < 0 || evID >= (int)info.size()) {
          std::cout << "Cluster MC label eventID out of range" << std::endl;
          continue;
        }
        if (trackID < 0 || trackID >= (int)info[evID].size()) {
          std::cout << "Cluster MC label trackID out of range" << std::endl;
          continue;
        }
        info[evID][trackID].lab = lab; // seems redundant but we are going to copy these info and loosing the nice evt/tr_id ordering
        const CompClusterExt& c = (*clusArr)[iClus];
        auto layer = gman->getLayer(c.getSensorID());
        info[evID][trackID].clusters |= 1 << layer;
        info[evID][trackID].rofs[layer] = rof;
      }
    }
  }
  //return;

  for (size_t evt{0}; evt < info.size(); ++evt) {
    auto& evInfo = info[evt];
    int ntrackable{0};
    int nusable{0};
    for (auto& part : evInfo) {
      if (part.clusters & (1 << 0) && part.clusters & (1 << 1) && part.clusters & (1 << 2)) {
        ++ntrackable;
        if (part.rofs[0] > -1 && part.rofs[0] == part.rofs[1] && part.rofs[1] == part.rofs[2]) {
          ++nusable;
          part.canContribToVertex = true;
          rofinfo[part.rofs[0]].parts.push_back(part);
          int trackID, evID, srcID;
          bool fake;
          part.lab.get(trackID, evID, srcID, fake);
          rofinfo[part.rofs[0]].simVerts[evID] = simVerts[evID];
          globalInfo.simVertInfo[evID].nCells++;  // increment number of contributors
        }
      }
    }
  }

  // Reco vertices processing
  for (int frame = 0; frame < recTree->GetEntriesFast(); frame++) { // Vertices frames
    if (!recTree->GetEvent(frame)) {
      continue;
    }
    std::cout << "Reco frame " << frame << " with size of recVerROFArr: " << recVerROFArr->size() << "\n";
    // loop on rof records
    int contLabIdx{0};
    for (size_t iRecord{0}; iRecord < recVerROFArr->size(); ++iRecord) {
      auto& rec = recVerROFArr->at(iRecord);
      auto verStartIdx = rec.getFirstEntry(), verSize = rec.getNEntries();
      rofinfo[iRecord].id = iRecord;
      rofinfo[iRecord].vertLabels.resize(verSize);
      int vertCounter{0};
      for (int iVertex{verStartIdx}; iVertex < verStartIdx + verSize; ++iVertex, ++vertCounter) {
        auto vert = recVerArr->at(iVertex);
        rofinfo[iRecord].recoVerts.push_back(vert);
        for (int ic{0}; ic < vert.getNContributors(); ++ic, ++contLabIdx) {
          rofinfo[iRecord].vertLabels[vertCounter].push_back(recLabelsArr->at(contLabIdx));
          // std::cout << "Pushed " << rofinfo[iRecord].vertLabels[vertCounter].back() << " at position " << rofinfo[iRecord].vertLabels[vertCounter].size() << std::endl;
        }
        auto label = getMainLabel(rofinfo[iRecord].vertLabels[vertCounter]);
        int eventID = label.getEventID();
        if (!label.isSet()) {
          eventID = -1;
        }
        globalInfo.recoVertInfo.push_back(
          RecoVertInfo{ eventID, vert.getNContributors(), vert.getX(), vert.getY(), vert.getZ(), 0 } );
      }
    }
  }
  // Epilog
  LOGP(info, "ROF inspection summary");
  size_t nvt{0}, nevts{0}, nrofSimFilled{0}, nrofRecoFilled{0};
  float addeff{0}, addPur{0};
  if (dumprof < 0) {
    for (size_t iROF{0}; iROF < rofinfo.size(); ++iROF) {
      auto& rof = rofinfo[iROF];
      nvt += rof.recoVerts.size();
      nevts += rof.simVerts.size();
      rof.uniqeff();
      if (rof.eventIds.size()) {
        addeff += rof.recoeff;
        nrofSimFilled++;
      }
      if (rof.recoVerts.size()) {
        addPur += rof.purity;
        nrofRecoFilled++;
      }
      // if (iROF) rof.print();
    }
  } else {
    rofinfo[dumprof].uniqeff();
    rofinfo[dumprof].print();
    addeff += rofinfo[dumprof].recoeff;
    addPur += rofinfo[dumprof].purity;
    nvt += rofinfo[dumprof].recoVerts.size();
    nevts += rofinfo[dumprof].simVerts.size();
  }

  int nReco = 0;
  for (const auto& [evtID, simVert] : globalInfo.simVertInfo) {
    nReco += simVert.wasRecod;
    globalInfo.globalRecoInfo.nZeros += simVert.nCells == 0;
    globalInfo.globalRecoInfo.nLowMult += simVert.nCells < 16;  // define low mult as <16 here
    globalInfo.globalRecoInfo.nRecoLowMult += (simVert.nCells < 16) * simVert.wasRecod;
  }
  globalInfo.globalRecoInfo.globalEff = (double)nReco / simVerts.size();
  globalInfo.globalRecoInfo.globalPur = (double)addPur / nrofRecoFilled;

  LOGP(info, "Summary:");
  LOGP(info, "Found {} vertices in {} usable out of {} simulated", nvt, nevts, simVerts.size());
  LOGP(info, "Average good vertexing efficiency: {}%", (addeff / (float)nrofSimFilled) * 100);
  LOGP(info, "Fraction of vertices reconstructed (not necessarily same ROF): {}%",
    globalInfo.globalRecoInfo.globalEff * 100);
  LOGP(info, "Average total vertexing purity: {}%", globalInfo.globalRecoInfo.globalPur * 100);    
  LOGP(info, "{} not reconstructed low multiplicity events, {} of which are zero.",
       (globalInfo.globalRecoInfo.nLowMult - globalInfo.globalRecoInfo.nRecoLowMult),
       globalInfo.globalRecoInfo.nZeros );
  
  std::cout << "Size of globalInfo: " << globalInfo.simVertInfo.size() << " and " <<
               globalInfo.recoVertInfo.size() << std::endl;

  // close all files again
  kineFile->Close();
  itsHitFile->Close();
  tracFile->Close();
  clusFile->Close();

  return globalInfo;
}

void CheckLowMultVertices(
  const int dumprof = -1, const float min=0.1, const float max=0.1, const float step_size = 0.05,
  const int nBatches = -1, std::string parent_dir = "exp100-0-1-0_05/", std::string tracfile = "o2trac_its.root",
  std::string clusfile = "o2clus_its.root", std::string kinefile = "sgn_1_Kine.root")
{
  // if nBatches is unset (-1), then use just parent dir. If set, the parent dir contains
  // a number (nBatches) of datasets, each just numbered from 0 to nBatches+1
  int nDatasets = (nBatches > 0) ? nBatches : 1;
  
  for (int batchNo = 0; batchNo < nDatasets; batchNo++) {

    std::vector<GlobalRecoInfo> globalRecoInfo;  // depends on lowMultBeamDistCut, i.e. will change
    GlobalInfo finalInfo;  // only for the largest (best) lowMultBeamDistCut
    unsigned int iterCount = 0;
    for (float curr=min; curr <= max; curr+=step_size, iterCount++ ) {
      std::string path = parent_dir;
      // add batch number to parent dir if we have batch system
      if (nBatches > 0) path += std::to_string(batchNo) + "/";

      path += "lowMultBeamDistCut-" + std::to_string(iterCount) + "/";
      GlobalInfo verticesInfo = CheckVerticesSingle(dumprof, path, tracfile, clusfile, kinefile);

      std::cout << "Size of globalInfo: " << verticesInfo.simVertInfo.size() << " and " <<
                verticesInfo.recoVertInfo.size() << std::endl;

      // go to next batch if nothing was reco'd. This means that the dataset is "corrupt".
      // If zero reconstructed vertices is expected as possible, then comment the next line.
      if (verticesInfo.recoVertInfo.size() == 0) break;

      globalRecoInfo.push_back(verticesInfo.globalRecoInfo);
      float eps = 1e-6;
      if ((max - curr - eps) < 0) finalInfo = verticesInfo;
    }
    // now create the output file and the corresponding tuples
    std::string outDir = parent_dir;
    if (nBatches > 0) outDir += std::to_string(batchNo) + "/";
    TFile* recoOutFile = new TFile((outDir + "reco_metadata.root").c_str(), "RECREATE", "Metadata");
  
    TNtuple* globalRecoInfoTuple = new TNtuple("lowMultRecoInfo", "Information wrt lowMultBeamDistCut",
                                            "lowMultBeamDistCut:nLowMult:nRecoLowMult:nZeroTrackEvents:"
                                            "eff:purity");

    TNtuple* simVertInfoTuple = new TNtuple("simVertInfo", "Information on simulated vertices",
                                            "eventID:nContributors:nCells:nTracks:wasReco:x:y:z");
    TNtuple* trackInfoTuple = new TNtuple("trackInfo", "Information on simulated tracks",
                                          "trackID:nHits:pT:totP:maxLayer");
    TNtuple* recoVertInfoTuple = new TNtuple("recoVertInfo", "Information on reconstructed vertices",
                                            "eventID:nContributors:x:y:z");

    // First fill averages wrt lowMultBeamDistCut
    for ( int iLMBDC = 0; iLMBDC < globalRecoInfo.size(); iLMBDC++ ) {
      GlobalRecoInfo vInfo = globalRecoInfo[iLMBDC];
      globalRecoInfoTuple->Fill( min + step_size*iLMBDC, vInfo.nLowMult, vInfo.nRecoLowMult,
                                 vInfo.nZeros, vInfo.globalEff, vInfo.globalPur );
    }
    // Now fill totals for the largest (best) lowMultBeamDistCut
    for (const auto& [evtID, sVInfo] : finalInfo.simVertInfo) {
      simVertInfoTuple->Fill( evtID, sVInfo.nContributors, sVInfo.nCells, sVInfo.trackInfo.size(),
                              sVInfo.wasRecod, sVInfo.x, sVInfo.y, sVInfo.z);
      
      // for each simVert, we have a lot of tracks/hits to fill...
      for (auto nHitIt = sVInfo.trackInfo.begin(); nHitIt != sVInfo.trackInfo.end(); nHitIt++) {
        trackInfoTuple->Fill( (double) nHitIt->first, (double) nHitIt->second.nHits,
                              nHitIt->second.pT, nHitIt->second.totMomentum,
                              nHitIt->second.maxLayer );
      }
    }
    // lastly, fill reco's info
    for (auto rVinfo : finalInfo.recoVertInfo ) {
      recoVertInfoTuple->Fill( rVinfo.eventID, rVinfo.nContributors,
                               rVinfo.x, rVinfo.y, rVinfo.z );
    }
    // write the tuples and file
    // globalRecoInfoTuple->Write();
    // simVertInfoTuple->Write();
    // recoVertInfoTuple->Write();
    // trackInfoTuple->Write();

    recoOutFile->Write();
    recoOutFile->Close();
  }
  
}