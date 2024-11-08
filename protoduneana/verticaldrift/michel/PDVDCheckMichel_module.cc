////////////////////////////////////////////////////////////////////////
// Class:       PDVDCheckMichel

// Plugin Type: analyzer (art v3_05_01)
// File:        PDVDCheckMichel_module.cc
//
// Generate simple summary tree from tracks and hits to check sim and reco
// Follow the same logic as in CosmicsdQdx
//
// Generated at 2/06/2023 by Thibaut Houdy
////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include <utility>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "art_root_io/TFileService.h"

// LArSoft includes
#include "lardataobj/Simulation/SimEnergyDeposit.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/Exceptions.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/ArtDataHelper/TrackUtils.h"

#include "TPolyMarker3D.h"
#include "TPolyMarker.h"
#include "TPolyLine.h"
#include "TPolyLine3D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TArc.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TLine.h"
#include "TLorentzVector.h"

#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"

// ROOT
#include "TTree.h"
#include "TLegend.h"

using std::cerr;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::pair;

//
namespace pdvdana {
  class PDVDCheckMichel;
}

class pdvdana::PDVDCheckMichel : public art::EDAnalyzer {
public:
  explicit PDVDCheckMichel(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PDVDCheckMichel(PDVDCheckMichel const&) = delete;
  PDVDCheckMichel(PDVDCheckMichel&&) = delete;
  PDVDCheckMichel& operator=(PDVDCheckMichel const&) = delete;
  PDVDCheckMichel& operator=(PDVDCheckMichel&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  int      fLogLevel;
  int      fFlagBasics = 1;
  int      fFlagInfos = 2;
  int      fFlagWarning = 3;
  int      fFlagDetails = 4;
// int      fFlagDebug = 5;

  string   fHitModuleLabel = "gaushit"; // or "hitpdune"
  string   fTrackModuleLabel = "pandoraTrack"; 
  string   fSpacePointLabel="reco3d"; // or "pandora"
  string   fSimuModuleLabel = "IonAndScint::G4Stage2"; //or "largeant:LArG4DetectorServicevolTPCActive"
  string   fClusterModuleLabel = "pandora"; // or "Reco3D"
  string   fGeneratorTag = "generator";
  float    fTrackMinLen;
  float    fMichelRadiusSphere=10;
  int      ev_num;

  int michel_counter=0, muon_michel_counter=0;

  unsigned fTotalTracks;
  unsigned fSelectedTracks;

  // summary tree
  TTree *ftrackTree;
  TTree *fhitTree;

  // Michel electron trees
  TTree *fmichelTree;
  // TTree *fcaloTree;
  
  protoana::ProtoDUNETruthUtils        truthUtil;
  protoana::ProtoDUNETrackUtils        trackUtil;
  protoana::ProtoDUNEPFParticleUtils    pfpUtil;

  //
  unsigned fEventNum;
  unsigned fTrackId;
  unsigned fNtpcs, fNplanes;
  unsigned ftrackEndTick;

  unsigned fMichTrackId;
  unsigned fMichNum;
  float fMothertrackLen;
  float fMichDistanceMuonEnd;
  vector<float> fMichX;
  vector<float> fMichY;
  vector<float> fMichZ;

  int fMichelPDG;
  float fMichelEne;
  bool fIsMichelInside;
  
  float fCaloResidualRange;
  float fCalodQdX;
  float fCalodEdX;
  float fCaloTrackLength;
  bool fCaloIsMuonMich;

//Saving track characteristics
  float ftrackLen;
  float ftrackDx;
  float ftrackDy;
  float ftrackDz;
  float ftrackStartX;
  float ftrackStartY;
  float ftrackStartZ;
  int ftrackEndWire;
  float ftrackEndX;
  float ftrackStartTime;
  float ftrackEndTime;
  float ftrackEndY;
  float ftrackEndZ;
  int ftrackStartWire;
  int ftrackCheatPDG;  
  float ftracktheta;
  float ftrackphi;
  float ftracknorm;

//Detector properties
  float fUsTimeFromTick;
  float fElectronVelocity;
  float fPitchCollection;


  float fhitsX;
  float fhitsY;
  float fhitsZ;
  float fhitscharge;
  float fhitsenergy;
  unsigned fNhits;


  vector<int> fOffsetWireID_u;
  vector<int> fOffsetWireID_v;
  vector<int> fOffsetWireID_z;
  
  float fgeoXmin = 1e6; 
  float fgeoXmax =-1e6; 
  float fgeoYmin = 1e6; 
  float fgeoYmax =-1e6; 
  float fgeoZmin = 1e6; 
  float fgeoZmax =-1e6; 
  bool  ftrackCheatIsMichMuon,ftrackIsInside, ftrackCheatIsMichMuon_Al;

  // detector geometry
  const geo::Geometry* fGeom;
  const string myname = "pdvdana::PDVDCheckMichel::analyze: ";

  void checkTrackCharacs(TVector3 l_track, vector<float> &track_char );
  bool IsThisPointInsideTPC(TVector3 PosPoint);
  bool IsThisTrackEnteringTPC(TVector3 PosPointA, TVector3 PosPointB);
  bool IsMichelMichel(simb::MCParticle l_particle, simb::MCParticle l_mother_particle);
  bool IsMuonMichel(simb::MCParticle l_particle);
  void DrawRectangle(TPolyLine *Rectangle, float a_1, float a_2, float b_1, float b_2);  
  pair<int,int> GetWireFromGeoPoint(geo::Point_t localpoint);


  void SetHitStyle(TPolyMarker *l_hits, int l_size, int l_color, int l_style);
  void SetPointStyle(TPolyMarker3D* l_pm, int l_size, int l_color, int l_style);
  void Drawing3D_AddingSpacePoints(TCanvas* l_c, vector<vector<double>> l_s, int l_size, int l_color, int l_style);
  void Drawing3D_HitsAndTracks(TCanvas* l_c,  vector<vector<double>> l_d);
  void Drawing2D_AddingSpacePoints(TCanvas* l_c, vector<vector<double>> l_s, int l_size, int l_color, int l_style);
  void Drawing2D_HitsAndTracks(TCanvas* l_c, vector<vector<double>> l_t, vector<vector<double>> l_v);
  void Drawing2D_yz_HitsAndTracks(TCanvas* l_c, vector<vector<double>> l_t, vector<vector<double>> l_v);
  void Drawing2D_yz_AddingSpacePoints(TCanvas* l_c, vector<vector<double>> l_s, int l_size, int l_color, int l_style);
  void DrawCube(TPolyLine3D *Cube, float x_min, float y_min, float z_min, float x_max, float y_max, float z_max );
  unsigned GetWireOffset(unsigned plane_id, unsigned tpc_id);

};


pdvdana::PDVDCheckMichel::PDVDCheckMichel(fhicl::ParameterSet const& p)
  : EDAnalyzer{p} ,
  fLogLevel( p.get< int >("LogLevel") ),
  fHitModuleLabel( p.get< std::string >("HitModuleLabel") ),
  fTrackModuleLabel( p.get< std::string >("TrackModuleLabel") ),
  fSpacePointLabel( p.get< std::string >("SpacePointLabel") ),
  fSimuModuleLabel( p.get< std::string >("SimuModuleLabel") ),
  fClusterModuleLabel( p.get< std::string >("ClusterModuleLabel") ),
  fGeneratorTag( p.get<std::string>("GeneratorTag") ),
  fTrackMinLen( p.get< float  >("TrackMinLen") ),
  fMichelRadiusSphere(p.get< float  >("MichelRadiusSphere") )
  {
    fGeom    = &*art::ServiceHandle<geo::Geometry>();
  }

//
void pdvdana::PDVDCheckMichel::analyze(art::Event const& ev)
{

  ev_num++;

  if (fLogLevel>=fFlagBasics)
    cout<<myname<< "Entering into the event: " << ev_num << endl;

  art::InputTag sim_tag(fSimuModuleLabel); //"IonAndScint::G4Stage2" or "largeant:LArG4DetectorServicevolTPCActive"
  
  art::Handle<std::vector<sim::SimEnergyDeposit>> simEnergyHandle;
  ev.getByLabel(fSimuModuleLabel, simEnergyHandle);

  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;

  vector<vector<double>> TracksSpacePoints;
  vector<vector<double>> MuTrackStartEnd;
  vector<vector<double>> MuMichSpacePoints;
  vector<vector<double>> AllHitsPosition;
  vector<vector<double>> AllHitsPosition_WT;

  vector<vector<double>> CheatPositionTrueElectron;
  vector<vector<double>> mich_candidate_interest;  


  vector<int> MichelElectronParticleTrackID;
  vector<int> MuonParticleTrackID;
  vector<int> trackID_identified_as_muons;

  int event_muon_counter, event_michel_electron_counter;
  event_michel_electron_counter=0;
  event_muon_counter=0;

  TH1F *histo_distance_center= new TH1F("h_dist_center", "dist [mm]",200 ,0 ,100 );
  TH1F *histo_energy_michel_depo= new TH1F("h_energy_center","dist [mm]",2000 ,0 ,10000 );


  // get tracks
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(ev);

  art::InputTag hit_tag(fHitModuleLabel); //"gaushit" or "hitpdune"
  art::InputTag trck_tag(fTrackModuleLabel); //"pandoraTrack"
  art::InputTag spacepoint_tag(fSpacePointLabel); //"pandora" or "reco3d"

  string fcaloLabel = "pandoraGnocalo";
  art::InputTag calo_tag(fcaloLabel);

  auto Tracks = ev.getValidHandle<vector<recob::Track>>(trck_tag);
  auto hitHandle = ev.getValidHandle<std::vector<recob::Hit>>(hit_tag);
  auto PFParticles = ev.getValidHandle<std::vector<recob::PFParticle>>("pandora");
  std::vector<art::Ptr<recob::Hit>> all_hits;
  art::fill_ptr_vector(all_hits, hitHandle);
  
  fTotalTracks += Tracks->size();
  fEventNum = ev_num;
  if (fLogLevel>=fFlagInfos){
    cout<<myname<<all_hits.size()<<" hits found"<<endl;
    cout<<myname<<Tracks->size()<<" tracks found"<<endl;
    cout<<myname<<PFParticles->size()<<" PFParticules found"<<endl;
    cout<<myname<< "All primary pfParticles = " <<  pfpUtil.GetNumberPrimaryPFParticle(ev,"pandora") << endl;
  }

  // Get the hits associated with the space points
  const art::FindManyP<recob::SpacePoint> fmsph(hitHandle, ev, spacepoint_tag);
  const art::FindManyP<recob::Track> tfh(hitHandle, ev, trck_tag);
  const art::FindManyP<recob::Hit> fmtht(Tracks, ev, trck_tag); 
  if (!fmsph.isValid()) {
    throw cet::exception("LArPandoraShowerCheatingAlg")
      << "Spacepoint and hit association not valid. Stopping.";
  }    
  if (!tfh.isValid()) {
    throw cet::exception("LArPandoraShowerCheatingAlg")
      << "Track and hit association not valid. Stopping.";
  }

  if (fLogLevel>=fFlagDetails)
    cout<<myname<<"Entering tracks loop"<<endl;

  //----------------------------------------------------------------------------------------------
  //RECO - 1. Loops on tracks
  //Target is to identify reco track as decaying muons within the TPC. The output is the location of the muon
  //decay spot to look for deposits from michel electrons there
  for (unsigned itrk = 0; itrk < Tracks->size(); ++itrk) {

    vector<int> daughterid_being_michel;
    daughterid_being_michel.clear();
    vector<float> track_charac(6); 

    const recob::Track& track = Tracks->at(itrk);
    ftrackLen    = track.Length();
    bool bool_is_this_track_a_decaying_muon = false;

    if (ftrackLen<fTrackMinLen)
      continue;

    trackID_identified_as_muons.push_back(track.ID());

    const simb::MCParticle *MCParticle_from_track = truthUtil.GetMCParticleFromRecoTrack(clockData, track, ev, fTrackModuleLabel);
    
    //To check later : how many tracks are eliminated at this step
    if(!MCParticle_from_track) 
      continue;


    bool HasElectron = false, HasNuMu = false, HasNuE = false;

    //If the particle generating the reconstructed track is a muon
    if(abs(MCParticle_from_track->PdgCode())==13){
      event_muon_counter++;
      //If the particle generating the reconstructed track is a muon and has generated daughters particles
      if((MCParticle_from_track->NumberDaughters())>0){       
        for (int ii=0; ii<(MCParticle_from_track->NumberDaughters());++ii){
          const simb::MCParticle* daughter1 = pi_serv->TrackIdToParticle_P((MCParticle_from_track->Daughter(ii)));
          if(!daughter1) 
            continue;
          int d_pdg1 = abs(daughter1->PdgCode());
         // double KE = (daughter1->E() - daughter1->Mass())*1e3; //Kinetic energy in MeV
          //if the daugther is an electron
          if (d_pdg1 == 11){
            HasElectron = true;
            if(IsMichelMichel(*daughter1, *MCParticle_from_track)){
              daughterid_being_michel.push_back(ii);
              for(size_t pos_i=0; pos_i<daughter1->NumberTrajectoryPoints(); pos_i++)
                CheatPositionTrueElectron.push_back({daughter1->Position(pos_i).X(),daughter1->Position(pos_i).Y(),daughter1->Position(pos_i).Z()});
            }
          }
          else if (d_pdg1 == 14) HasNuMu = true;
          else if (d_pdg1 == 12) HasNuE = true;
        }
      }
    }

    // Filling the track tree from basics from Pandora
    ftrackStartX  = track.Start().X();
    ftrackStartY  = track.Start().Y();
    ftrackStartZ  = track.Start().Z();
    ftrackEndX    = track.End().X();
    ftrackEndY    = track.End().Y();
    ftrackEndZ    = track.End().Z();
    TVector3 pos_init(ftrackStartX, ftrackStartY, ftrackStartZ);
    TVector3 pos_final(ftrackEndX, ftrackEndY, ftrackEndZ);
    checkTrackCharacs( pos_final-pos_init, track_charac );
    ftrackCheatPDG = MCParticle_from_track->PdgCode();
    ftrackCheatIsMichMuon_Al = (HasElectron and HasNuMu) and HasNuE;
    ftrackCheatIsMichMuon = IsMuonMichel(*MCParticle_from_track);
    ftrackIsInside = IsThisTrackEnteringTPC(pos_init, pos_final) and (IsThisPointInsideTPC(pos_final) or IsThisPointInsideTPC(pos_init));
    ftrackDx = track_charac[0];
    ftrackDy = track_charac[1];
    ftrackDz = track_charac[2];
    ftracknorm = track_charac[3];
    ftracktheta = track_charac[4];
    ftrackphi = track_charac[5];
    fTrackId     = 0;//track.key();

    // Extract track TIME
    //Loop over all the hits of the current track to find the min and the max of the time associated with the track
    std::vector<art::Ptr<recob::Hit>> this_track_hits=fmtht.at(itrk); 
    std::vector<float> hitpeakT;
    for(size_t idh1=0; idh1<this_track_hits.size(); idh1++){
      hitpeakT.push_back(this_track_hits[idh1]->PeakTime());
    }
    float max=-999999;
    float min=-999999;
    min=*min_element(hitpeakT.begin(),hitpeakT.end());
    max=*max_element(hitpeakT.begin(),hitpeakT.end());

    //Associating correctly the time of the start or the end of the track
    if(abs(ftrackStartX)>abs(ftrackEndX)){
      ftrackStartTime = min;
      ftrackEndTime = max;
    }
    else{
      ftrackStartTime = max;
      ftrackEndTime = min;
    }

    // Extract track WIRES
    geo::Point_t point_start_track = geo::Point_t(ftrackStartX, ftrackStartY, ftrackStartZ);
    geo::Point_t point_end_track = geo::Point_t(ftrackEndX, ftrackEndY, ftrackEndZ);     
    pair <int, int> trackStartWireTPC = GetWireFromGeoPoint(point_start_track);
    pair <int, int> trackEndWireTPC = GetWireFromGeoPoint(point_end_track);

    ftrackStartWire = trackStartWireTPC.second;
    ftrackEndWire = trackEndWireTPC.second;

    ftrackTree->Fill();
    //Track being candidate for a good michel muon track
    if(ftrackCheatIsMichMuon && ftrackIsInside){
      bool_is_this_track_a_decaying_muon = true;
      MuTrackStartEnd.push_back({ftrackStartX, ftrackStartY, ftrackStartZ, ftrackEndX, ftrackEndY, ftrackEndZ, 1.0});
      
      double KE_max = -0.1;

      //Filling the Michel electron tree
      for(auto electron_id:daughterid_being_michel){
        fMichX.clear();
        fMichY.clear();
        fMichZ.clear();
        const simb::MCParticle* daughter1 = pi_serv->TrackIdToParticle_P((MCParticle_from_track->Daughter(electron_id)));
        double KE = (daughter1->E() - daughter1->Mass())*1e3 ; //Kinetic energy in MeV
        fMichTrackId = daughter1->TrackId();
        fMothertrackLen = ftrackLen;
        for(size_t pos_i=0; pos_i<daughter1->NumberTrajectoryPoints(); pos_i++){
          fMichX.push_back(daughter1->Position(pos_i).X());
          fMichY.push_back(daughter1->Position(pos_i).Y());
          fMichZ.push_back(daughter1->Position(pos_i).Z());
        }
        fMichDistanceMuonEnd = sqrt(pow(daughter1->Position().Y()-ftrackEndY,2.0)+pow(daughter1->Position().Z()-ftrackEndZ,2.0));
        fMichNum = ev_num;
        fMichelPDG = daughter1->PdgCode();
        fMichelEne = KE;
        fmichelTree->Fill();
        if(KE>KE_max)
          KE_max=KE;
        michel_counter++;
        event_michel_electron_counter++;
        muon_michel_counter++;
      }       
      
      if(IsThisPointInsideTPC(pos_init))
        mich_candidate_interest.push_back({ftrackStartX, ftrackStartY, ftrackStartZ, ftrackStartTime, KE_max});        
      if(IsThisPointInsideTPC(pos_final))
        mich_candidate_interest.push_back({ftrackEndX, ftrackEndY, ftrackEndZ, ftrackEndTime, KE_max});      
    } //Track being candidate for being NOT a good michel muon track
    else{
      bool_is_this_track_a_decaying_muon = false;
      MuTrackStartEnd.push_back({ftrackStartX, ftrackStartY, ftrackStartZ, ftrackEndX, ftrackEndY, ftrackEndZ, 0.});
    }

    std::vector<anab::Calorimetry> calovector = trackUtil.GetRecoTrackCalorimetry(track, ev, fTrackModuleLabel, fcaloLabel);
    for (auto & calo : calovector){

      if (calo.PlaneID().Plane != 2) //only collection plane
        continue;

      for (size_t ihit = 0; ihit < calo.dQdx().size(); ++ihit){ //loop over hits
        fCaloResidualRange = calo.ResidualRange()[ihit];
        fCalodQdX = calo.dQdx()[ihit];
        fCalodEdX = calo.dEdx()[ihit];
        fCaloTrackLength = ftrackLen;
        fCaloIsMuonMich = bool_is_this_track_a_decaying_muon;
        // fcaloTree->Fill();       
        const auto &primtrk_pos=(calo.XYZ())[ihit];
        float local_hit_x = primtrk_pos.X();
        float local_hit_y = primtrk_pos.Y();
        float local_hit_z = primtrk_pos.Z();        
        TracksSpacePoints.push_back({local_hit_x, local_hit_y, local_hit_z});   
        } //End of hits loop from 1 calo
      
    }//End of calo loop coming from 1 track

  }// End of Tracks loop

  if (fLogLevel>=fFlagDetails)
    cout<<myname<<"Leaving tracks loop, scan of the hits done"<<endl;

  vector<float> lvec_hitsX, lvec_hitsY, lvec_hitsZ, lvec_hitscharge, lvec_energy;
  vector<unsigned> lvec_NumHits;

  lvec_hitsX.clear();
  lvec_hitsY.clear();
  lvec_hitsZ.clear();
  lvec_hitscharge.clear();
  lvec_NumHits.clear();
  lvec_energy.clear();
  for(size_t l_id =0; l_id<mich_candidate_interest.size(); l_id++){
    lvec_hitsX.push_back(0);
    lvec_hitsY.push_back(0); 
    lvec_hitsZ.push_back(0); 
    lvec_hitscharge.push_back(0);
    lvec_NumHits.push_back(0);
    lvec_energy.push_back(0);
  }

  //RECO - 2. Loops on hits
  //Target is to look for hits around the end point of tracks and integrate hits not associated with long tracks
  for (auto hit : all_hits) {
  
    //only checking hits from the collection plane
    if(!(hit->WireID().Plane ==2))
      continue;


    // Get the hit wire and time
    geo::WireID fWire;
    fWire = hit->WireID();
    auto const Zcent = fGeom->WireIDToWireGeo( fWire ).GetCenter();
    float Hit_Z = Zcent.Z();   
    AllHitsPosition_WT.push_back({Hit_Z, hit->PeakTime()});
  
    //checking if these hits are closed to local center (end or beginning of a Pandora track)
    for(size_t l_id =0; l_id<mich_candidate_interest.size(); l_id++){
      auto local_candidate = mich_candidate_interest[l_id];
      TVector3 local_center_xyz(local_candidate[0],local_candidate[1], local_candidate[2]);
      double local_center_time = local_candidate[3];
      // double local_center_wire = local_candidate[4];
      double local_KE = local_candidate[4];
      std::vector<art::Ptr<recob::Track>> trackfromhit = tfh.at(hit.key());

      //If a hit is identified as belonging to track
      bool thishit_is_from_muon_track = false;
      if(trackfromhit.size()>0){
        //loop over the different tracks the current hit belong to
        for(auto ltrack:trackfromhit){
          //loop over the track that are identified as potential muon (length> fTrackMinLen)
          for(int local_track_id_muon:trackID_identified_as_muons){
            if(ltrack->ID()==local_track_id_muon){
              thishit_is_from_muon_track=true;
              break;
            }
          }
          if(thishit_is_from_muon_track)
            break;
        }
      }
      if(thishit_is_from_muon_track)
        continue;
      

      //This is calculating the distance from identified michel electron center (muon decay)
    
      float distance_in_X = ((local_center_time-hit->PeakTime())*fUsTimeFromTick)*fElectronVelocity;
      float distance_in_Z = (local_center_xyz[2]-Hit_Z);
      float distance_WT_centre = sqrt(pow(distance_in_X,2.0)+pow(distance_in_Z,2.0));
      // double distance_centre = sqrt(pow(local_center_xyz[1]-local_hit_pos[1],2.0)+pow(local_center_xyz[2]-local_hit_pos[2],2.0));
      histo_distance_center->Fill(distance_WT_centre);
      if(distance_WT_centre < fMichelRadiusSphere){
        lvec_energy[l_id]=local_KE;
        lvec_hitsX[l_id]+=distance_in_X+local_center_xyz[0];
        lvec_hitsY[l_id]+=local_center_xyz[1];
        lvec_hitsZ[l_id]+=distance_in_Z+local_center_xyz[2];
        lvec_hitscharge[l_id]+=hit->SummedADC();
        histo_energy_michel_depo->Fill(hit->SummedADC());
        lvec_NumHits[l_id]++;
        MuMichSpacePoints.push_back({distance_in_X+local_center_xyz[0], local_center_xyz[1],distance_in_Z+local_center_xyz[2]});
      }
    }

 // Get hit position
    std::vector<art::Ptr<recob::SpacePoint>> sps = fmsph.at(hit.key());

    if(sps.size() < 1)
        continue; 
    
    art::Ptr<recob::SpacePoint> sp = sps.front();
    TVector3 local_hit_pos(sp->XYZ()[0], sp->XYZ()[1], sp->XYZ()[2]);

    AllHitsPosition.push_back({local_hit_pos[0], local_hit_pos[1], local_hit_pos[2]});

  }// End of the hits loop


//---------------------------------------------------------------------------------------------------------------------
//--------------------------------------ANALYSIS OF THE CANDIDATES-----------------------------------------------------

  art::InputTag itag2("wclsdatavd","gauss");
  auto wires = ev.getValidHandle<std::vector<recob::Wire> >(itag2);
  vector<TH2F*> wireData_histo;

  for(size_t l_id =0; l_id<mich_candidate_interest.size(); l_id++){
    
    auto michel_candidate = mich_candidate_interest[l_id];

    //Estimating the energy of each Michel candidate
    //Checking if around this center we found some hits and the SummedADC is not negative or null
    if(lvec_NumHits[l_id]>0 && lvec_hitscharge[l_id]>0){
      fhitsX = lvec_hitsX[l_id]/lvec_NumHits[l_id];
      fhitsY = lvec_hitsY[l_id]/lvec_NumHits[l_id];
      fhitsZ = lvec_hitsZ[l_id]/lvec_NumHits[l_id];
      fNhits = lvec_NumHits[l_id];
      fhitsenergy = lvec_energy[l_id];
      fhitscharge=lvec_hitscharge[l_id];
      fhitTree->Fill();
    }

    geo::Point_t localvertex = geo::Point_t(michel_candidate[0], michel_candidate[1], michel_candidate[2]);
    int localticktime = int(michel_candidate[3]);        

    //Creating the 2D plots for wire and time around each candidate
    geo::Point_t limit_inf = geo::Point_t(localvertex.X(), localvertex.Y(), localvertex.Z()-fMichelRadiusSphere);
    geo::Point_t limit_sup = geo::Point_t(localvertex.X(), localvertex.Y(), localvertex.Z()+fMichelRadiusSphere);

    pair<int, int> firstWire = GetWireFromGeoPoint(limit_inf);
    pair<int, int> lastWire = GetWireFromGeoPoint(limit_sup);
    if(firstWire.first != lastWire.first)
    {
      if(fLogLevel>=fFlagWarning)
        cout<<myname<<"Problem with plotting on a window time/wire if wire from different TPC"<<endl;
      continue;
    }

    int nselecTicks = 2000; // 500 ns per tick
    int nselecWires = abs(firstWire.second - lastWire.second); //0.5 cm pitch
    int firstTimeTick = int(localticktime-nselecTicks/2);
    int lastTimeTick = firstTimeTick+nselecTicks;

   // bool b_touched = false;

    size_t ncal = wires->size();

    /* Now do the loop over the wires. */
    TH2F* local_wire_histo = new TH2F("histo_WT", "histo_WT", ncal, 0, ncal, 100, 0, 13000);
    cout<<myname<<"TEST CANVAS Creating histo :"<< nselecWires<<"  "<<  firstWire.second<<"  "<< lastWire.second<<"  "<< nselecTicks<<"  "<< firstTimeTick<<"  "<< lastTimeTick<<"  "<<endl;

    if(wires){
      cout<<myname<<"ANALYSIS : Entering into the recob::wires"<<endl;
      for (size_t ichan=0;ichan<ncal; ++ichan){
        cout<<myname<<"TEST here there is a channel ("<<ichan<<" ; "<<(wires->at(ichan)).Channel()<<"  )" <<endl;
        vector<float> wire_signal = (wires->at(ichan)).Signal();
        cout<<myname<<"   comparing size of vectors : "<<wire_signal.size()<<"  "<<ncal<<"  "<<ichan<<endl;
        for (size_t iadc=0;iadc<wire_signal.size(); ++iadc){
          if(wire_signal[iadc]!=0)
            cout<<myname<<"TEST here there is signal ("<<iadc<<" ; "<<(wires->at(ichan)).Channel()<<") -> "<<wire_signal[iadc]<<endl;
          local_wire_histo->Fill((wires->at(ichan)).Channel(),iadc,wire_signal[iadc]);
        }
      }
    }
    wireData_histo.push_back(local_wire_histo);
    cout<<myname<<"TEST looking at wireData "<< wireData_histo.size()<<"  "<< local_wire_histo[10][10]<<endl;

/*    for (int iwire = 0; iwire < nselecWires; iwire++) {

    b_touched = false;
      //looking for wire signal corresponding to the actual iwire (iwire+firstWire.second)



    for (auto & wire : * wires){
        
        if (wire.View() != 2)
          continue;
        
        if (b_touched)
          continue;

        std::vector<geo::WireID> wireids = fGeom->ChannelToWire(wire.Channel());

        if(int(wireids[0].Wire) == int(iwire+firstWire.second) && !b_touched){ 
          cout<<myname<<"TEST WHICH WIRE "<< wireids[0].Wire <<" ---  "<< wire.Channel()<<" ---  "<< iwire+firstWire.second<<endl;
          b_touched=true;
          for(int i_time = 0; i_time < nselecTicks; i_time++){
            if((wire.Signal()).at(i_time+firstTimeTick)!=0)
              cout<<myname<<"TEST WIRE SIGNAL "<< (wire.Signal()).at(i_time+firstTimeTick) <<endl;
            wireData_histo[l_id]->SetBinContent(iwire,i_time,(wire.Signal()).at(i_time+firstTimeTick));
          }
          break;
        }
      }//end of the wire signal loop

    }// end of the loop to find corresponding wires*/
  
  } //End of the identified centers loop

  
//-----------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------
  
  if (fLogLevel>=fFlagDetails)
    cout<<myname<<"End of the analysis - Start of the plotting section"<<endl;


  if (fLogLevel>=fFlagDetails || ev_num<5){

    art::ServiceHandle<art::TFileService> tfs;
    gStyle->SetOptStat(0);
    //Drawing 3D tracks
    TString canvasName_3D = Form("canvas3D_%i", ev_num);
    TCanvas* canvas_3D = tfs->make<TCanvas>(canvasName_3D, canvasName_3D);

    //Drawing all the space points
    // MuTrackStartEnd -> blue thick line if mich muon candidate, grey ok line if not mich muon candidate 
    Drawing3D_HitsAndTracks(canvas_3D, MuTrackStartEnd);
    // AllHitsPosition -> markersize 1, color : rouille
    Drawing3D_AddingSpacePoints(canvas_3D,AllHitsPosition, 0.5, 45, 8);
    // True position of electrons in the simulation -> markersize 1, color : pink
    Drawing3D_AddingSpacePoints(canvas_3D, CheatPositionTrueElectron, 1, 6, 8);  
    // Position of mich electron candidates -> markersize 2, color : green
    Drawing3D_AddingSpacePoints(canvas_3D, MuMichSpacePoints, 1, 3, 4);  
    canvas_3D->Write(canvasName_3D);
    cout<<myname<<"Done with plotting 3D"<<endl;
    TString canvasName_2D = Form("canvas2D_%i", ev_num);
    TCanvas* canvas_2D = tfs->make<TCanvas>(canvasName_2D, canvasName_2D);

    //Drawing all the space points in 3 planes (xy,yz,xz)
    // MuTrackStartEnd -> blue thick line if mich muon candidate, grey ok line if not mich muon candidate 
    Drawing2D_HitsAndTracks(canvas_2D, MuTrackStartEnd, mich_candidate_interest);
    // AllHitsPosition -> markersize 0.5, color : grey
    Drawing2D_AddingSpacePoints(canvas_2D, AllHitsPosition, 0.5, 45, 8);
    // True position of electrons in the simulation -> markersize 0.5, color : pink
    Drawing2D_AddingSpacePoints(canvas_2D, CheatPositionTrueElectron, 0.5, 6 , 8);
    // Position of mich electron candidates -> markersize 2, color : green
    Drawing2D_AddingSpacePoints(canvas_2D, MuMichSpacePoints, 1, 3, 8);
    canvas_2D->Write(canvasName_2D);
    
    TString canvasName_2D_yz = Form("canvas2D_%i_yz", ev_num);
    TCanvas* canvas_2D_yz = tfs->make<TCanvas>(canvasName_2D_yz, canvasName_2D_yz);
    //Drawing all the space points in 3 planes (xy,yz,xz)
    // MuTrackStartEnd -> blue thick line if mich muon candidate, grey ok line if not mich muon candidate 
    Drawing2D_yz_HitsAndTracks(canvas_2D_yz, MuTrackStartEnd, mich_candidate_interest);
    // AllHitsPosition -> markersize 0.5, color : grey
    Drawing2D_yz_AddingSpacePoints(canvas_2D_yz, AllHitsPosition, 0.5, 45, 8);
    // True position of electrons in the simulation -> markersize 0.5, color : pink
    Drawing2D_yz_AddingSpacePoints(canvas_2D_yz, CheatPositionTrueElectron, 0.5, 6 , 8);
    // Position of mich electron candidates -> markersize 2, color : green
    Drawing2D_yz_AddingSpacePoints(canvas_2D_yz, MuMichSpacePoints, 1, 3, 8);
    canvas_2D_yz->Write(canvasName_2D_yz);
    cout<<myname<<"Done with plotting 2D"<<endl;

    //wireData_histo[0]->Write("histo_hardcoded");
    /*int mich_num=1;
    for(auto fhisto : wireData_histo){
      if(!fhisto)
        break;
      TString histo_name = Form("histo_candidate_%i_%i", ev_num, mich_num);
      cout<<myname<<"Plotting "<< histo_name<<endl;
      fhisto->Write(histo_name);
      mich_num++;
    }*/

    TString histo_distance_center_name = Form("histo_distance_center_%i", ev_num);
    histo_distance_center->Write(histo_distance_center_name);
    TString histo_energy_michel_depo_name = Form("histo_energy_michel_depo%i", ev_num);
    histo_energy_michel_depo->Write(histo_energy_michel_depo_name);


  }

  if (fLogLevel>=fFlagDetails){
    std::cout << myname<<" REPORT - in this event - # muons : " << event_muon_counter <<"  # michel electrons : " 
    << event_michel_electron_counter << " for a total of : "<< michel_counter <<" michel e"<<std::endl;
  }

  TracksSpacePoints.clear();
  MuTrackStartEnd.clear();
  MuMichSpacePoints.clear();
  AllHitsPosition.clear();
  CheatPositionTrueElectron.clear();
  mich_candidate_interest.clear(); 
  MichelElectronParticleTrackID.clear();
  MuonParticleTrackID.clear();
  trackID_identified_as_muons.clear();
  wireData_histo.clear();


} 
//
void pdvdana::PDVDCheckMichel::beginJob()
{

  ev_num = 0;
  fTotalTracks    = 0;
  fSelectedTracks = 0;
  fNplanes = fGeom->Nplanes();
  fNtpcs = fGeom->NTPC();
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clockData);
  fUsTimeFromTick  = 1e-3*sampling_rate(clockData); //from ns to us
  fElectronVelocity   = detProp.DriftVelocity(); //velocity in cm/us
  fPitchCollection = 5.; //(in mm)

  fOffsetWireID_u.clear();
  fOffsetWireID_v.clear();
  fOffsetWireID_z.clear();

  unsigned int mem_nwires_u = 0;
  unsigned int mem_nwires_v = 0;
  unsigned int mem_nwires_z = 0;
  
  if(fLogLevel>=fFlagBasics){
    cout<<myname<<"  DETECTOR GEOMETRY  : "<<endl;
    cout<<myname<<"number of planes :  "<<fNplanes<<endl;
    cout<<myname<<"number of tpc :  "<<fNtpcs<<endl;
    cout<<myname<<"pitch in collection plane used here : "<<fPitchCollection<<endl;
    cout<<myname<<"drift velocity : "<< fElectronVelocity <<"  with a tick time in us :  "<<fUsTimeFromTick << endl;
  }
  
  for(unsigned t_tpc_id=0;t_tpc_id<fNtpcs;t_tpc_id++){

    geo::TPCID tpcid{0, t_tpc_id};
    geo::PlaneID const uplane_id{tpcid, geo::View_t::kU};
    geo::PlaneID const vplane_id{tpcid, geo::View_t::kV};
    geo::PlaneID const zplane_id{tpcid, geo::View_t::kZ};

    unsigned int nwires_u = fGeom->Nwires(uplane_id);
    unsigned int nwires_v = fGeom->Nwires(vplane_id);
    unsigned int nwires_z = fGeom->Nwires(zplane_id);

    fOffsetWireID_u.push_back(mem_nwires_u);
    fOffsetWireID_v.push_back(mem_nwires_v);
    fOffsetWireID_z.push_back(mem_nwires_z);

    mem_nwires_u += nwires_u;
    mem_nwires_v += nwires_v;
    mem_nwires_z += nwires_z;

    if(fLogLevel>=fFlagDetails){
      std::cout <<myname<< "  TPC " << t_tpc_id << " center: ("<< fGeom->TPC(tpcid).GetCenter().X()      << "," << fGeom->TPC(tpcid).GetCenter().Y()      << ","<< fGeom->TPC(tpcid).GetCenter().Z() << ")"
                                 << " box:  ["  << fGeom->TPC(tpcid).BoundingBox().MinX() << "," << fGeom->TPC(tpcid).BoundingBox().MaxX() << "]" 
                                        << "["  << fGeom->TPC(tpcid).BoundingBox().MinY() << "," << fGeom->TPC(tpcid).BoundingBox().MaxY() << "]"
                                        << "["  << fGeom->TPC(tpcid).BoundingBox().MinZ() << "," << fGeom->TPC(tpcid).BoundingBox().MaxZ() << "]" << std::endl;
      cout<<myname<< "Number of wires : "<<nwires_u<<"  "<<nwires_v <<"  "<<nwires_z<<endl;
    }

    if(fgeoXmin > fGeom->TPC(tpcid).BoundingBox().MinX()) fgeoXmin = fGeom->TPC(tpcid).BoundingBox().MinX();
    if(fgeoXmax < fGeom->TPC(tpcid).BoundingBox().MaxX()) fgeoXmax = fGeom->TPC(tpcid).BoundingBox().MaxX();
    if(fgeoYmin > fGeom->TPC(tpcid).BoundingBox().MinY()) fgeoYmin = fGeom->TPC(tpcid).BoundingBox().MinY();
    if(fgeoYmax < fGeom->TPC(tpcid).BoundingBox().MaxY()) fgeoYmax = fGeom->TPC(tpcid).BoundingBox().MaxY();
    if(fgeoZmin > fGeom->TPC(tpcid).BoundingBox().MinZ()) fgeoZmin = fGeom->TPC(tpcid).BoundingBox().MinZ();
    if(fgeoZmax < fGeom->TPC(tpcid).BoundingBox().MaxZ()) fgeoZmax = fGeom->TPC(tpcid).BoundingBox().MaxZ();

  }
    if(fLogLevel>=fFlagBasics){
      std::cout <<myname<< "  Total LAr actif volume is : ["<<fgeoXmin<<"; " <<fgeoXmax<<"] in X, "
      <<"["<<fgeoYmin<<"; " <<fgeoYmax<<"] in Y, "
      <<"["<<fgeoZmin<<"; " <<fgeoZmax<<"] in Z"
      <<endl;
    }

  // init summary tree
  art::ServiceHandle<art::TFileService> tfs;

  fmichelTree = tfs->make<TTree>("michel","Michel candidates");
  fmichelTree->Branch("EventNum", &fMichNum);
  fmichelTree->Branch("TrackId", &fMichTrackId);
  fmichelTree->Branch("MotherTrackLen",   &fMothertrackLen);
  fmichelTree->Branch("PosX",  &fMichX);
  fmichelTree->Branch("PosY",  &fMichY);
  fmichelTree->Branch("PosZ",  &fMichZ);
  fmichelTree->Branch("PDG",  &fMichelPDG);
  fmichelTree->Branch("Energy", &fMichelEne);
  fmichelTree->Branch("DistanceMuonEnd", &fMichDistanceMuonEnd);
  fmichelTree->Branch("IsInside", &fIsMichelInside);

  ftrackTree = tfs->make<TTree>("tracks","Check tracks");
  ftrackTree->Branch("EventNum", &fEventNum, "EventNum/i");
  ftrackTree->Branch("TrackId", &fTrackId, "TrackId/i");
  ftrackTree->Branch("TrackLen",   &ftrackLen,   "TrackLen/F");
  ftrackTree->Branch("Dx",  &ftrackDx,  "Dx/F");
  ftrackTree->Branch("Dy",  &ftrackDy,  "Dy/F");
  ftrackTree->Branch("Dz",  &ftrackDz,  "Dz/F");
  ftrackTree->Branch("StartX",  &ftrackStartX,  "StartX/F");
  ftrackTree->Branch("StartY",  &ftrackStartY,  "StartY/F");
  ftrackTree->Branch("StartZ",  &ftrackStartZ,  "StartZ/F");  
  ftrackTree->Branch("StartWire",  &ftrackStartWire);
  ftrackTree->Branch("StartTime",  &ftrackStartTime);
  ftrackTree->Branch("EndX",  &ftrackEndX,  "EndX/F");
  ftrackTree->Branch("EndY",  &ftrackEndY,  "EndY/F");
  ftrackTree->Branch("EndZ",  &ftrackEndZ,  "EndZ/F");
  ftrackTree->Branch("EndWire",  &ftrackEndWire);
  ftrackTree->Branch("EndTime",  &ftrackEndTime);
  ftrackTree->Branch("theta",&ftracktheta,  "theta/F");
  ftrackTree->Branch("phi",  &ftrackphi,  "phi/F");
  ftrackTree->Branch("norm", &ftracknorm,  "norm/F");
  ftrackTree->Branch("trackCheatPDG", &ftrackCheatPDG);
  ftrackTree->Branch("trackCheatIsMichMuon_Al", &ftrackCheatIsMichMuon_Al);
  ftrackTree->Branch("trackCheatIsMichMuon", &ftrackCheatIsMichMuon);
  ftrackTree->Branch("trackIsInside", &ftrackIsInside);


  fhitTree  = tfs->make<TTree>("hits","Check reconstruction");
  fhitTree->Branch("X", &fhitsX, "X/F");
  fhitTree->Branch("Y", &fhitsY, "Y/F");
  fhitTree->Branch("Z", &fhitsZ, "Z/F");
  fhitTree->Branch("Nhits", &fNhits);
  fhitTree->Branch("charge", &fhitscharge);
  fhitTree->Branch("energy", &fhitsenergy);
  fhitTree->Branch("EventNum", &fEventNum, "EventNum/i");
  fhitTree->Branch("TrackId",  &fTrackId, "TrackId/i");


}

//
void pdvdana::PDVDCheckMichel::endJob()
{

  if(fLogLevel>=fFlagBasics){
    cout<<myname<<"tracks processed total     : "<<fTotalTracks<<endl;
    cout<<myname<<"tracks processed selected  : "<<fSelectedTracks<<endl;
  }
}
void pdvdana::PDVDCheckMichel::DrawCube(TPolyLine3D *Cube, float x_min, float y_min, float z_min, float x_max, float y_max, float z_max ){
  Cube->SetPoint(0,x_min, y_min, z_min);
  Cube->SetPoint(1,x_max, y_min, z_min);
  Cube->SetPoint(2,x_max, y_max, z_min);
  Cube->SetPoint(3,x_min, y_max, z_min);
  Cube->SetPoint(4,x_min, y_min, z_min);
  Cube->SetPoint(5,x_min, y_min, z_max);
  Cube->SetPoint(6,x_max, y_min, z_max);
  Cube->SetPoint(7,x_max, y_min, z_min);
  Cube->SetPoint(8,x_max, y_max, z_min);
  Cube->SetPoint(9,x_max, y_max, z_max);
  Cube->SetPoint(10,x_min, y_max, z_max);
  Cube->SetPoint(11,x_min, y_min, z_max);
  Cube->SetPoint(12,x_min, y_max, z_max);    
  Cube->SetPoint(13,x_min, y_max, z_min);
  Cube->SetPoint(14,x_min, y_max, z_max);
  Cube->SetPoint(15,x_max, y_max, z_max);
  Cube->SetPoint(16,x_max, y_min, z_max);
  Cube->SetLineColor(12);
  Cube->SetLineWidth(3);
  return;
}

void pdvdana::PDVDCheckMichel::SetPointStyle(TPolyMarker3D* l_PolyMarker, int l_size, int l_color, int l_style){
  l_PolyMarker->SetMarkerSize(l_size);    
  l_PolyMarker->SetMarkerColor(l_color);
  l_PolyMarker->SetMarkerStyle(l_style);
}
void pdvdana::PDVDCheckMichel::SetHitStyle(TPolyMarker* l_PolyMarker, int l_size, int l_color, int l_style){
  l_PolyMarker->SetMarkerSize(l_size);    
  l_PolyMarker->SetMarkerColor(l_color);
  l_PolyMarker->SetMarkerStyle(l_style);
}
void pdvdana::PDVDCheckMichel::Drawing3D_AddingSpacePoints(TCanvas* canvas_3D, vector<vector<double>> l_SpacePoints, int l_size=1, int l_color=2, int l_style=8){

  TPolyMarker3D* SpacePoints = new TPolyMarker3D(l_SpacePoints.size());
  
  int actual_point=0;
  for (auto const& sp : l_SpacePoints) {
      SpacePoints->SetPoint(actual_point, sp[0], sp[1], sp[2]);
      actual_point++;
  }

  SetPointStyle(SpacePoints, l_size, l_color, l_style);    
  canvas_3D->cd();
  SpacePoints->Draw();
}


void pdvdana::PDVDCheckMichel::Drawing3D_HitsAndTracks(TCanvas* canvas_3D, vector<vector<double>> l_trackStartEnd){
  //Drawing 3D tracks

  TH3F *axes = new TH3F("axes", "3D Hits and tracks distribution; X; Y; Z", 1, -400, 400, 1, -400, 400, 1, -10, 400);
  axes->SetDirectory(0);
  canvas_3D->cd();
  axes->Draw();

  //Drawing the TPC boxes (to be understood in a pandora/LArSoft TPC definition)
  for(unsigned t_tpc_id=0;t_tpc_id<fNtpcs;t_tpc_id++){
    geo::TPCID tpcid{{0}, t_tpc_id};
    TPolyLine3D* TPC_3D = new TPolyLine3D(17);
    DrawCube(TPC_3D, fGeom->TPC(tpcid).BoundingBox().MinX(), fGeom->TPC(tpcid).BoundingBox().MinY(), fGeom->TPC(tpcid).BoundingBox().MinZ(), 
            fGeom->TPC(tpcid).BoundingBox().MaxX(), fGeom->TPC(tpcid).BoundingBox().MaxY(), fGeom->TPC(tpcid).BoundingBox().MaxZ());
    TPC_3D->SetLineWidth(1);
    TPC_3D->SetLineColor(29);
    canvas_3D->cd();
    TPC_3D->Draw();
  }

  //Drawing the active detector 
  TPolyLine3D* Detector3D = new TPolyLine3D(17);
  DrawCube(Detector3D, fgeoXmin, fgeoYmin, fgeoZmin, fgeoXmax, fgeoYmax, fgeoZmax);
  canvas_3D->cd();
  Detector3D->Draw();
  
  //Drawing the tracks
  for (auto const& sp : l_trackStartEnd) {
    TPolyLine3D* t_line = new TPolyLine3D(2);
    t_line->SetPoint(0, sp[0], sp[1], sp[2]);
    t_line->SetPoint(1, sp[3], sp[4], sp[5]);
    //Here the tracks that are michel muon candidates
    if(sp[6]>0.5){
      t_line->SetLineColor(4);
      t_line->SetLineWidth(2);
    }//Here the tracks that are not michel muon candidates
    else{
      t_line->SetLineColor(15);
      t_line->SetLineWidth(1);
    }
    t_line->Draw();
  }

  //Adding a legend associated with tracks (mich muon candidate)
  TH1 *h_leg_dm = new TH1F("", "", 1, 0, 1);
  h_leg_dm->SetLineColor(4);
  h_leg_dm->SetLineWidth(2);
  
  //Adding a legend associated with tracks (not mich muon candidate)
  TH1 *h_leg_ndm = new TH1F("", "", 1, 0, 1);
  h_leg_ndm->SetLineColor(15);
  h_leg_ndm->SetLineWidth(1);

  //Adding a legend associated with electron position from simulation
  TH1F *h_leg_sp_mc_elec = new TH1F("", "", 1, 0, 1);
  h_leg_sp_mc_elec->SetMarkerStyle(8);
  h_leg_sp_mc_elec->SetMarkerSize(1);    
  h_leg_sp_mc_elec->SetMarkerColor(6);

  //Adding a legend associated with space points that are electron michel candidate
  TH1F *h_leg_sp_michmuon = new TH1F("", "", 1, 0, 1);
  h_leg_sp_michmuon->SetMarkerStyle(4);
  h_leg_sp_michmuon->SetMarkerSize(1);    
  h_leg_sp_michmuon->SetMarkerColor(3);

  //Adding a legend associated with space points that are associated with conduction hit
  TH1F *h_leg_sp_allhits = new TH1F("", "", 1, 0, 1);
  h_leg_sp_allhits->SetMarkerStyle(8);
  h_leg_sp_allhits->SetMarkerSize(1);    
  h_leg_sp_allhits->SetMarkerColor(15);

  auto legend = new TLegend(0.1,0.75,0.4,0.9);
  //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
  legend->AddEntry(h_leg_sp_allhits,"Reco SpacePoint of all hits","p");
  legend->AddEntry(h_leg_ndm,"Reco Tracks (sel : lenght), ","l");
  legend->AddEntry(h_leg_dm,"Reco tracks identified as decaying muons (cheat)","l");
  legend->AddEntry(h_leg_sp_mc_elec,"True position of michel electrons (cheat)","p");
  legend->AddEntry(h_leg_sp_michmuon,"Reco SP of michel electron candidates","p");
  legend->Draw();

  return;
} 
void pdvdana::PDVDCheckMichel::Drawing2D_AddingSpacePoints(TCanvas* canvas_2D, vector<vector<double>> l_SpacePoints, int l_size=1, int l_color=2, int l_style=8){

  vector<TPolyMarker*> l_SpacePoints_v(3);  
  for (int i =0; i<3; i++){
    l_SpacePoints_v[i] = new TPolyMarker(l_SpacePoints.size());
    int actual_point=0;
    int idx = int(i/2);
    for (auto const& sp : l_SpacePoints) {
        l_SpacePoints_v[i]->SetPoint(actual_point, sp[idx], sp[int(i+1-idx)]);
        actual_point++;
    }
  canvas_2D->cd(i+1);
  SetHitStyle(l_SpacePoints_v[i], l_size, l_color, l_style);      
  l_SpacePoints_v[i]->Draw();
  }
  return;
}


//Drawing 2D tracks
void pdvdana::PDVDCheckMichel::Drawing2D_HitsAndTracks(TCanvas* canvas2D_yz, vector<vector<double>> l_trackStartEnd, vector<vector<double>> l_mich_centers){


  vector<TPolyLine*> Rectangle(3);
  vector<TPolyLine*> Rectangle_tpc(3);
  vector<TPolyLine*> Tracks_lines(3);

  for (int i =0; i<3; i++){
      Rectangle[i] = new TPolyLine(5);
      Rectangle_tpc[i] = new TPolyLine(5);
  }

  DrawRectangle(Rectangle[0], fgeoXmin, fgeoXmax, fgeoYmin, fgeoYmax);
  DrawRectangle(Rectangle[1], fgeoXmin, fgeoXmax, fgeoZmin, fgeoZmax);
  DrawRectangle(Rectangle[2], fgeoYmin, fgeoYmax, fgeoZmin, fgeoZmax);

  vector<TH2F*> axos_2D(4);
  axos_2D[0] = new TH2F("axes_1", "Y(X); X; Y", 1, -400, 400, 1, -400, 400);
  axos_2D[1] = new TH2F("axes_2", "Z(X); X; Z", 1, -400, 400, 1, -10, 400);
  axos_2D[2] = new TH2F("axes_3", "Z(Y); Y; Z", 1, -400, 400, 1, -10, 400);
  axos_2D[3] = new TH2F("axes_4", "", 1, 0, 1, 1, 0, 1);

  canvas2D_yz->Divide(2,2);
  for(int i=0; i<3; i++){
      canvas2D_yz->cd(i+1);
      axos_2D[i]->Draw();
      Rectangle[i]->Draw();
  }

  for(unsigned t_tpc_id=0;t_tpc_id<fNtpcs;t_tpc_id++){
    geo::TPCID tpcid{{0}, t_tpc_id};
    DrawRectangle(Rectangle_tpc[0], fGeom->TPC(tpcid).BoundingBox().MinX(), fGeom->TPC(tpcid).BoundingBox().MaxX(), 
      fGeom->TPC(tpcid).BoundingBox().MinY(), fGeom->TPC(tpcid).BoundingBox().MaxY()); 
    DrawRectangle(Rectangle_tpc[1], fGeom->TPC(tpcid).BoundingBox().MinX(), fGeom->TPC(tpcid).BoundingBox().MaxX(), 
      fGeom->TPC(tpcid).BoundingBox().MinZ(), fGeom->TPC(tpcid).BoundingBox().MaxZ());
    DrawRectangle(Rectangle_tpc[2], fGeom->TPC(tpcid).BoundingBox().MinY(), fGeom->TPC(tpcid).BoundingBox().MaxY(), 
      fGeom->TPC(tpcid).BoundingBox().MinZ(), fGeom->TPC(tpcid).BoundingBox().MaxZ()); 
    
    for(int i=0; i<3; i++){
      Rectangle_tpc[i]->SetLineWidth(1);
      Rectangle_tpc[i]->SetLineColor(29);
      canvas2D_yz->cd(i+1);
      Rectangle_tpc[i]->Draw();
     }
  }

  //Drawing circle around the muon decay
  for(int i=0; i<3; i++){

      int idx = int(i/2);
      for (auto const& spp : l_trackStartEnd) {
        Tracks_lines[i] = new TPolyLine(2);
        Tracks_lines[i]->SetPoint(0, spp[idx],spp[int(i+1-idx)]);
        Tracks_lines[i]->SetPoint(1, spp[idx+3],spp[int(i+1-idx)+3]);
        if(spp[6]>0.5){
          Tracks_lines[i]->SetLineColor(4);
          Tracks_lines[i]->SetLineWidth(2);
        }
        else{
          Tracks_lines[i]->SetLineColor(15);
          Tracks_lines[i]->SetLineWidth(1);
        }
        canvas2D_yz->cd(i+1);
        Tracks_lines[i]->Draw();
    }

      for (auto const& sp : l_mich_centers) {
        TArc *arc = new TArc(sp[idx],sp[int(i+1-idx)],fMichelRadiusSphere);
        arc->SetLineColor(kRed);
        arc->SetLineWidth(2);
        arc->SetFillStyle(0);
        canvas2D_yz->cd(i+1);
        arc->Draw();
    }
  }

  canvas2D_yz->cd(4);
  axos_2D[3]->Draw();

  //Adding a legend associated with mich muon track candidate
  TH1 *h_leg_michmuon = new TH1F("", "", 1, 0, 1);
  h_leg_michmuon->SetLineColor(4);
  h_leg_michmuon->SetLineWidth(2);
  
  //Adding a legend associated with track that are not mich muon track candidate
  TH1 *h_leg_muon = new TH1F("", "", 1, 0, 1);
  h_leg_muon->SetLineColor(15);
  h_leg_muon->SetLineWidth(1);

  //Adding a legend associated with mich muon centers
  TH1 *h_leg_michcenter = new TH1F("", "", 1, 0, 1);
  h_leg_michcenter->SetLineColor(kRed);
  h_leg_michcenter->SetLineWidth(1);

  //Adding a legend associated with electrons from simulation coming from a muon
  TH1F *h_leg_sp_mc_elec = new TH1F("", "", 1, 0, 1);
  h_leg_sp_mc_elec->SetMarkerStyle(8);
  h_leg_sp_mc_elec->SetMarkerSize(0.5);    
  h_leg_sp_mc_elec->SetMarkerColor(6);

  //Adding a legend associated with space points that are electron michel candidate
  TH1F *h_leg_sp_michmuon = new TH1F("", "", 1, 0, 1);
  h_leg_sp_michmuon->SetMarkerStyle(4);
  h_leg_sp_michmuon->SetMarkerSize(2);    
  h_leg_sp_michmuon->SetMarkerColor(3);

  //Adding a legend associated with space points that are associated with conduction hit
  TH1F *h_leg_sp_allhits = new TH1F("", "", 1, 0, 1);
  h_leg_sp_allhits->SetMarkerStyle(8);
  h_leg_sp_allhits->SetMarkerSize(0.5);    
  h_leg_sp_allhits->SetMarkerColor(45);

  auto legend = new TLegend(0.1,0.2,0.9,0.9);
 //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
  legend->AddEntry(h_leg_sp_allhits,"Reco SpacePoint of all hits","p");
  legend->AddEntry(h_leg_michcenter,"Michel candidate centers","l");
  legend->AddEntry(h_leg_michmuon,"Reco tracks identified as decaying muons (cheat)","l");
  legend->AddEntry(h_leg_muon,"Reco tracks not identified as decaying muons","l");
  legend->AddEntry(h_leg_sp_mc_elec,"True position of michel electrons (cheat)","p");
  legend->AddEntry(h_leg_sp_michmuon,"Reco SP of michel electron candidates","p");
  legend->Draw();

  return;
}
void pdvdana::PDVDCheckMichel::Drawing2D_yz_AddingSpacePoints(TCanvas* canvas_2D, vector<vector<double>> l_SpacePoints, int l_size=1, int l_color=2, int l_style=8){

  TPolyMarker* l_SpacePoints_v = new TPolyMarker(l_SpacePoints.size());
  int actual_point=0;
  for (auto const& sp : l_SpacePoints) {
    l_SpacePoints_v->SetPoint(actual_point, sp[1], sp[2]);
    actual_point++;
  }
  
  canvas_2D->cd();
  SetHitStyle(l_SpacePoints_v, l_size, l_color, l_style);      
  l_SpacePoints_v->Draw();
  
  return;
}


//Drawing 2D tracks
void pdvdana::PDVDCheckMichel::Drawing2D_yz_HitsAndTracks(TCanvas* canvas2D_yz, vector<vector<double>> l_trackStartEnd, vector<vector<double>> l_mich_centers){


  TPolyLine* Rectangle = new TPolyLine(5);
  TPolyLine* Rectangle_tpc = new TPolyLine(5);

  DrawRectangle(Rectangle, fgeoYmin, fgeoYmax, fgeoZmin, fgeoZmax);

  TH2F* axos_2D = new TH2F("axes_1", "Y(Z); Y; Z", 1, -400, 400, 1, -10, 400);
  canvas2D_yz->cd();
  axos_2D->Draw();
  Rectangle->Draw();

  for(unsigned t_tpc_id=0;t_tpc_id<fNtpcs;t_tpc_id++){
    geo::TPCID tpcid{{0}, t_tpc_id};
    DrawRectangle(Rectangle_tpc, fGeom->TPC(tpcid).BoundingBox().MinY(), fGeom->TPC(tpcid).BoundingBox().MaxY(), fGeom->TPC(tpcid).BoundingBox().MinZ(), fGeom->TPC(tpcid).BoundingBox().MaxZ()); 
    Rectangle_tpc->SetLineWidth(1);
    Rectangle_tpc->SetLineColor(29);
    canvas2D_yz->cd();
    Rectangle_tpc->Draw();   
  }

  // Drawing the tracks from start and end
  for (auto const& spp : l_trackStartEnd) {
    TPolyLine* Tracks_lines = new TPolyLine(2);
    Tracks_lines->SetPoint(0, spp[1],spp[2]);
    Tracks_lines->SetPoint(1, spp[4],spp[5]);
    if(spp[6]>0.5){
      Tracks_lines->SetLineColor(4);
      Tracks_lines->SetLineWidth(2);
    }
    else{
      Tracks_lines->SetLineColor(15);
      Tracks_lines->SetLineWidth(1);
    }
    canvas2D_yz->cd();
    Tracks_lines->Draw();
  }
  //Drawing circle around the muon decay
  for (auto const& sp : l_mich_centers) {
    TArc *arc = new TArc(sp[1],sp[2],fMichelRadiusSphere);
    arc->SetLineColor(kRed);
    arc->SetLineWidth(2);
    arc->SetFillStyle(0);
    canvas2D_yz->cd();
    arc->Draw();
  }
  
  canvas2D_yz->cd();
  //Adding a legend associated with mich muon track candidate
  TH1 *h_leg_michmuon = new TH1F("", "", 1, 0, 1);
  h_leg_michmuon->SetLineColor(4);
  h_leg_michmuon->SetLineWidth(2);
  
  //Adding a legend associated with track that are not mich muon track candidate
  TH1 *h_leg_muon = new TH1F("", "", 1, 0, 1);
  h_leg_muon->SetLineColor(15);
  h_leg_muon->SetLineWidth(1);

  //Adding a legend associated with mich muon centers
  TH1 *h_leg_michcenter = new TH1F("", "", 1, 0, 1);
  h_leg_michcenter->SetLineColor(kRed);
  h_leg_michcenter->SetLineWidth(1);

  //Adding a legend associated with electrons from simulation coming from a muon
  TH1F *h_leg_sp_mc_elec = new TH1F("", "", 1, 0, 1);
  h_leg_sp_mc_elec->SetMarkerStyle(8);
  h_leg_sp_mc_elec->SetMarkerSize(0.5);    
  h_leg_sp_mc_elec->SetMarkerColor(6);

  //Adding a legend associated with space points that are electron michel candidate
  TH1F *h_leg_sp_michmuon = new TH1F("", "", 1, 0, 1);
  h_leg_sp_michmuon->SetMarkerStyle(8);
  h_leg_sp_michmuon->SetMarkerSize(1);    
  h_leg_sp_michmuon->SetMarkerColor(3);

  //Adding a legend associated with space points that are associated with conduction hit
  TH1F *h_leg_sp_allhits = new TH1F("", "", 1, 0, 1);
  h_leg_sp_allhits->SetMarkerStyle(8);
  h_leg_sp_allhits->SetMarkerSize(0.5);    
  h_leg_sp_allhits->SetMarkerColor(45);

  auto legend = new TLegend(0.1,0.75,0.4,0.9);
  //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
  legend->AddEntry(h_leg_sp_allhits,"Reco SpacePoint of all hits","p");
  legend->AddEntry(h_leg_michcenter,"Michel candidate centers","l");
  legend->AddEntry(h_leg_michmuon,"Reco tracks identified as decaying muons (cheat) ","l");
  legend->AddEntry(h_leg_muon,"Reco tracks not identified as decaying muons","l");
  legend->AddEntry(h_leg_sp_mc_elec,"True position of michel electrons (cheat)","p");
  legend->AddEntry(h_leg_sp_michmuon,"Reco SP of michel electron candidates","p");
  legend->Draw();
  return;
}

pair<int,int> pdvdana::PDVDCheckMichel::GetWireFromGeoPoint(geo::Point_t localpoint){
  
  geo::TPCID tpc_s = fGeom->FindTPCAtPosition(localpoint);
  int tpc_id = tpc_s.TPC;
  int wire_id = -9999;

  if(! (tpc_s.isValid)){
    if (fLogLevel >= fFlagWarning){
      cout << myname << "---!!--- TPCs are not valid : " << endl;
      cout << myname << "---!!--- TPCs ID: " << tpc_s.TPC << endl;
    }
  }
  else{
    geo::WireID g_wireID_s;
    try{
      g_wireID_s = fGeom->NearestWireID(localpoint, geo::PlaneID(0, tpc_s.TPC, 2));
    }
    catch(geo::InvalidWireError const& e_1) {
      g_wireID_s = e_1.suggestedWireID(); // pick the closest valid wire
    }
    wire_id = g_wireID_s.Wire;
  }
/*  if (wire_id!=-9999){
      wire_id += GetWireOffset(2, tpc_s.TPC);
  }*/
  pair<int,int> ReturnWireTPC (tpc_id, wire_id);
  return ReturnWireTPC;

}


bool pdvdana::PDVDCheckMichel::IsMichelMichel(simb::MCParticle l_particle, simb::MCParticle l_mother_particle){

  TString ParticleCreationProcess = (TString) l_particle.Process();
  
  if( (abs(l_particle.PdgCode())==11 and (l_particle.Mother()>0)) //electron not being a primal particle
    and (ParticleCreationProcess.Contains("Decay") and fabs(l_mother_particle.PdgCode())==13)) //the mother particle is a decaying muon
      return true;

  return false;

}
bool pdvdana::PDVDCheckMichel::IsMuonMichel(simb::MCParticle l_particle){

  bool t_inside = false;
  bool t_ismuon = false;

  TString ParticleDeadProcess = (TString) l_particle.EndProcess();
  TVector3 pos_init = l_particle.Position().Vect();
  TVector3 pos_final = l_particle.EndPosition().Vect();

  if(IsThisTrackEnteringTPC(pos_init, pos_final) and (IsThisPointInsideTPC(pos_final) or IsThisPointInsideTPC(pos_init)))
    t_inside = true;

  if(fabs(l_particle.PdgCode()) == 13 and (ParticleDeadProcess.Contains("Decay"))){
    if((l_particle.NumberDaughters())>0)
      t_ismuon=true;   
    
  }
  
  
  return (t_ismuon and t_inside);
}

void pdvdana::PDVDCheckMichel::DrawRectangle(TPolyLine *Rectangle, float a_1, float a_2, float b_1, float b_2){    
  Rectangle->SetPoint(0,a_1, b_1);
  Rectangle->SetPoint(1,a_1, b_2);    
  Rectangle->SetPoint(2,a_2, b_2);
  Rectangle->SetPoint(3,a_2, b_1);
  Rectangle->SetPoint(4,a_1, b_1);
  Rectangle->SetLineColor(12);
  Rectangle->SetLineWidth(3);
  return;
}
// check the characteristics of a track -> length in x, y, z, total length and theta/phi angles
void pdvdana::PDVDCheckMichel::checkTrackCharacs( TVector3 l_track, vector<float> &TrackCharac ){
  
  double norm = l_track.Mag();
  float Dx = l_track.X();
  float Dy = l_track.Y();
  float Dz = l_track.Z();

  if (l_track.X()>0){
    Dx = -Dx;
    Dz = -Dz;
  }

  TrackCharac[0] = Dx;
  TrackCharac[1] = Dy;
  TrackCharac[2] = Dz;
  TrackCharac[3] = norm;
  TrackCharac[4] = 0;
  TrackCharac[5] = 0;

  if (norm<1)
    return;

  // Determination of the angle
  double theta = 180 - acos(abs(Dx)/norm)*180/3.1415;
  double phi = 0;

  if(Dy >= 0){
    if (Dz>=0)
      phi = atan(abs(Dz)/abs(Dy))*180/3.1415 ;
    else{
      phi = 270+atan(abs(Dy)/abs(Dz))*180/3.1415;
    }
  }
  else{
    if (Dz<0){
      phi= atan(abs(Dz)/abs(Dy))*180/3.1415+180;
    }
    else{
      phi= 90+atan(abs(Dy)/abs(Dz))*180/3.1415;
    }
  }
  TrackCharac[4] = theta;
  TrackCharac[5] = phi;
  return;
}


bool pdvdana::PDVDCheckMichel::IsThisTrackEnteringTPC(TVector3 PointA, TVector3 PointB){

  //From A and B, check if a point in the segment AB is inside the TPC
  double alpha = (PointB.Y()*PointA.Z()-PointA.Y()*PointB.Z())/(PointA.X()*PointB.Y()-PointB.X()*PointA.Y());
  double beta = (PointA.Z()*PointB.X()-PointB.Z()*PointA.X())/(PointA.Y()*PointB.X()-PointB.Y()*PointA.X());
  float var_z=0;
  for (float var_x=fgeoXmin; var_x<fgeoXmax; var_x=var_x+1){
    for (float var_y=fgeoXmin; var_y<fgeoXmax; var_y=var_y+1){
      var_z=alpha*var_x+beta*var_y;
      if((var_z<fgeoZmin) && (var_z>fgeoZmin)){
        return false;
      }
    }
  }
  return true;
}

//From A check if the point is inside the TPC
bool pdvdana::PDVDCheckMichel::IsThisPointInsideTPC(TVector3 PointA){
  
  bool l_isin = false;
  if((PointA.Y() > (fgeoYmin+fMichelRadiusSphere) && PointA.Y()<(fgeoYmax-fMichelRadiusSphere)) &&
  (PointA.Z() > (fgeoZmin+fMichelRadiusSphere) && PointA.Z()<(fgeoZmax-fMichelRadiusSphere))){
    l_isin = true;
  }

  return l_isin;
}

unsigned pdvdana::PDVDCheckMichel::GetWireOffset(unsigned plane_id, unsigned tpc_id){

  if (plane_id==0)
    return fOffsetWireID_u[tpc_id];
  else if (plane_id==1)
    return fOffsetWireID_v[tpc_id];
  else if(plane_id==2)
    return fOffsetWireID_z[tpc_id];

  return 0;

}


DEFINE_ART_MODULE(pdvdana::PDVDCheckMichel)
