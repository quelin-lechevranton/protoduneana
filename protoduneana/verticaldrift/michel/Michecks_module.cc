////////////////////////////////////////////////////////////////////////
// Class:       Michecks
// Plugin Type: analyzer
// File:        Michecks_module.cc
//
// Generated at Wed Nov  6 10:31:28 2024 by Jeremy Quelin Lechevranton
////////////////////////////////////////////////////////////////////////

// Art includes
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
#include "art_root_io/TFileService.h"

// LArSoft includes
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "larcore/Geometry/Geometry.h"

#include "larcorealg/Geometry/Exceptions.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/ArtDataHelper/TrackUtils.h"

// ProtoDUNE includes
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"

// ROOT includes
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph2D.h"

// std includes
#include <iostream>
#include <vector>
#include <string>
#include <iterator>

using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::vector;
using std::string;
using std::pair;


namespace ana {
  class Michecks;
}


class ana::Michecks : public art::EDAnalyzer {
public:
    explicit Michecks(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    Michecks(Michecks const&) = delete;
    Michecks(Michecks&&) = delete;
    Michecks& operator=(Michecks const&) = delete;
    Michecks& operator=(Michecks&&) = delete;

    // Required functions.
    void analyze(art::Event const& e) override;

    // Selected optional functions.
    void beginJob() override;
    void endJob() override;

private:

    bool v; //verbosity

    TGraph2D* gMuon;
    TGraph2D* gMichel;


    // Products
    vector<vector<string>> vvsProducts;
    art::InputTag   tag_mcp,
                    tag_sed,
                    tag_wir,
                    tag_hit,
                    tag_clu,
                    tag_trk,
                    tag_spt;


    const geo::Geometry* fGeom;
    protoana::ProtoDUNETruthUtils truthUtil;

};


ana::Michecks::Michecks(fhicl::ParameterSet const& p)
    : EDAnalyzer{p},
    v(p.get<bool>("Verbosity",false)),
    vvsProducts(p.get<vector<vector<string>>>("Products"))
{
    fGeom = &*art::ServiceHandle<geo::Geometry>();

    for (vector<string> prod : vvsProducts) {

        const string    label   =prod[0],
                        instance=prod[1],
                        object  =prod[2],
                        process =prod[3];

        const auto  tag = art::InputTag(label,instance);

        if (object == "simb::MCParticle")           tag_mcp = tag;
        else if (object == "sim::SimEnergyDeposit") tag_sed = tag;
        else if (object == "recob::Hit")            tag_hit = tag;
        else if (object == "recob::Wire")           tag_wir = tag;
        else if (object == "recob::Cluster")        tag_clu = tag;
        else if (object == "recob::Track")          tag_trk = tag;
        else if (object == "recob::SpacePoint")     tag_spt = tag;
    }

    TGraph2D* gMuon = new TGraph2D();
    gMuon->SetTitle("Truly Decaying Muon Track Points");    
    gMuon->SetMarkerStyle(20);
    gMuon->SetMarkerSize(0.5);
    gMuon->SetMarkerColor(kOrange-3);


    TGraph2D* gMichel = new TGraph2D();
    gMichel->SetTitle("True Michel Electron MCParticle Points");
    gMichel->SetMarkerStyle(20);
    gMichel->SetMarkerSize(0.5);
    gMichel->SetMarkerColor(kAzure+7);

}

void ana::Michecks::beginJob()
{
    // Implementation of optional member function here.
}


void ana::Michecks::analyze(art::Event const& e)
{

    if (v) cout << "evt#" << e.id().event() << "\r" << flush;

    // auto const vh_sed = e.getValidHandle<vector<sim::SimEnergyDeposit>>(tag_sed);
    // vector<art::Ptr<sim::SimEnergyDeposit>> vp_sed;
    // art::fill_ptr_vector(vp_sed,vh_sed);

    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;

    // what is this for?
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);

    // art::FindManyP<recob::Hit>          const fmp_sed_hit(vh_sed,e,tag_hit);

    // auto const vh_hit = e.getValidHandle<vector<recob::Hit>>(tag_hit);
    // vector<art::Ptr<recob::Hit>> vp_hit;
    // art::fill_ptr_vector(vp_hit,vh_hit);

    // art::FindManyP<recob::Wire>         const fmp_hit_wir(vh_hit,e,tag_wir);
    // art::FindManyP<recob::Cluster>      const fmp_hit_clu(vh_hit,e,tag_clu);
    // art::FindManyP<recob::Track>        const fmp_hit_trk(vh_hit,e,tag_trk);
    // art::FindManyP<recob::SpacePoint>   const fmp_hit_spt(vh_hit,e,tag_spt);

    auto const vh_trk = e.getValidHandle<vector<recob::Track>>(tag_trk);
    auto const vh_mcp = e.getValidHandle<vector<simb::MCParticle>>(tag_mcp);

    for (recob::Track const & trk : *vh_trk) {

        if (v) cout << "\ttrk#" << trk.ID()+1 << "\r" << flush;

        // Get the MCParticle associated to the track
        simb::MCParticle const * mcp = truthUtil.GetMCParticleFromRecoTrack(clockData, trk, e, tag_trk.label());

        // if (!mcp) continue;

        // bool has_michel = false;
        // bool has_nue = false;
        // bool has_numu = false;

        if (abs(mcp->PdgCode()) != 13) continue;

        if (v) cout << "\t\tMuon found" << endl;

        if (mcp->NumberDaughters() < 1) continue;

        int i_mcp_michel = -1;

        for (int i_dau=0; i_dau < mcp->NumberDaughters(); i_dau++) {

            // IS IT THE SAME THING ???
            //simb::MCParticle const * mcp_dau = pi_serv->TrackIdToParticle_P(mcp->Daughter(i_dau));
            simb::MCParticle const& mcp_dau = vh_mcp->at(mcp->Daughter(i_dau));

            // if (!mcp_dau) continue;

            if (abs(mcp_dau.PdgCode()) == 11) continue;
            if (mcp_dau.Process() != "Decay") continue;
            
            if (v) cout << "\t\tMichel found at dau#" << i_dau << endl;

            i_mcp_michel = mcp->Daughter(i_dau);


            // switch (abs(mcp_dau.PdgCode())) {
            //     case 11:
            //         if (mcp_dau.Process() != "Decay") break;
            //         has_michel = true;
            //         i_mcp_michel = mcp->Daughter(i_dau);
            //         break;
            //     case 12:
            //         has_nue = true; break;
            //     case 14:
            //         has_numu = true; break;
            //     default: break;
            // }
            
            // if (abs(mcp_dau.PdgCode()) == 11) {
            //     // ((TString) mcp_dau.Process).Contains("Decay")
            //     // mcp_dau.Process().find("Decay") != string::npos
            //     if (mcp_dau.Process() != "Decay") continue;
            //     has_michel = true;
            //     i_mcp_michel = mcp->Daughter(i_dau);
            // }
            // else if (abs(mcp_dau.PdgCode()) == 12) has_nue = true;
            // else if (abs(mcp_dau.PdgCode()) == 14) has_numu = true;
        }

        if (i_mcp_michel < 0) continue;

        simb::MCParticle const& mcp_michel = vh_mcp->at(i_mcp_michel);


        for (size_t i_ppt=0; i_ppt < mcp_michel.NumberTrajectoryPoints(); i_ppt++) {
            gMichel->AddPoint(
                mcp_michel.Vx(i_ppt),
                mcp_michel.Vy(i_ppt),
                mcp_michel.Vz(i_ppt)
            );
        } // end loop over michel points

        for (size_t i_tpt= trk.FirstPoint(); i_tpt < trk.LastPoint(); i_tpt++) {
            if (!trk.HasValidPoint(i_tpt)) continue;
            
            gMuon->AddPoint(
                trk.LocationAtPoint(i_tpt).X(),
                trk.LocationAtPoint(i_tpt).Y(),
                trk.LocationAtPoint(i_tpt).Z()
            );
        } // end loop over track points
    } // end loop over tracks
} // end analyze

void ana::Michecks::endJob()
{
    art::ServiceHandle<art::TFileService> tfs;
    TCanvas* c1 = tfs->make<TCanvas>("c1","c1");

    c1->cd();
    gMuon->Draw("p");
    gMichel->Draw("same p");
    c1->Write();
}

// pair<int,int> ana::Michecks::GetWireTPC(geo::Point_t *pt) {
//     geo::TPCID tpc = fGeom->FindTPCAtPosition(*pt);
//     geo::WireID wire = fGeom->NearestWireID(*pt, tpc);
//     return {tpc.TPC, wire.Wire};
// }

DEFINE_ART_MODULE(ana::Michecks)
