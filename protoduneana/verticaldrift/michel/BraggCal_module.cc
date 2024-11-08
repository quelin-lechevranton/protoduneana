/*
 * Class: BraggCal
 * Plugin Type: analyzer
 * Art version: 3.14.04
 * LArSoft version: v09_91_04d00
 * File: BraggCal_module.cc
 * 
 * Analyzing Calorimetic information from tracks
 * to identify decaying muons through Bragg peak
 * 
 * Author: Jérémy Quélin Lechevranton
*/


// Art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
// #include "larsim/MCCheater/BackTrackerService.h"
// #include "larsim/MCCheater/PhotonBackTrackerService.h"
// #include "larsim/MCCheater/ParticleInventoryService.h"
// #include "larsim/MCCheater/BackTracker.h"
// #include "larsim/MCCheater/BackTrackerService.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"

// DUNE includes
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaTrackUtils.h"

// ROOT includes
#include <TTree.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH2F.h>

// std includes
#include <vector>
#include <iterator> 
#include <string>

// using namespace std;
using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::vector;
using std::string;
using std::pair;

namespace ana { class BraggCal; }

class ana::BraggCal : public art::EDAnalyzer {
public:

    explicit BraggCal(fhicl::ParameterSet const& fcl);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    BraggCal(BraggCal const &) = delete;
    BraggCal(BraggCal &&) = delete;
    BraggCal & operator = (BraggCal const &) = delete;
    BraggCal & operator = (BraggCal &&) = delete;
    
    void analyze(art::Event const & evt) override; 
    void beginJob() override;
    void endJob()   override;

private:

    bool v; //verbosity
    vector<vector<string>> vvsProducts;

    TTree* tTree;
    TH2F* hdEdx0;
    TH2F* hdEdx1;
    vector<float> vfdEdx;
    vector<float> vfResRange;
    vector<float> vfdEdx1;
    vector<float> vfResRange1;

    art::InputTag   tag_pfp,
                    tag_clu,
                    tag_trk,
                    tag_cal,
                    tag_shw,
                    tag_spt,
                    tag_hit;

    float fBraggLength;
    size_t iCalMin;
    size_t iThreshold;
    size_t iNormalize;
    float fdEdx_MIP;
    float fdEdx_min_ratio;
    float fAvgdEdxMin;
    float fAvgdEdxMax;
    float fBraggMin;
    
};

ana::BraggCal::BraggCal(fhicl::ParameterSet const & fcl) :
    EDAnalyzer{fcl},
    v(fcl.get<bool>("Verbose")),
    vvsProducts(fcl.get<vector<vector<string>>>("Products")),
    fBraggLength(fcl.get<float>("BraggLength")),
    iCalMin(fcl.get<size_t>("CalMin")),
    iThreshold(fcl.get<size_t>("Threshold")),
    iNormalize(fcl.get<size_t>("Normalize")),
    fdEdx_MIP(fcl.get<float>("dEdx_MIP")),
    fdEdx_min_ratio(fcl.get<float>("dEdx_min_ratio")),
    fAvgdEdxMin(fcl.get<float>("AvgdEdxMin")),
    fAvgdEdxMax(fcl.get<float>("AvgdEdxMax")),
    fBraggMin(fcl.get<float>("BraggMin"))
{


    art::ServiceHandle<art::TFileService> tfs;

    tTree = tfs->make<TTree>("Bragg","Bragg");
    tTree->Branch("vfdEdx",&vfdEdx);
    tTree->Branch("vfResRange",&vfResRange);
    tTree->Branch("vfdEdx1",&vfdEdx1);
    tTree->Branch("vfResRange1",&vfResRange1);

    hdEdx0 = new TH2F(
        "hdEdx0",
        "before bragg selection;residual range (cm);dE/dx (MeV/cm)",
        100,0,300,
        50,0,5
    );
    hdEdx1 = new TH2F(
        "hdEdx1",
        "after bragg selection;residual range (cm);dE/dx (MeV/cm)",
        100,0,300,
        50,0,5
    );

    cout << hdEdx0->GetEntries() << endl;
    cout << hdEdx1->GetEntries() << endl;

    for (vector<string> prod : vvsProducts) {

        string  label=prod[0],
                instance=prod[1],
                object=prod[2],
                process=prod[3];

        if (object == "recob::PFParticle")          tag_pfp= art::InputTag(label,instance);
        else if (object == "recob::Cluster")        tag_clu= art::InputTag(label,instance);
        else if (object == "recob::Track")          tag_trk= art::InputTag(label,instance);
        else if (object == "anab::Calorimetry")     tag_cal= art::InputTag(label,instance);
        else if (object == "recob::Shower")         tag_shw= art::InputTag(label,instance);
        else if (object == "recob::SpacePoint")     tag_spt= art::InputTag(label,instance);
        else if (object == "recob::Hit")            tag_hit= art::InputTag(label,instance);
    } // end prod loop
} // end constructor()

void ana::BraggCal::beginJob() {}

void ana::BraggCal::analyze(const art::Event & evt) {

    if (v) cout << "evt#" << evt.id().event()+1 << "\r" << flush;

    auto const vh_pfp = evt.getValidHandle<vector<recob::PFParticle>>(tag_pfp);
    art::FindManyP<recob::Track>        const fmp_pfp_trk(vh_pfp,evt,tag_trk);

    auto const vh_trk = evt.getValidHandle<vector<recob::Track>>(tag_trk);
    art::FindManyP<anab::Calorimetry>   const fmp_trk_cal(vh_trk,evt,tag_cal);

    vector<art::Ptr<recob::PFParticle>> vp_pfp;
    art::fill_ptr_vector(vp_pfp,vh_pfp);
    
    for (art::Ptr<recob::PFParticle> const & p_pfp : vp_pfp) {

        vector<art::Ptr<recob::Track>> const vp_trk = fmp_pfp_trk.at(p_pfp.key());

        if (vp_trk.empty()) continue;

        art::Ptr<recob::Track> const p_trk = vp_trk[0];
        if (v) cout << "\ttrk#" << p_trk->ID()+1 << "\r" << flush; 

        // is inside ?
        // geometry for fcl ?

        float X_first = p_trk->LocationAtPoint(p_trk->FirstValidPoint()).X();
        float X_last = p_trk->LocationAtPoint(p_trk->LastValidPoint()).X();
        bool upright = X_first > X_last;

        vector<art::Ptr<anab::Calorimetry>> const vp_cal = fmp_trk_cal.at(p_trk.key());

        int i_vp_cal = -1;
        while (vp_cal[++i_vp_cal]->PlaneID().Plane != geo::kW);
        art::Ptr<anab::Calorimetry> const & p_cal = vp_cal[i_vp_cal];

        size_t n_cal = p_cal->ResidualRange().size();

        if (n_cal < iCalMin) continue;
        if (v) cout << "\t\tn cal: " << n_cal << endl;

        float avg_dEdx = 0;
        for (size_t i_cal=0; i_cal < n_cal; i_cal++) {
            float res_range;
            float dEdx = p_cal->dEdx()[i_cal];
            if (upright) {
                res_range = p_cal->ResidualRange()[i_cal];
            } else {
                res_range = p_cal->Range() - p_cal->ResidualRange()[i_cal];
            }
            avg_dEdx += dEdx;

            vfdEdx.push_back(dEdx);
            vfResRange.push_back(res_range);
            hdEdx0->Fill(res_range,dEdx);
        }
        avg_dEdx /= n_cal;

        if (fAvgdEdxMin > avg_dEdx || avg_dEdx > fAvgdEdxMax) continue;
        if (v) cout << "\t\tavg dEdx: " << avg_dEdx << " MeV/cm" << endl;

        float dEdx_min = 0;
        switch (iThreshold) {
            case 1: dEdx_min = fdEdx_MIP * fdEdx_min_ratio; break;
            case 2: dEdx_min = avg_dEdx * fdEdx_min_ratio; break;
            default: ;
        }

        float bragg_int=0;
        size_t n_bragg_int=0;
        size_t n_cal_bragg=0;
        if (upright) {
            while (
                p_cal->ResidualRange()[n_cal_bragg++] < fBraggLength 
                && n_cal_bragg < n_cal
            );

            for (size_t i_cal=0; i_cal < n_cal_bragg; i_cal++) {
                float dEdx = p_cal->dEdx()[i_cal]; 
                if (dEdx < dEdx_min) continue;
                bragg_int += dEdx;
                n_bragg_int ++;
            }
        } else {
            while (
                p_cal->ResidualRange()[n_cal-1-n_cal_bragg++] > p_cal->Range() - fBraggLength
                && n_cal_bragg < n_cal-1
            );

            for (size_t i_cal=n_cal-n_cal_bragg; i_cal < n_cal; i_cal++) {
                float dEdx = p_cal->dEdx()[i_cal]; 
                if (dEdx < dEdx_min) continue;
                bragg_int += dEdx;
                n_bragg_int ++;
            }
        }

        // divide by n_bragg_int ??
        switch (iNormalize) {
            case 0: bragg_int /= fdEdx_MIP;
            case 1: bragg_int /= avg_dEdx;
        }

        if (bragg_int < fBraggMin) continue;
        if (v) cout << "\t\tbragg int: " << bragg_int << " MIP" << endl;
        if (!v) cout << " \u2192 decaying muon at evt#" << evt.id().event() << " trk#" << p_trk->ID() << endl;

        for (size_t i_cal=0; i_cal < n_cal; i_cal++) {
            float res_range;
            float dEdx = p_cal->dEdx()[i_cal];
            if (upright) {
                res_range = p_cal->ResidualRange()[i_cal];
            } else {
                res_range = p_cal->Range() - p_cal->ResidualRange()[i_cal];
            }

            vfdEdx1.push_back(dEdx);
            vfResRange1.push_back(res_range);
            hdEdx1->Fill(res_range,dEdx);
        }
    } // end pfp loop

    tTree->Fill();
} // end analyze()

void ana::BraggCal::endJob() {

    art::ServiceHandle<art::TFileService> tfs;
    TCanvas* cdEdx = tfs->make<TCanvas>("cdEdx","cdEdx");

    cdEdx->cd();
    cdEdx->Divide(2,1);

    cdEdx->cd(1);
    gPad->SetLogz();
    hdEdx0->SetMinimum(1);
    hdEdx0->Draw("colz");

    cdEdx->cd(2);
    gPad->SetLogz();
    hdEdx1->SetMinimum(1);
    hdEdx1->Draw("colz");

    cdEdx->Write("cdEdx");

}

DEFINE_ART_MODULE(ana::BraggCal)
