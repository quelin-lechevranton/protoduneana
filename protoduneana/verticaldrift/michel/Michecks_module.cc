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
#include "TLorentzVector.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TEllipse.h"

// std includes
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <iterator>

using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::vector;
using std::string;
using std::stringstream;
using std::pair;
using std::min;
using std::max;
using std::stoi;


namespace ana {
  class Michecks;
  struct Binning {
        int n;
        double min, max;
  };
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

    // Utilities
    art::ServiceHandle<art::TFileService> tfs;

    const geo::Geometry* fGeom;
    const detinfo::DetectorPropertiesService* fDetProp;
    const detinfo::DetectorClocksService* fDetClocks;

    protoana::ProtoDUNETruthUtils truthUtil;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    
    //verbosity
    bool v; 

    vector<vector<Binning>> bin_ch;
    ana::Binning bin_tick;

    TGraph2D* gMuon;
    TGraph2D* gMichel;
    // TH1F* hWS;
    // TH1F* hSignal;
    // TH1F* hWire;
    // TH2F* hWT;

    vector<vector<TH2F*>> hHT; // [plane][section]
    vector<vector<TGraph*>> gMichelHits; // [plane][section]
    // vector<vector<TH1F*>> hTrueE; // [plane][section]
    TH1F* hTrueE;
    TH1F* hRecoE;



    // Products
    vector<vector<string>> vvsProducts;
    art::InputTag   tag_mcp,
                    tag_sed,
                    tag_wir,
                    tag_hit,
                    tag_clu,
                    tag_trk,
                    tag_spt;


    float fMichelTimeRadius, // in µs
          fMichelSpaceRadius, // in cm
          fTrackLengthCut; // in cm

    bool IsInVolume(double x, double y, double z, double eps);
    bool IsInVolume(float x, float y, float z, float eps);
    bool IsInVolume(TLorentzVector const& V);
    bool IsInVolume(geo::Point_t const& P);
    bool IsUpright(recob::Track const& T);
    int ChannelSection(unsigned int ch);
};


ana::Michecks::Michecks(fhicl::ParameterSet const& p)
    : EDAnalyzer{p},
    v(p.get<bool>("Verbosity",false)),
    vvsProducts(p.get<vector<vector<string>>>("Products")),

    fMichelTimeRadius(p.get<float>("MichelTimeRadius")), //in µs
    fMichelSpaceRadius(p.get<float>("MichelSpaceRadius")), //in cm

    fTrackLengthCut(p.get<float>("LengthCut")) // in cm
{

    // Basic Utilities
    fGeom = &*art::ServiceHandle<geo::Geometry>();
    fDetProp = &*art::ServiceHandle<detinfo::DetectorPropertiesService>();    
    fDetClocks = &*art::ServiceHandle<detinfo::DetectorClocksService>();

    // Retrieving product tags
    for (vector<string> prod : vvsProducts) {

        const string    process     =prod[0],
                        label       =prod[1],
                        instance    =prod[2],
                        type        =prod[3];

        const auto  tag = art::InputTag(label,instance);

        if      (type == "simb::MCParticle")        tag_mcp = tag;
        else if (type == "sim::SimEnergyDeposit")   tag_sed = tag;
        else if (type == "recob::Hit")              tag_hit = tag;
        else if (type == "recob::Wire")             tag_wir = tag;
        else if (type == "recob::Cluster")          tag_clu = tag;
        else if (type == "recob::Track")            tag_trk = tag;
        else if (type == "recob::SpacePoint")       tag_spt = tag;
    }

    int section_n = fGeom->Nchannels()/12;

    // Binning of the Frame Channel-Tick for which plane and section
    bin_ch = vector<vector<Binning>>(3, vector<Binning>(4));
    for (int plane=0; plane<3; plane++) {
        for (int section=0; section<4; section++) {
            bin_ch[plane][section].n    = section_n;
            bin_ch[plane][section].min  = (3*section + plane) * section_n;
            bin_ch[plane][section].max  = (3*section + plane + 1) * section_n - 1;
        }
    }
    bin_tick.n = 1500;                  bin_tick.min = 0.;  bin_tick.max = 6000.;


    gMuon = new TGraph2D();
    gMuon->SetTitle("Truly Decaying Muon Track Points");    
    gMuon->SetMarkerStyle(20);
    gMuon->SetMarkerSize(0.5);
    gMuon->SetMarkerColor(kOrange-3);


    gMichel = new TGraph2D();
    gMichel->SetTitle("True Michel Electron MCParticle Points");
    gMichel->SetMarkerStyle(20);
    gMichel->SetMarkerSize(0.5);
    gMichel->SetMarkerColor(kAzure+7);

    // hWS = new TH1F("hWS","#Wire per channel",4,0,4);
    // hSignal = new TH1F("hSignal","Signal",60,0,60);


    // Initializing TH2F Channel-Tick for each plane and section
    hHT = vector<vector<TH2F*>>(3, vector<TH2F*>(4));
    for (int plane=0; plane<3; plane++) {
        for (int section=0; section<4; section++) {
            stringstream name; name << "hHT" << plane << section;
            hHT[plane][section] = new TH2F(
                name.str().c_str(),"",
                bin_tick.n, bin_tick.min, bin_tick.max,
                bin_ch[plane][section].n, bin_ch[plane][section].min, bin_ch[plane][section].max
            );
        }
    }

    // Initializing TGraph Channel-Tick for each plane and section
    gMichelHits = vector<vector<TGraph*>>(3, vector<TGraph*>(4));
    for (int plane=0; plane<3; plane++) {
        for (int section=0; section<4; section++) {
            gMichelHits[plane][section] = new TGraph();
        }
    }

    // Initializing TH1F True Energy for each plane and section
    // hTrueE = vector<vector<TH1F*>>(3, vector<TH1F*>(4));
    // for (int plane=0; plane<3; plane++) {
    //     for (int section=0; section<4; section++) {
    //         stringstream name; name << "hE" << plane << section;
    //         hTrueE[plane][section] = new TH1F(
    //             name.str().c_str(),"",
    //             100, 0, 100
    //         );
    //     }
    // }
    hTrueE = new TH1F("hTrueE","True Michel Electron Energy",100,0,100);
    hRecoE = new TH1F("hRecoE","Reconstructed Michel Electron Energy",100,0,100);
}

void ana::Michecks::beginJob()
{
    auto const clockData = fDetClocks->DataForJob();
    auto const detProp = fDetProp->DataForJob(clockData);

    if (v) {
        geo::CryostatID cryoid{0};
        cout << "\033[94m" << "BeginJob: Detector dimensions===================================================" << "\033[0m" << endl
             << "Number of channels: " << fGeom->Nchannels() << endl
             << "Number of ticks: " << "???" << endl
             << "Cryostat coordinates: " << fGeom->Cryostat(cryoid).Min() << " - " << fGeom->Cryostat(cryoid).Max() << endl
             << "\tvolume: " << fGeom->Cryostat(cryoid).Width() << " x " << fGeom->Cryostat(cryoid).Height() << " x " << fGeom->Cryostat(cryoid).Length() << endl;
        geo::TPCID tpc0{cryoid, 0}, tpcN{cryoid, fGeom->NTPC(cryoid)-1};
        cout << "\tactive coordinates: " << fGeom->TPC(tpc0).Min() << " - " << fGeom->TPC(tpcN).Max() << endl
             << "\tactive volume: " << (fGeom->TPC(tpc0).Max() - fGeom->TPC(tpcN).Min()).X() << " x " << (fGeom->TPC(tpc0).Max() - fGeom->TPC(tpcN).Min()).Y() << " x " << (fGeom->TPC(tpc0).Max() - fGeom->TPC(tpcN).Min()).Z() << endl
             << "TPCs:" << endl; 
        for (unsigned int i_tpc=0; i_tpc < fGeom->NTPC(); i_tpc++) {
            geo::TPCID tpcid{cryoid, i_tpc};
            cout << "\tTPC#" << i_tpc << "\tcoordinates: " << fGeom->TPC(tpcid).Min() << " - " << fGeom->TPC(tpcid).Max() << endl
                 << "\t\tvolume: " << fGeom->TPC(tpcid).Width() << " x " << fGeom->TPC(tpcid).Height() << " x " << fGeom->TPC(tpcid).Length() << endl
                 << "\t\tactive volume: " << fGeom->TPC(tpcid).ActiveWidth() << " x " << fGeom->TPC(tpcid).ActiveHeight() << " x " << fGeom->TPC(tpcid).ActiveLength() << endl
                 << "\t\tPlanes:";
            for (unsigned int i_plane=0; i_plane < fGeom->Nplanes(); i_plane++) {
                geo::PlaneID planeid{tpcid, i_plane};
                vector<string> plane_names = {"U","V","W"};
                cout << "\t" << plane_names[i_plane] << " #Wires: " << fGeom->Nwires(planeid);
            } // end loop over Planes
            cout << endl;
        } // end loop over TPCs
        cout << "\033[94m" << "End of detector dimensions======================================================" << "\033[0m" << endl;
    }
}


void ana::Michecks::analyze(art::Event const& e)
{
    // what is this for?
    auto const clockData = fDetClocks->DataFor(e);
    // auto const detProp = fDetProp->DataFor(e,clockData);



    if (v) cout << "evt#" << e.id().event() << "\r" << flush;

    auto const vh_trk = e.getValidHandle<vector<recob::Track>>(tag_trk);
    vector<art::Ptr<recob::Track>> vp_trk;
    art::fill_ptr_vector(vp_trk, vh_trk);

    auto const vh_mcp = e.getValidHandle<vector<simb::MCParticle>>(tag_mcp);

    // art::FindManyP<recob::Hit> const fmp_trk_hit(vh_trk, e, tag_trk);

    // vector<int> decaying_muon_trackids,
                // michel_trackids;    
    
    // struct HitPos {
    //     float tick;
    //     float channel;
    // };

    // vector<vector<HitPos>> michel_hits(3);
    // vector<geo::Point_t> pt_mu_decays;

    for (art::Ptr<recob::Track> const& p_trk : vp_trk) {

        if (v) cout << "\ttrk#" << p_trk->ID()+1 << "\r" << flush;

        if (p_trk->ID() != (int) p_trk.key()) {
            cerr << endl << "\033[91m" << "Track ID and key do not match" << "\033[0m" << endl;
            continue;
        }

        if (v) cout << "\t\ttrk2mcp...";
        simb::MCParticle const * mcp = truthUtil.GetMCParticleFromRecoTrack(clockData, *p_trk, e, tag_trk.label());

        // if (!mcp) continue;

        if (v) cout << "\tis muon...";
        if (abs(mcp->PdgCode()) != 13) continue;

        // if (v) cout << "\t\tMuon found with #dau:" << mcp.NumberDaughters() << endl;

        // if (mcp->NumberDaughters() < 1) continue;

        if (v) cout << "\tloop over dau...";
        for (int i_dau=0; i_dau < mcp->NumberDaughters(); i_dau++) {
            
            simb::MCParticle const * mcp_dau = pi_serv->TrackIdToParticle_P(mcp->Daughter(i_dau));    

            // if (!mcp_dau) continue;

            if (v) cout << "\tis decaying...";
            if (abs(mcp_dau->PdgCode()) != 11) continue;
            if (mcp_dau->Process() != "Decay") continue;


            if (v) cout << endl << "\t\tMichel at dau#" << i_dau << " at (" << mcp_dau->Vx(0) << "," << mcp_dau->Vy(0) << "," << mcp_dau->Vz(0) << ")";
            if (!IsInVolume(mcp_dau->Position(0))) {
                if (v) cout << "\033[91m" << " outside" << "\033[0m" << endl;
                continue;
            }
            if (v) cout << "\033[92m" << " inside" << "\033[0m";

            hTrueE->Fill(mcp_dau->E());
            // pt_mu_decays.emplace_back(mcp_dau->EndPosition().Vect());

            vector<const recob::Hit*> hits_michel = truthUtil.GetMCParticleHits(clockData, *mcp_dau, e, tag_hit.label(), false);

            if (v) cout << "w/ " << hits_michel.size() << " hits" << endl;

            // float average_time = 0;
            // float average_channel = 0;
            // vector<float> average_tick_per_plane(3,0);
            // vector<float> average_channel_per_plane(3,0);
            double RecoE = 0;
            if (v) cout << "\t\tlooping over michel hits...";
            for (recob::Hit const* hit : hits_michel) {

                int plane = hit->View();
                int section = ChannelSection(hit->Channel());
                gMichelHits.at(plane).at(section)->AddPoint(hit->PeakTime(), hit->Channel());

                RecoE += hit->SummedADC(); 

                // average_tick_per_plane.at(plane) += hit->PeakTime();
                // average_channel_per_plane.at(plane) += hit->Channel();
            }
            if (v) cout << " done" << endl;

            hRecoE->Fill(RecoE);

            // for (int plane=0; plane<3; plane++) {

            //     average_tick_per_plane.at(plane) /= hits_michel.size();
            //     average_channel_per_plane.at(plane) /= hits_michel.size();

            //     // michel_hits.at(plane).push_back({average_tick_per_plane.at(plane), average_channel_per_plane.at(plane)});

            //     // gMichelHits.at(plane)->AddPoint(average_tick_per_plane.at(plane), average_channel_per_plane.at(plane));
            // }
            
        } // end loop over daughters

        // vector<art::Ptr<recob::Hit>> vp_hit = fmp_trk_hit.at(p_trk.key());




        // for (size_t i_tpt= trk.FirstPoint(); i_tpt < trk.LastPoint(); i_tpt++) {
        //     if (!trk.HasValidPoint(i_tpt)) continue;

        //     if (v) cout << "\t\t\tAdding muon point#" << i_tpt+1 
        //                 << "/" << trk.NPoints() 
        //                 << ":(" << trk.LocationAtPoint(i_tpt).X() << ","
        //                 << trk.LocationAtPoint(i_tpt).Y() << ","
        //                 << trk.LocationAtPoint(i_tpt).Z() << ")"
        //                 << "\r" << flush;
        //     //             << endl;
            
        //     gMuon->AddPoint(
        //         trk.LocationAtPoint(i_tpt).X(),
        //         trk.LocationAtPoint(i_tpt).Y(),
        //         trk.LocationAtPoint(i_tpt).Z()
        //     );
        // } // end loop over track points

    } // end loop over tracks

    auto const & vh_hit = e.getValidHandle<vector<recob::Hit>>(tag_hit);
    vector<art::Ptr<recob::Hit>> vp_hit;
    art::fill_ptr_vector(vp_hit, vh_hit);

    art::FindManyP<recob::Track> fmp_hit2trk(vh_hit, e, tag_trk);


    if (v) cout << "\tlooping over all hits...";
    for (art::Ptr<recob::Hit> const& p_hit : vp_hit) {

        // geo::WireID wireid = fGeom->ChannelToWire(hit.Channel()).at(0);
        int plane = p_hit->View();
        int section = ChannelSection(p_hit->Channel());

        hHT.at(plane).at(section)->Fill(p_hit->PeakTime(), p_hit->Channel());

        // if (p_hit->View() != 2) continue;

        // vector<art::Ptr<recob::Track>> vp_trk = fmp_hit2trk.at(p_hit.key());

        // bool BelongsToTrack = false;
        // for (art::Ptr<recob::Track> const& p_trk : vp_trk) {
        //     BelongsToTrack = p_trk->Length() < fTrackLengthCut;
        //     if (BelongsToTrack) break;
        // }
        // if (BelongsToTrack) continue;
    } // end loop over hits
    if (v) cout << " done" << endl;
} // end analyze

void ana::Michecks::endJob()
{

    if (v) cout << "\033[94m" << "EndJob: Drawing=================================================================" << "\033[0m" << endl;

    auto const clockData = fDetClocks->DataForJob();
    // auto const detProp = fDetProp->DataForJob(clockData);




    TCanvas* c3d = tfs->make<TCanvas>("c3d","Decay Muon and Michel Electron");

    c3d->cd();
    gMuon->Draw("p");
    gMichel->Draw("same p");
    c3d->Write();

    // TCanvas* c2 = tfs->make<TCanvas>("cwc","Wire content");

    // c2->Divide(2,2);
    // c2->cd(1);
    // gPad->SetLogy();
    // hWS->Draw();
    // c2->cd(2);
    // gPad->SetLogy();
    // hSignal->Draw();
    // c2->cd(3);
    // hWire->Draw();
    // c2->Write();

    // TCanvas* c3 = tfs->make<TCanvas>("cwt","Wire Time TH2F");

    // c3->cd();
    // gPad->SetLogz();
    // hWT->Draw("colz");
    // // hWT->ProjectionX("px",i,i+1)->Draw("same");
    // c3->Write();

    if (v) cout << "sampling rate: " << detinfo::sampling_rate(clockData) << endl;

    TCanvas* cU = tfs->make<TCanvas>("cU","U Plane");
    TCanvas* cV = tfs->make<TCanvas>("cV","V Plane");
    TCanvas* cW = tfs->make<TCanvas>("cW","W Plane");
    vector<TCanvas*> cPlanes = {cU,cV,cW};


    for (int plane=0; plane<3; plane++) {

        if (v) cout << "plane#" << plane << endl;
        cPlanes.at(plane)->Divide(2,2);

        for (int section=0; section<4; section++) {

            if (v) cout << "\tsection#" << section << endl;
            cPlanes.at(plane)->cd(section+1);

            stringstream title;
            title << "Plane " << plane << " Section " << section << ";Tick;Channel";
            gPad->DrawFrame(
                bin_tick.min, bin_ch[plane][section].min, 
                bin_tick.max, bin_ch[plane][section].max, 
                title.str().c_str()
            );
            gPad->SetLogz();

            hHT.at(plane).at(section)->Draw("colz same");

            if (v) cout << "\t\tDrawn Hit/Time" << endl;

            TGraph* g = gMichelHits.at(plane).at(section);

            for (int i=0; i<g->GetN(); i++) {

                float michel_tick = g->GetPointX(i);
                float michel_channel = g->GetPointY(i);
                float tick_radius = fMichelTimeRadius / (detinfo::sampling_rate(clockData) * 1e-3);
                float channel_radius = fMichelSpaceRadius / 0.5; // pitch in cm WHERE TO FIND IT ?

                auto e = new TEllipse(
                    michel_tick,
                    michel_channel,
                    tick_radius,
                    channel_radius
                );
                e->SetFillStyle(0);
                e->SetLineWidth(1);
                e->SetLineColor(kRed);
                e->Draw();

                if (v) cout << "\t\tDrawn ellipse at (" 
                            << michel_tick << "," << michel_channel << ") " 
                            << "of radii (" << tick_radius << "," << channel_radius << ")" << endl;
            } // end loop over michel hits
        } // end loop over sections
        cPlanes.at(plane)->Write();
    } // end loop over planes

    TCanvas* cE = tfs->make<TCanvas>("cE","Michel Electron Energy");
    cE->Divide(2,1);
    cE->cd(1);
    hTrueE->Draw();
    cE->cd(2);
    hRecoE->Draw();
    cE->Write();

    if (v) cout << "\033[94m" << "End of Drawing==================================================================" << "\033[0m" << endl;
}

// pair<int,int> ana::Michecks::GetWireTPC(geo::Point_t *pt) {
//     geo::TPCID tpc = fGeom->FindTPCAtPosition(*pt);
//     geo::WireID wire = fGeom->NearestWireID(*pt, tpc);
//     return {tpc.TPC, wire.Wire};
// }

bool ana::Michecks::IsInVolume(double x, double y, double z, double eps) {

    geo::CryostatID cryoid{0};
    geo::TPCID tpc0{cryoid, 0}, tpcN{cryoid, fGeom->NTPC(cryoid)-1};
    return (x-eps > fGeom->TPC(tpc0).MinX() and x+eps < fGeom->TPC(tpcN).MaxX() and
            y-eps > fGeom->TPC(tpc0).MinY() and y+eps < fGeom->TPC(tpcN).MaxY() and
            z-eps > fGeom->TPC(tpc0).MinZ() and z+eps < fGeom->TPC(tpcN).MaxZ());   
}
bool ana::Michecks::IsInVolume(float x, float y, float z, float eps) {
    return IsInVolume(double(x), double(y), double(z), double(eps));
}
bool ana::Michecks::IsInVolume(TLorentzVector const& V) {
    return IsInVolume(V.X(), V.Y(), V.Z(), 0.);
}
bool ana::Michecks::IsInVolume(geo::Point_t const& P) {
    return IsInVolume(P.X(), P.Y(), P.Z(), 0.);
}
bool ana::Michecks::IsUpright(recob::Track const& T) {
    return T.Start().X() > T.End().X();
}
int ana::Michecks::ChannelSection(unsigned int ch) {
    return ch / (fGeom->Nchannels()/4) ;
}

DEFINE_ART_MODULE(ana::Michecks)
