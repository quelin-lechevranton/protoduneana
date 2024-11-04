/*
 * Based on Dumper module by Laura Perez Molina
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
// #include "Math/GenVector/LorentzVector.h"
// #include "Math/GenVector/PositionVector3D.h"
// #include <TVector3.h>

// std includes
#include <vector>
#include <iterator> 
#include <string>

using namespace std;

namespace ana { class OmegaDumper; }

class ana::OmegaDumper : public art::EDAnalyzer {
public:

    explicit OmegaDumper(fhicl::ParameterSet const& fcl);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    OmegaDumper(OmegaDumper const &) = delete;
    OmegaDumper(OmegaDumper &&) = delete;
    OmegaDumper & operator = (OmegaDumper const &) = delete;
    OmegaDumper & operator = (OmegaDumper &&) = delete;
    
    // Required functions.
    void analyze(art::Event const & evt) override; 
    void beginJob() override;
    // void endJob()   override;
    void reset();

private:

    // vector<string> fTrees;
    bool fIsTruth;
    bool fIsReco;
    vector<vector<string>> fProducts;

    TTree *fTruth;
    TTree *fReco;

    art::InputTag   tag_tru,
                    tag_prt,
                    tag_dep,
                    tag_pfp,
                    tag_clu,
                    tag_trk,
                    tag_cal,
                    tag_shw,
                    tag_spt,
                    tag_hit;

    //Events
    size_t fEvent, fRun, fSubRun;

    //MCTruth
    size_t fNPrt;

    //MCParticle
    vector<int> fPrtPdg;
    vector<int> fPrtMomID;
    vector<size_t>  fPrtNPt,
                    fPrtNDau;
    vector<vector<double>>  fPrtDauID, //should be int
                            fPrtX,
                            fPrtY,
                            fPrtZ,
                            fPrtT,
                            fPrtPx, 
                            fPrtPy, 
                            fPrtPz, 
                            fPrtP, 
                            fPrtE;

    //SimEnergyDeposit
    vector<size_t>  fPrtNDep;
    vector<vector<double>>  fDepPdg, //should be int
                            fDepX,
                            fDepY,
                            fDepZ,
                            fDepT,
                            fDepE;

    //PFParticle
    size_t fNPfp;
    vector<int> fPfpTrkID,
                fPfpShwID;
    vector<size_t>  fPfpNClu,
                    fPfpNSpt,
                    fPfpNHit;
    vector<vector<double>>  fPfpCluKey, //should be ?
                            fPfpSptKey;

    //Track
    size_t fNTrk;
    vector<size_t>  fTrkKey,
                    fTrkNHit,
                    fTrkNPt;
    vector<double>  fTrkLength;
    vector<vector<double>>  fTrkPtX, //should be Double_t
                            fTrkPtY,
                            fTrkPtZ,
                            fTrkDirX,
                            fTrkDirY,
                            fTrkDirZ,
                            fTrkHitKey;

    //Calorimetry
    vector<size_t>  fTrkCalPlane,   //should be unsigned int
                    fTrkCalNPt;
    vector<float>   fTrkCalRange,
                    fTrkCalKineticE;
    vector<vector<double>>  fTrkCaldEdx, //should be float
                            fTrkCaldQdx, //should be float
                            fTrkCalResRange; //should be float

    //Shower
    size_t fNShw;
    vector<int> fShwID; //useless?

    //Cluster
    vector<size_t>  fCluKey,
                    fCluNHit,
                    fCluPlane; //should be unsigned int
    vector<float>   fCluIntegral,
                    fCluSumADC,
                    fCluWidth;
    vector<vector<double>>  fCluHitKey; //should be ?

    //SpacePoint
    vector<size_t>  fSptKey,
                    fSptNHit;
    vector<double>  fSptX, //should be Double32_t
                    fSptY,
                    fSptZ;
    vector<vector<double>>  fSptHitKey; //should be ?

    //Hit
    size_t fNHit;
    vector<size_t>  fHitKey,
                    fHitNTrk,
                    fHitNSpt,
                    fHitNClu,  
                    fHitPlane; //should be unsigned int
    vector<float>   fHitSumADC, //should be float
                    fHitIntegral; //should be float
    vector<vector<double>>  fHitTrkKey, //should be ?
                            fHitSptKey,
                            fHitCluKey;

};

ana::OmegaDumper::OmegaDumper(fhicl::ParameterSet const & fcl) :
    EDAnalyzer{fcl} 
{
    art::ServiceHandle<art::TFileService> tfs;
    fIsTruth =  fcl.get<int>("Truth",0);
    fIsReco =   fcl.get<int>("Reco",0);
    fProducts = fcl.get<vector<vector<string>>>("Products");

    for (vector<string> prod : fProducts) {

        string  label=prod[0],
                instance=prod[1],
                object=prod[2],
                process=prod[3];

        if (fIsTruth) {
            if (object == "simb::MCTruth")              tag_tru= art::InputTag(label,instance);
            else if (object == "simb::MCParticle")      tag_prt= art::InputTag(label,instance);
            else if (object == "sim::SimEnergyDeposit") tag_dep= art::InputTag(label,instance);
        }

        if (fIsReco) {
            if (object == "recob::PFParticle")          tag_pfp= art::InputTag(label,instance);
            else if (object == "recob::Cluster")        tag_clu= art::InputTag(label,instance);
            else if (object == "recob::Track")          tag_trk= art::InputTag(label,instance);
            else if (object == "anab::Calorimetry")     tag_cal= art::InputTag(label,instance);
            else if (object == "recob::Shower")         tag_shw= art::InputTag(label,instance);
            else if (object == "recob::SpacePoint")     tag_spt= art::InputTag(label,instance);
            else if (object == "recob::Hit")            tag_hit= art::InputTag(label,instance);
        }
    } //end prod loop

    if (fIsTruth) {

        fTruth = tfs->make<TTree>("Truth","Truth");

        //Events
        fTruth->Branch("fEvent",        &fEvent);
        fTruth->Branch("fRun",          &fRun);
        fTruth->Branch("fSubrun",       &fSubRun);

        //MCTruth

        //MCParticle
        fTruth->Branch("fNPrt",         &fNPrt);
        fTruth->Branch("fPrtPdg",       &fPrtPdg); 
        fTruth->Branch("fPrtMomID",     &fPrtMomID); 
        fTruth->Branch("fPrtNDau",      &fPrtNDau); 
        fTruth->Branch("fPrtDauID",     &fPrtDauID); 
        fTruth->Branch("fPrtNPt",       &fPrtNPt); 
        fTruth->Branch("fPrtX",         &fPrtX); 
        fTruth->Branch("fPrtY",         &fPrtY); 
        fTruth->Branch("fPrtZ",         &fPrtZ); 
        //fTruth->Branch("fPrtT",         &fPrtT); 
        //fTruth->Branch("fPrtPx",        &fPrtPx); 
        //fTruth->Branch("fPrtPy",        &fPrtPy); 
        //fTruth->Branch("fPrtPz",        &fPrtPz); 
        fTruth->Branch("fPrtP",         &fPrtP); 
        fTruth->Branch("fPrtE",         &fPrtE); 

        //SimEnergyDeposit
        fTruth->Branch("fPrtNDep",      &fPrtNDep);
        fTruth->Branch("fDepPdg",       &fDepPdg); 
        fTruth->Branch("fDepX",         &fDepX); 
        fTruth->Branch("fDepY",         &fDepY); 
        fTruth->Branch("fDepZ",         &fDepZ); 
        //fTruth->Branch("fDepT",         &fDepT); 
        fTruth->Branch("fDepE",         &fDepE); 

    } //end Truth tree

    if (fIsReco) {

        fReco = tfs->make<TTree>("Reco","Reco");

        //Events
        fReco->Branch("fEvent",         &fEvent);
        fReco->Branch("fRun",           &fRun);
        fReco->Branch("fSubrun",        &fSubRun);

        //PFParticle
        fReco->Branch("fNPfp",          &fNPfp);
        fReco->Branch("fPfpTrkID",      &fPfpTrkID);
        //fReco->Branch("fPfpShwID",      &fPfpShwID);
        //fReco->Branch("fPfpNClu",       &fPfpNClu);
        fReco->Branch("fPfpNSpt",       &fPfpNSpt);
        fReco->Branch("fPfpNHit",       &fPfpNHit);
        //fReco->Branch("fPfpCluKey",     &fPfpCluKey);
        fReco->Branch("fPfpSptKey",     &fPfpSptKey);

        //Track
        fReco->Branch("fNTrk",          &fNTrk);
        fReco->Branch("fTrkKey",        &fTrkKey);
        fReco->Branch("fTrkNHit",       &fTrkNHit);
        fReco->Branch("fTrkLength",     &fTrkLength);
        fReco->Branch("fTrkNPt",        &fTrkNPt);
        fReco->Branch("fTrkPtX",        &fTrkPtX);
        fReco->Branch("fTrkPtY",        &fTrkPtY);
        fReco->Branch("fTrkPtZ",        &fTrkPtZ);
        //fReco->Branch("fTrkDirX",       &fTrkDirX);
        //fReco->Branch("fTrkDirY",       &fTrkDirY);
        //fReco->Branch("fTrkDirZ",       &fTrkDirZ);
        fReco->Branch("fTrkHitKey",     &fTrkHitKey);

        //Calorimetry
        fReco->Branch("fTrkCalPlane",   &fTrkCalPlane);
        fReco->Branch("fTrkCalRange",   &fTrkCalRange);
        //fReco->Branch("fTrkCalKineticE",&fTrkCalKineticE);
        fReco->Branch("fTrkCalNPt",     &fTrkCalNPt);
        fReco->Branch("fTrkCaldEdx",    &fTrkCaldEdx);
        fReco->Branch("fTrkCaldQdx",    &fTrkCaldQdx);
        fReco->Branch("fTrkCalResRange",&fTrkCalResRange);

        //Shower
        fReco->Branch("fNShw",          &fNShw);
        fReco->Branch("fShwID",         &fShwID);

        //Cluster
        //fReco->Branch("fCluKey",        &fCluKey);
        //fReco->Branch("fCluNHit",       &fCluNHit);
        //fReco->Branch("fCluPlane",      &fCluPlane);
        //fReco->Branch("fCluIntegral",   &fCluIntegral);
        //fReco->Branch("fCluSumADC",     &fCluSumADC);
        //fReco->Branch("fCluWidth",      &fCluWidth);
        //fReco->Branch("fCluHitKey",     &fCluHitKey);

        //SpacePoint
        fReco->Branch("fSptKey",        &fSptKey);
        fReco->Branch("fSptNHit",       &fSptNHit);
        fReco->Branch("fSptX",          &fSptX);
        fReco->Branch("fSptY",          &fSptY);
        fReco->Branch("fSptZ",          &fSptZ);
        fReco->Branch("fSptHitKey",     &fSptHitKey);

        //Hit
        fReco->Branch("fNHit",          &fNHit);
        fReco->Branch("fHitKey",        &fHitKey);
        fReco->Branch("fHitNTrk",       &fHitNTrk);
        fReco->Branch("fHitNSpt",       &fHitNSpt);
        //fReco->Branch("fHitNClu",       &fHitNClu);
        fReco->Branch("fHitPlane",      &fHitPlane);
        fReco->Branch("fHitSumADC",     &fHitSumADC);
        fReco->Branch("fHitIntegral",   &fHitIntegral);
        fReco->Branch("fHitTrkKey",     &fHitTrkKey);
        fReco->Branch("fHitSptKey",     &fHitSptKey);
        //fReco->Branch("fHitCluKey",     &fHitCluKey);

    } //end Reco tree

} //end OmegaDumper()

void ana::OmegaDumper::beginJob() {} //end beginJob()

void ana::OmegaDumper::analyze(const art::Event & evt) {
    reset();
    fEvent  = evt.id().event(); 
    fRun    = evt.id().run();
    fSubRun = evt.id().subRun();

    //Truth Tree
    if (fIsTruth) {

    //simb::MCTruth
    // if (!tag_tru.empty()) {
    // art::ValidHandle<vector<simb::MCTruth>> const vh_tru = evt.getValidHandle<vector<simb::MCTruth>>(tag_tru);

    // for (simb::MCTruth const & tru : *vh_tru) {

    // } //end vh_tru loop
    // } //end simb::MCTruth

    //simb::MCParticle
    if (!tag_prt.empty()) {
    art::ValidHandle<vector<simb::MCParticle>> const vh_prt = evt.getValidHandle<vector<simb::MCParticle>>(tag_prt);
    
    fNPrt=vh_prt->size();

    for (simb::MCParticle const & prt : *vh_prt) {

        fPrtPdg.        push_back(prt.PdgCode());
        fPrtMomID.      push_back(prt.Mother()-1); //IDs start at 1, we want 0 for vector indexing
        fPrtNDau.       push_back(prt.NumberDaughters());

        vector<double> tpPrtDauID;
        for (int i_dau=0; i_dau < prt.NumberDaughters(); i_dau++) {
            tpPrtDauID.push_back(prt.Daughter(i_dau)-1); //IDs start at 1, we want 0 for vector indexing
        }
        fPrtDauID.   push_back(tpPrtDauID);

        fPrtNPt.        push_back(prt.NumberTrajectoryPoints());

        vector<double> tpPrtX,tpPrtY,tpPrtZ,tpPrtT,tpPrtPx,tpPrtPy,tpPrtPz,tpPrtP,tpPrtE;
        for (size_t i_ppt=0; i_ppt < prt.NumberTrajectoryPoints(); i_ppt++) {
            tpPrtX.     push_back(prt.Vx(i_ppt));
            tpPrtY.     push_back(prt.Vy(i_ppt));
            tpPrtZ.     push_back(prt.Vz(i_ppt));
            tpPrtT.     push_back(prt.T(i_ppt));
            tpPrtPx.    push_back(prt.Px(i_ppt));
            tpPrtPy.    push_back(prt.Py(i_ppt));
            tpPrtPz.    push_back(prt.Pz(i_ppt));
            tpPrtP.     push_back(prt.P(i_ppt));
            tpPrtE.     push_back(prt.E(i_ppt));
        }
        fPrtX.          push_back(tpPrtX);
        fPrtY.          push_back(tpPrtY);
        fPrtZ.          push_back(tpPrtZ);
        fPrtT.          push_back(tpPrtT);
        fPrtPx.         push_back(tpPrtPx);
        fPrtPy.         push_back(tpPrtPy);
        fPrtPz.         push_back(tpPrtPz);
        fPrtP.          push_back(tpPrtP);
        fPrtE.          push_back(tpPrtE);

        //sim::SimEnergyDeposit
        if (!tag_dep.empty()) {
        art::ValidHandle<vector<sim::SimEnergyDeposit>> const vh_dep = evt.getValidHandle<vector<sim::SimEnergyDeposit>>(tag_dep);

        vector<double> tpDepPdg,tpDepX,tpDepY,tpDepZ,tpDepT,tpDepE;
        for (sim::SimEnergyDeposit const & dep : *vh_dep ) {

            if (dep.TrackID()!=prt.TrackId()) continue;

            tpDepPdg.   push_back(dep.PdgCode());
            tpDepX.     push_back(dep.X());
            tpDepY.     push_back(dep.Y());
            tpDepZ.     push_back(dep.Z());
            tpDepT.     push_back(dep.T());
            tpDepE.     push_back(dep.E());

        } //end vh_dep loop
        fDepPdg.        push_back(tpDepPdg);
        fDepX.          push_back(tpDepX);
        fDepY.          push_back(tpDepY);
        fDepZ.          push_back(tpDepZ);
        fDepT.          push_back(tpDepT);
        fDepE.          push_back(tpDepE);

        fPrtNDep.       push_back(tpDepPdg.size());

        } //sim::SimEnergyDeposit
    } //end vh_prt loop
    } //end simb::MCParticle

    fTruth->Fill();
    } //end Truth

    //Reco Tree
    if (fIsReco) {

    //PFParticle
    auto const vh_pfp = evt.getValidHandle<vector<recob::PFParticle>>(tag_pfp);
    
    art::FindManyP<recob::Track>        const fmp_pfp_trk(vh_pfp,evt,tag_trk);
    art::FindManyP<recob::Shower>       const fmp_pfp_shw(vh_pfp,evt,tag_shw);
    art::FindManyP<recob::Cluster>      const fmp_pfp_clu(vh_pfp,evt,tag_clu);
    art::FindManyP<recob::SpacePoint>   const fmp_pfp_spt(vh_pfp,evt,tag_spt);
    
    auto const vh_trk = evt.getValidHandle<vector<recob::Track>>(tag_trk);
    art::FindManyP<anab::Calorimetry>   const fmp_trk_cal(vh_trk,evt,tag_cal);
    art::FindManyP<recob::Hit>          const fmp_trk_hit(vh_trk,evt,tag_trk);

    auto const vh_clu = evt.getValidHandle<vector<recob::Cluster>>(tag_clu);
    art::FindManyP<recob::Hit>          const fmp_clu_hit(vh_clu,evt,tag_clu);
    
    auto const vh_spt = evt.getValidHandle<vector<recob::SpacePoint>>(tag_spt);
    art::FindManyP<recob::Hit>          const fmp_spt_hit(vh_spt,evt,tag_spt);

    vector<art::Ptr<recob::PFParticle>> vp_pfp;
    art::fill_ptr_vector(vp_pfp,vh_pfp);
    
    fNPfp = vp_pfp.size();

    for (art::Ptr<recob::PFParticle> const & p_pfp : vp_pfp) {

        //Track from PFParticle
        vector<art::Ptr<recob::Track>> const vp_trk = fmp_pfp_trk.at(p_pfp.key());

        if (vp_trk.empty()) {

            fPfpTrkID.      push_back(-1);

        }
        else {
            art::Ptr<recob::Track> const p_trk = vp_trk[0];

            fPfpTrkID.      push_back(p_trk->ID());

            fTrkKey.        push_back(p_trk.key());
            fTrkLength.     push_back(p_trk->Length());

            vector<double>  tpTrkPtX,tpTrkPtY,tpTrkPtZ,
                            tpTrkDirX,tpTrkDirY,tpTrkDirZ;
            for (size_t i_tpt = p_trk->FirstPoint(); i_tpt < p_trk->LastPoint(); i_tpt++) {
                if (!p_trk->HasValidPoint(i_tpt)) {continue;}

                tpTrkPtX.   push_back(p_trk->LocationAtPoint(i_tpt).X());
                tpTrkPtY.   push_back(p_trk->LocationAtPoint(i_tpt).Y());
                tpTrkPtZ.   push_back(p_trk->LocationAtPoint(i_tpt).Z());
                tpTrkDirX.  push_back(p_trk->MomentumVectorAtPoint(i_tpt).X());
                tpTrkDirY.  push_back(p_trk->MomentumVectorAtPoint(i_tpt).Y());
                tpTrkDirZ.  push_back(p_trk->MomentumVectorAtPoint(i_tpt).Z());
            }
            fTrkPtX.        push_back(tpTrkPtX);
            fTrkPtY.        push_back(tpTrkPtY);
            fTrkPtZ.        push_back(tpTrkPtZ);
            fTrkDirX.       push_back(tpTrkDirX);
            fTrkDirY.       push_back(tpTrkDirY);
            fTrkDirZ.       push_back(tpTrkDirZ);

            fTrkNPt.        push_back(tpTrkPtX.size());

            //Calorimetry from Track
            vector<art::Ptr<anab::Calorimetry>> const vp_cal = fmp_trk_cal.at(p_trk.key());

            for (art::Ptr<anab::Calorimetry> const & p_cal : vp_cal) {
                if (p_cal->PlaneID().Plane != geo::kW) continue; 

                fTrkCalPlane.   push_back(p_cal->PlaneID().Plane);
                fTrkCalRange.   push_back(p_cal->Range());
                fTrkCalKineticE.push_back(p_cal->KineticEnergy());
                fTrkCalNPt.     push_back(p_cal->dEdx().size());
                // fTrkCaldEdx.    push_back(p_cal->dEdx());
                // fTrkCaldQdx.    push_back(p_cal->dQdx());
                // fTrkCalResRange.push_back(p_cal->ResidualRange());

                vector<double> tpTrkCaldEdx,tpTrkCaldQdx,tpTrkCalResRange;
                for (double dEdx : p_cal->dEdx()) {tpTrkCaldEdx.push_back(dEdx);}
                for (double dQdx : p_cal->dQdx()) {tpTrkCaldQdx.push_back(dQdx);}
                for (double ResR : p_cal->ResidualRange()) {tpTrkCalResRange.push_back(ResR);}
                fTrkCaldEdx.    push_back(tpTrkCaldEdx);
                fTrkCaldQdx.    push_back(tpTrkCaldQdx);
                fTrkCalResRange.push_back(tpTrkCalResRange);
            } //end v_cal loop

            //Hit from Tracks
            vector<art::Ptr<recob::Hit>> const vp_hit = fmp_trk_hit.at(p_trk.key());

            fTrkNHit.       push_back(vp_hit.size());

            vector<double> tpTrkHitKey;
            for (art::Ptr<recob::Hit> p_hit : vp_hit) {
                tpTrkHitKey.    push_back(p_hit.key());
            } //end vp_hit loop
            fTrkHitKey.push_back(tpTrkHitKey);

        } //end IsTrack condition

        //Shower
        vector<art::Ptr<recob::Shower>> const vp_shw = fmp_pfp_shw.at(p_pfp.key());

        if (vp_shw.empty()) {
            
            fPfpShwID.      push_back(-1);
        }
        else {
            art::Ptr<recob::Shower> const p_shw = vp_shw[0];

            fPfpShwID.      push_back(p_shw->ID());

            fShwID.         push_back(p_shw->ID());
        } //end IsShower condition

        //Cluster from PFParticle
        vector<art::Ptr<recob::Cluster>> const vp_clu = fmp_pfp_clu.at(p_pfp.key());

        fPfpNClu.           push_back(vp_clu.size());

        vector<double> tpPfpCluKey;
        for (art::Ptr<recob::Cluster> p_clu : vp_clu) {
            tpPfpCluKey.    push_back(p_clu.key());

            fCluKey.        push_back(p_clu.key());
            fCluNHit.       push_back(p_clu->NHits());
            fCluPlane.      push_back(p_clu->Plane().Plane);
            fCluIntegral.   push_back(p_clu->Integral());
            fCluSumADC.     push_back(p_clu->SummedADC());
            fCluWidth.      push_back(p_clu->Width());

            //Hit from Cluster
            vector<art::Ptr<recob::Hit>> const vp_hit = fmp_clu_hit.at(p_clu.key());

            vector<double> tpCluHitKey;
            for (art::Ptr<recob::Hit> p_hit : vp_hit) {
                tpCluHitKey.push_back(p_hit.key());
            } //end vp_hit loop
            fCluHitKey.     push_back(tpCluHitKey);

        } //end vp_clu loop
        fPfpCluKey.         push_back(tpPfpCluKey);

        //SpacePoint from PFParticle
        vector<art::Ptr<recob::SpacePoint>> const vp_spt = fmp_pfp_spt.at(p_pfp.key());

        fPfpNSpt.           push_back(vp_spt.size());

        size_t tpPfpNHit=0;
        vector<double> tpPfpSptKey;
        for (art::Ptr<recob::SpacePoint> p_spt : vp_spt) {
            tpPfpSptKey.    push_back(p_spt.key());

            fSptKey.        push_back(p_spt.key());
            fSptX.          push_back(p_spt->position().X());
            fSptY.          push_back(p_spt->position().Y());
            fSptZ.          push_back(p_spt->position().Z());

            //Hit from SpacePoint
            vector<art::Ptr<recob::Hit>> const vp_hit = fmp_spt_hit.at(p_spt.key());

            fSptNHit.       push_back(vp_hit.size());
            tpPfpNHit+=vp_hit.size();

            vector<double> tpSptHitKey;
            for (art::Ptr<recob::Hit> p_hit : vp_hit) {
                tpSptHitKey.    push_back(p_hit.key());
            } //end vp_hit loop
            fSptHitKey.push_back(tpSptHitKey);

        } //end vp_spt loop
        fPfpSptKey.push_back(tpPfpSptKey);
        fPfpNHit.push_back(tpPfpNHit);
    } //end vp_pfp loop
    fNTrk=fTrkLength.size();
    fNShw=fShwID.size();
    //end PFParticle

    //Hit
    auto const vh_hit = evt.getValidHandle<vector<recob::Hit>>(tag_hit);
    art::FindManyP<recob::Track>        const fmp_hit_trk(vh_hit,evt,tag_trk);
    art::FindManyP<recob::SpacePoint>   const fmp_hit_spt(vh_hit,evt,tag_spt);
    art::FindManyP<recob::Cluster>      const fmp_hit_clu(vh_hit,evt,tag_clu);

    vector<art::Ptr<recob::Hit>> vp_hit;
    art::fill_ptr_vector(vp_hit,vh_hit);

    fNHit=vp_hit.size();

    for (art::Ptr<recob::Hit> p_hit : vp_hit) {
        fHitKey.        push_back(p_hit.key());
        fHitPlane.      push_back(p_hit->View());
        fHitSumADC.     push_back(p_hit->SummedADC());
        fHitIntegral.   push_back(p_hit->Integral());

        try {
            vector<art::Ptr<recob::Track>> const vp_trk = fmp_hit_trk.at(p_hit.key());

            fHitNTrk.       push_back(vp_trk.size());

            vector<double> tpHitTrkKey;
            for (art::Ptr<recob::Track> p_trk : vp_trk) {
                tpHitTrkKey.push_back(p_trk.key());
            } //end vp_spt loop
            fHitTrkKey.     push_back(tpHitTrkKey);
        } 
        catch (...) {
            fHitNTrk.       push_back(0);
        } //end try/catch
        try {
            vector<art::Ptr<recob::SpacePoint>> const vp_spt = fmp_hit_spt.at(p_hit.key());

            fHitNSpt.       push_back(vp_spt.size());

            vector<double> tpHitSptKey;
            for (art::Ptr<recob::SpacePoint> p_spt : vp_spt) {
                tpHitSptKey.push_back(p_spt.key());
            } //end vp_spt loop
            fHitSptKey.     push_back(tpHitSptKey);
        } 
        catch (...) {
            fHitNSpt.       push_back(0);
        } //end try/catch
        try {
            vector<art::Ptr<recob::Cluster>> const vp_clu = fmp_hit_clu.at(p_hit.key());

            fHitNClu.       push_back(vp_hit.size());

            vector<double> tpHitCluKey;
            for (art::Ptr<recob::Cluster> p_clu : vp_clu) {
                tpHitCluKey.push_back(p_clu.key());
            } //end vp_hit loop
            fHitCluKey.     push_back(tpHitCluKey);
        }
        catch (...) {
            fHitNClu.       push_back(0);
        } //end try/catch
    } //end vp_hit loop
    //end recob::Hit
    
    fReco->Fill();
    } //end Reco

} //end analyze()

void ana::OmegaDumper::reset() {
    
    // Events
    fEvent=0;
    fRun=0;
    fSubRun=0;

    //MCTruth
    fNPrt=0;

    //MCParticle
    fPrtPdg.clear();
    fPrtMomID.clear();
    fPrtNDau.clear();
    fPrtDauID.clear();
    fPrtNPt.clear();
    fPrtX.clear();
    fPrtY.clear();
    fPrtZ.clear();
    fPrtT.clear();
    fPrtPx.clear();
    fPrtPy.clear();
    fPrtPz.clear();
    fPrtP.clear();
    fPrtE.clear();

    //SimEnergyDeposit
    fPrtNDep.clear();
    fDepPdg.clear();
    fDepX.clear();
    fDepY.clear();
    fDepZ.clear();
    fDepT.clear();
    fDepE.clear();

    //PFParticle
    fNPfp=0;
    fPfpTrkID.clear();
    fPfpShwID.clear();
    fPfpNClu.clear();
    fPfpNSpt.clear();
    fPfpNHit.clear();
    fPfpCluKey.clear();
    fPfpSptKey.clear();

    //Track
    fNTrk=0;
    fTrkKey.clear();
    fTrkNHit.clear();
    fTrkNPt.clear();
    fTrkLength.clear();
    fTrkPtX.clear();
    fTrkPtY.clear();
    fTrkPtZ.clear();
    fTrkHitKey.clear();

    //Calorimetry
    fTrkCalPlane.clear();
    fTrkCalRange.clear();
    fTrkCalKineticE.clear();
    fTrkCalNPt.clear();
    fTrkCaldEdx.clear();
    fTrkCaldQdx.clear();
    fTrkCalResRange.clear();

    //Shower
    fNShw=0;
    fShwID.clear();

    //Cluster
    fCluKey.clear();
    fCluNHit.clear();
    fCluPlane.clear();
    fCluIntegral.clear();
    fCluSumADC.clear();
    fCluWidth.clear();
    fCluHitKey.clear();

    //SpacePoint
    fSptKey.clear();
    fSptNHit.clear();
    fSptX.clear();
    fSptY.clear();
    fSptZ.clear();
    fSptHitKey.clear();

    //Hit
    fNHit=0;
    fHitKey.clear();
    fHitNTrk.clear();
    fHitNSpt.clear();
    fHitNClu.clear();
    fHitPlane.clear();
    fHitSumADC.clear();
    fHitIntegral.clear();
    fHitSptKey.clear();
    fHitCluKey.clear();
    fHitTrkKey.clear();

} //end reset()

DEFINE_ART_MODULE(ana::OmegaDumper)
