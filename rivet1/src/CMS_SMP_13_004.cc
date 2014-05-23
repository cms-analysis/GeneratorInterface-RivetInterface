// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Math/MathUtils.hh"

/// @todo Include more projections as required, e.g. ChargedFinalState, FastJets, ZFinder...

namespace Rivet {


  class CMS_SMP_13_004 : public Analysis {
  public:
    /// Constructor
    CMS_SMP_13_004()
      : Analysis("CMS_SMP_13_004")
    {    }

  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      /// @todo Initialise and register projections here
      const FinalState fs(-MAXRAPIDITY,MAXRAPIDITY);
      addProjection(fs, "FS");
      
      vector<pair<PdgId,PdgId> > vidsZ;
      vidsZ.push_back(make_pair(ELECTRON, POSITRON));
      vidsZ.push_back(make_pair(MUON, ANTIMUON));

      FinalState fsZ(-MAXRAPIDITY,MAXRAPIDITY);
      InvMassFinalState invfsZ(fsZ, vidsZ, 60*GeV, 120*GeV);
      addProjection(invfsZ, "INVFSZ");
 
      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(invfsZ); ///check genjet def in PAT
      addProjection(vfs, "VFS");
      addProjection(FastJets(vfs, FastJets::ANTIKT, 0.5), "Jets");

      UnstableFinalState ufs;
      addProjection(ufs, "UFS");
      
      ////example for the tutorial where fastjets are final state:
      //FastJets fj(fs, FastJets::ANTIKT, 0.5);
      //fj.useInvisibles();
      //addProjection(fj, "Jets");
      
      /// @todo Book histograms here, e.g.:
      // _h_XXXX = bookProfile1D(1, 1, 1);
      // _h_YYYY = bookHistogram1D(2, 1, 1);
      _h_Z_pT_normalised = bookHistogram1D(1, 1, 1);
      //_h_jet1_pT_normalised = bookHistogram1D(1, 1, 1);

    }

    bool ApplyElectronCutsForZee(double pt1, double pt2, double eta1, double eta2){
      bool isFid1 = ((fabs(eta1)<1.4442)||((fabs(eta1)>1.566)&&(fabs(eta1)<2.4)));
      bool isFid2 = ((fabs(eta2)<1.4442)||((fabs(eta2)>1.566)&&(fabs(eta2)<2.4)));
      if( isFid1 && isFid2 && pt1>20 && pt2 >20) return true;
      else return false;
    }
    
    bool ApplyMuonCutsForZmm(double pt1, double pt2, double eta1, double eta2){
      bool isFid1 = ((fabs(eta1)<2.4));
      bool isFid2 = ((fabs(eta2)<2.4));
      if( isFid1 && isFid2 && pt1>20 && pt2 >20) return true;
      else return false;
    } 

    vector<FourMomentum> match(vector<FourMomentum> Vec1,vector<FourMomentum> Vec2,double dr){
      vector<FourMomentum> matching;
      vector<FourMomentum> vector1=Vec1;
      map <double, FourMomentum > in_cone;
      foreach (FourMomentum vec2, Vec2) { 
	foreach (FourMomentum vec1, vector1) { 
	  double deltar = sqrt(pow((vec1.eta() - vec2.eta()),2) + pow((vec1.phi()-vec2.phi()),2));
	  if (deltar < dr) {
	    in_cone[deltar]=vec1;
	  }
	}
	if(in_cone.size()>0){
	FourMomentum bestmatch;
	bestmatch = in_cone.begin()->second;
	matching.push_back(bestmatch);
	//v.erase(std::remove(v.begin(), v.end(), 99), v.end()); 
	//vector1.erase(bestmatch);
	}	
      }      
      return matching;
    }
    
    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      /// @todo Do the event by event analysis here

      const InvMassFinalState& invMassFinalStateZ = applyProjection<InvMassFinalState>(event, "INVFSZ");
      const UnstableFinalState& ufs = applyProjection<UnstableFinalState>(event, "UFS");

      bool isZmm(false); bool isZee(false);
      bool isZ(false);
      isZ  = (!(invMassFinalStateZ.empty())); //&& invMassFinalStateW.empty());

      const ParticleVector&  ZDecayProducts =  invMassFinalStateZ.particles();

      if (ZDecayProducts.size() < 2) vetoEvent;
      
      double pt1=-9999.,  pt2=-9999.;
      double phi1=-9999., phi2=-9999.;
      double eta1=-9999., eta2=-9999.;
      
      if(isZ){
	pt1  = ZDecayProducts[0].momentum().pT();
	pt2  = ZDecayProducts[1].momentum().pT();
	eta1 = ZDecayProducts[0].momentum().eta();
	eta2 = ZDecayProducts[1].momentum().eta();
	phi1 = ZDecayProducts[0].momentum().phi();
	phi2 = ZDecayProducts[1].momentum().phi();
      }
      isZmm = isZ && ((fabs(ZDecayProducts[0].pdgId()) == 13) && (fabs(ZDecayProducts[1].pdgId()) == 13));
      isZee = isZ && ((fabs(ZDecayProducts[0].pdgId()) == 11) && (fabs(ZDecayProducts[1].pdgId()) == 11));
      
      if(!((isZmm||isZee))){
	cout << "vetoEvent" << endl;
	vetoEvent;
      }

      FourMomentum Z;
      Z = ZDecayProducts[0].momentum()+ZDecayProducts[1].momentum();	       	

      bool passBosonConditions = false;
      if(isZmm)passBosonConditions = ApplyMuonCutsForZmm(pt1,pt2,eta1,eta2);
      if(isZee)passBosonConditions = ApplyElectronCutsForZee(pt1,pt2,eta1,eta2);
      if(!passBosonConditions)vetoEvent;

      //Obtain the jets.
      vector<FourMomentum> jet_list;
      foreach (const Jet& j, applyProjection<FastJets>(event,"Jets").jetsByPt(25.0*GeV)) {
	const double jeta = j.momentum().eta();
	const double jphi = j.momentum().phi();
	const double jpt = j.momentum().pT();
	if (fabs(jeta) < 2.1) 
	  if(jpt>25){
	    if(isZee||isZmm){
	      if (deltaR(pt1, phi1, jeta, jphi) > 0.3 && deltaR(pt2, phi2, jeta, jphi) > 0.3)
		{jet_list.push_back(j.momentum());}
	      continue;
	      }	      	
	  }
      }

      //only events with 2jets
      if (jet_list.size()<2)vetoEvent;

      //Obtain the bhadrons
      vector<FourMomentum> Bhadron_list;
      foreach (const Particle& p, ufs.particles()) {
	int pdgid = abs(p.pdgId());
	if ((499<pdgid && pdgid<600) || (4999<pdgid && pdgid<6000)){
	  Bhadron_list.push_back(p.momentum());}
      }

	//matched the jets to the Hadrons
      vector<FourMomentum> matchedBjet_list;
      matchedBjet_list = match(jet_list,Bhadron_list,0.5);
	


      if((isZee||isZmm)&& jet_list.size()>1 && matchedBjet_list.size()>1){
	double Zpt = Z.pT()/GeV;
	//_h_Z_pT_normalised->fill(Zpt, weight);
	cerr << "_h_Z_pT_normalised fill:" << _h_Z_pT_normalised << endl;
	//double jet1pT = jet_list[0].pT();
	//_h_jet1_pT_normalised->fill(jet1pT, weight);
      
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      /// @todo Normalise, scale and otherwise manipulate histograms here

      // scale(_h_YYYY, crossSection()/sumOfWeights()); # norm to cross section

      double pT_integral = _h_Z_pT_normalised->sumBinHeights();
      cerr << "pT_integral" << pT_integral << endl;
      
      // normalize(_h_YYYY); # normalize to unity
      ////normalize(_h_Z_pT_normalised,1.0);
      cerr << "_h_Z_pT_normalised finalize:" << _h_Z_pT_normalised << endl;

      //normalize(_h_jet1_pT_normalised,1.0);      
    }

  private:

    // Data members like post-cuts event weight counters go here

  private:
    /// @name Histograms
    //AIDA::IProfile1D *_h_XXXX;
    //AIDA::IHistogram1D *_h_YYYY;
    //AIDA::IHistogram1D* _histJetMultZelec;
    //AIDA::IHistogram1D* _histJetMultZmu;    
  AIDA::IHistogram1D* _h_Z_pT_normalised;
    //AIDA::IHistogram1D* _h_jet1_pT_normalised;

  };

  // The hook for the plugin system
  //DECLARE_RIVET_PLUGIN(CMS_SMP_13_004);
  AnalysisBuilder<CMS_SMP_13_004> plugin_CMS_SMP_13_004;

}

//Use e.g. 'rivet-buildplugin RivetCMS_SMP_13_004.so CMS_SMP_13_004.cc' to compile the plugin
