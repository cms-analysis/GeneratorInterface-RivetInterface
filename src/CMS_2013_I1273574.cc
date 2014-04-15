// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include <cmath>

/// @todo Include more projections as required, e.g. ChargedFinalState, FastJets, ZFinder...

namespace Rivet {


  class CMS_2013_I1273574 : public Analysis {
  
  private:
  
    AIDA::IHistogram1D *_hist_DeltaPtRelSoft;
    AIDA::IHistogram1D *_hist_DeltaPhiSoft;
    AIDA::IHistogram1D *_hist_DeltaS;
    
    AIDA::IHistogram1D *_hist_DeltaPtRelSoft_norm;
    AIDA::IHistogram1D *_hist_DeltaPhiSoft_norm;
    AIDA::IHistogram1D *_hist_DeltaS_norm;

    AIDA::IHistogram1D *_hist_LeadingHardJetPt;
    AIDA::IHistogram1D *_hist_SubLeadingHardJetPt;
    AIDA::IHistogram1D *_hist_LeadingSoftJetPt;
    AIDA::IHistogram1D *_hist_SubLeadingSoftJetPt;

    AIDA::IHistogram1D *_hist_LeadingHardJetEta;
    AIDA::IHistogram1D *_hist_SubLeadingHardJetEta;
    AIDA::IHistogram1D *_hist_LeadingSoftJetEta;
    AIDA::IHistogram1D *_hist_SubLeadingSoftJetEta;


  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CMS_2013_I1273574()
      : Analysis("CMS_2013_I1273574")
    {
      /// @todo Set whether your finalize method needs the generator cross section
      setBeams(PROTON, PROTON);
      setNeedsCrossSection(true);
    }

    //@}


    /// Book histograms and initialise projections before the run
    void init() {

      /// @todo Initialise and register projections here
      const FinalState cnfs(-4.7,4.7);
      addProjection(FastJets(cnfs, FastJets::ANTIKT, 0.5), "Jets");

      _hist_DeltaPtRelSoft = bookHistogram1D(7,1,1);
      _hist_DeltaPhiSoft = bookHistogram1D(5,1,1);
      _hist_DeltaS = bookHistogram1D(3,1,1);

      _hist_DeltaPtRelSoft_norm = bookHistogram1D(8,1,1);
      _hist_DeltaPhiSoft_norm = bookHistogram1D(6,1,1);
      _hist_DeltaS_norm = bookHistogram1D(4,1,1);
    
      _hist_LeadingHardJetPt = bookHistogram1D(2,1,1);
      _hist_SubLeadingHardJetPt = bookHistogram1D(14,1,1);
      _hist_LeadingSoftJetPt = bookHistogram1D(10,1,1);
      _hist_SubLeadingSoftJetPt = bookHistogram1D(12,1,1);

      _hist_LeadingHardJetEta = bookHistogram1D(1,1,1);
      _hist_SubLeadingHardJetEta = bookHistogram1D(13,1,1);
      _hist_LeadingSoftJetEta = bookHistogram1D(9,1,1);
      _hist_SubLeadingSoftJetEta = bookHistogram1D(11,1,1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      
      /// @todo Do the event by event analysis here
      
      unsigned int num_jets = 0;
      unsigned int num_Softjets = 0;
      const Jets jets = applyProjection<FastJets>(event, "Jets").jetsByPt(20*GeV);
      
      if (jets.size() < 4) vetoEvent;
      
      Jets jetAnalysis;
      Jets jetSoftAnalysis;
      foreach (const Jet& j, jets) {
	
	if(fabs(j.momentum().eta()) < 4.7 && j.momentum().pT() >= 50) {
	
	  num_jets +=1;
	  jetAnalysis.push_back(j);
	  	  
	}

	if(fabs(j.momentum().eta()) < 4.7 && j.momentum().pT() >= 20) {
	  
	  num_Softjets +=1;
	  jetSoftAnalysis.push_back(j);
          	  
        }
      }

      if (num_jets >= 2 && num_Softjets == 4) {

	double dphisoft=-100;
	dphisoft = deltaPhi(jetSoftAnalysis[3].momentum().azimuthalAngle(),jetSoftAnalysis[2].momentum().azimuthalAngle());
	_hist_DeltaPhiSoft->fill(dphisoft,weight);
	_hist_DeltaPhiSoft_norm->fill(dphisoft,weight);

	double vecsumlightsoft=-100;		
	vecsumlightsoft=sqrt(pow(jetSoftAnalysis[2].momentum().px()+jetSoftAnalysis[3].momentum().px(),2)+pow(jetSoftAnalysis[2].momentum().py()+jetSoftAnalysis[3].momentum().py(),2));
	double SptSoft=(vecsumlightsoft/(fabs(jetSoftAnalysis[2].momentum().pT())+fabs(jetSoftAnalysis[3].momentum().pT())));
	_hist_DeltaPtRelSoft->fill(SptSoft, weight);
	_hist_DeltaPtRelSoft_norm->fill(SptSoft, weight);

	double DPtHard = ((jetAnalysis[0].momentum().px()+jetAnalysis[1].momentum().px())*(jetAnalysis[0].momentum().px()+jetAnalysis[1].momentum().px())+(jetAnalysis[0].momentum().py()+jetAnalysis[1].momentum().py())*(jetAnalysis[0].momentum().py()+jetAnalysis[1].momentum().py()));
	double DPtSoft = ((jetSoftAnalysis[2].momentum().px()+jetSoftAnalysis[3].momentum().px())*(jetSoftAnalysis[2].momentum().px()+jetSoftAnalysis[3].momentum().px())+(jetSoftAnalysis[2].momentum().py()+jetSoftAnalysis[3].momentum().py())*(jetSoftAnalysis[2].momentum().py()+jetSoftAnalysis[3].momentum().py()));
	
	double Px = (jetAnalysis[0].momentum().px()+jetAnalysis[1].momentum().px())*(jetSoftAnalysis[2].momentum().px()+jetSoftAnalysis[3].momentum().px());
	double Py = (jetAnalysis[0].momentum().py()+jetAnalysis[1].momentum().py())*(jetSoftAnalysis[2].momentum().py()+jetSoftAnalysis[3].momentum().py());  
	
	double p1p2_mag = sqrt(DPtHard)*sqrt(DPtSoft);
	double DeltaS = acos((Px+Py)/p1p2_mag);
	
	_hist_DeltaS->fill(DeltaS,weight);
	_hist_DeltaS_norm->fill(DeltaS,weight);

	_hist_LeadingHardJetPt->fill(jetAnalysis[0].momentum().pT(), weight);
	_hist_SubLeadingHardJetPt->fill(jetAnalysis[1].momentum().pT(), weight);
	_hist_LeadingSoftJetPt->fill(jetSoftAnalysis[2].momentum().pT(), weight);
	_hist_SubLeadingSoftJetPt->fill(jetSoftAnalysis[3].momentum().pT(), weight);

	_hist_LeadingHardJetEta->fill(jetAnalysis[0].momentum().eta(), weight);
        _hist_SubLeadingHardJetEta->fill(jetAnalysis[1].momentum().eta(), weight);
	_hist_LeadingSoftJetEta->fill(jetSoftAnalysis[2].momentum().eta(), weight);
	_hist_SubLeadingSoftJetEta->fill(jetSoftAnalysis[3].momentum().eta(), weight);

	}    

    }

    /// Normalise histograms etc., after the run
    void finalize() {

      double invlumi = crossSection()/picobarn/sumOfWeights(); //norm to cross section
      
      scale(_hist_DeltaPtRelSoft, invlumi);
      scale(_hist_DeltaPhiSoft, invlumi);
      scale(_hist_DeltaS, invlumi);

      normalize(_hist_DeltaPtRelSoft_norm);
      normalize(_hist_DeltaPhiSoft_norm);
      normalize(_hist_DeltaS_norm);

      scale(_hist_LeadingHardJetPt, invlumi);
      scale(_hist_SubLeadingHardJetPt, invlumi);
      scale(_hist_LeadingSoftJetPt, invlumi);
      scale(_hist_SubLeadingSoftJetPt, invlumi);

      scale(_hist_LeadingHardJetEta, invlumi);
      scale(_hist_SubLeadingHardJetEta, invlumi);
      scale(_hist_LeadingSoftJetEta, invlumi);
      scale(_hist_SubLeadingSoftJetEta, invlumi);

    }

  };

  // This global object acts as a hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2013_I1273574);


}
