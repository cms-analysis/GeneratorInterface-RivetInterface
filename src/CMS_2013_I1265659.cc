// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

using namespace std;

namespace Rivet {


  class CMS_2013_I1265659 : public Analysis {

  public:

    /// Constructor
    CMS_2013_I1265659()
      : Analysis("CMS_2013_I1265659")
    {    }
    
    /// Book histograms and initialise projections before the run
    void init() {

      const FastJets jets(FinalState(-10, 10, 0.0*GeV), FastJets::ANTIKT, 0.5);
      addProjection(jets, "Jets");
      

      _h_hTotD = bookHistogram1D(1, 1, 1);
      _h_hTotDF = bookHistogram1D(1, 1, 2);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();
      
	
	const Jets& jets = applyProjection<FastJets>(event, "Jets").jetsByPt(30.0*GeV);
	
	if (jets.size()<3) vetoEvent;


        FourMomentum jet1 = jets[0].momentum();
        FourMomentum jet2 = jets[1].momentum();
        FourMomentum jet3 = jets[2].momentum();
       
       
       	double phi2 = jet2.phi();
	double phi3 = jet3.phi();

	double eta1 = jet1.eta();
	double eta2 = jet2.eta();
	double eta3 = jet3.eta();

	double pT1 = jet1.perp();
	
	double dEta23 = eta3 - eta2;	
	double dPhi23 = phi3 - phi2;       
        if (dPhi23 > M_PI)  dPhi23 -= 2*M_PI;  
        if (dPhi23 < -M_PI) dPhi23 += 2*M_PI;  
	double R23 = sqrt(dPhi23*dPhi23 + dEta23*dEta23);

	FourMomentum DiJet = jet1 + jet2;
	if ( fabs(eta1) > 2.5 || fabs(eta2) > 2.5 || DiJet.mass() < 220.0 || R23 < 0.5 || R23 > 1.5 || pT1 < 100) {
	  vetoEvent;
	}

	double beta = fabs(atan2(dPhi23,fabs(eta2)/eta2*dEta23));
	
	if (fabs(eta2) < 0.8) {
	  _h_hTotD->fill(beta, weight);
	}
	if (fabs(eta2) > 0.8) {
	  _h_hTotDF->fill(beta, weight);
	}

    }
    /// Normalise histograms etc., after the run
    void finalize() {
      
      double width = _h_hTotD->axis().binUpperEdge(0)-_h_hTotD->axis().binLowerEdge(0);

      normalize(_h_hTotD, width);
      normalize(_h_hTotDF, width);

    }
    

  private:
    /// @name Histograms

    AIDA::IHistogram1D* _h_hTotD;
    AIDA::IHistogram1D* _h_hTotDF;    
    
  };




  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2013_I1265659);

}
