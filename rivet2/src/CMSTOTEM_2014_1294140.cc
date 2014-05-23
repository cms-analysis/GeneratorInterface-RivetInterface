// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  class CMSTOTEM_2014_1294140 : public Analysis {
  public:

    CMSTOTEM_2014_1294140()
      : Analysis("CMSTOTEM_2014_1294140") 
    {     }

    void init() {
      ChargedFinalState cfs(-7.0, 7.0, 0.0*GeV);
      addProjection(cfs, "CFS");

      _Nevt_after_cuts_or = 0;
      _Nevt_after_cuts_and = 0;
      _Nevt_after_cuts_orveto = 0;

      if(fuzzyEquals(sqrtS(), 8000*GeV, 1E-3)){
	_h_dNch_dEta_OR = bookHisto1D(1, 1, 1);
	_h_dNch_dEta_AND = bookHisto1D(2, 1, 1);
        _h_dNch_dEta_ORVETO = bookHisto1D(3, 1, 1);
      }
    }

    void analyze(const Event& event) {
      const double weight = event.weight();
      int count_plus = 0;
      int count_minus = 0;

      const ChargedFinalState& charged = applyProjection<ChargedFinalState>(event, "CFS");
      
      foreach (const Particle& p, charged.particles()) {
	if ( 5.3 < p.momentum().pseudorapidity() && p.momentum().pseudorapidity() < 6.5) count_plus++;
	if ( -6.5 < p.momentum().pseudorapidity() && p.momentum().pseudorapidity() < -5.3) count_minus++;
      }
   
      if ( count_plus > 0 || count_minus > 0 ) _Nevt_after_cuts_or += weight;    
      if ( count_plus > 0 && count_minus > 0 ) _Nevt_after_cuts_and += weight;
      if ( (count_plus > 0 && !count_minus > 0) || ( !count_plus > 0 && count_minus > 0 ) ) 
      		_Nevt_after_cuts_orveto += weight;
           
      foreach (const Particle& p, charged.particles()) {
        const double eta = p.momentum().eta();		
	if (count_plus > 0 || count_minus > 0) _h_dNch_dEta_OR->fill(fabs(eta), weight);
	if (count_plus > 0 && count_minus > 0) _h_dNch_dEta_AND->fill(fabs(eta), weight);
	if (( count_plus > 0 && !count_minus > 0 ) || ( !count_plus > 0 && count_minus > 0 )) _h_dNch_dEta_ORVETO->fill(fabs(eta), weight);	
      }
      
    }
    
    
    void finalize() {
      const double norm_or = 1.0/(2.*_Nevt_after_cuts_or);
      const double norm_and = 1.0/(2.*_Nevt_after_cuts_and);
      const double norm_orveto = 1.0/(2.*_Nevt_after_cuts_orveto);                                                                         
 
      scale(_h_dNch_dEta_OR, norm_or);
      scale(_h_dNch_dEta_AND, norm_and);
      scale(_h_dNch_dEta_ORVETO, norm_orveto);
    }


  private:

    Histo1DPtr _h_dNch_dEta_OR;
    Histo1DPtr _h_dNch_dEta_AND;
    Histo1DPtr _h_dNch_dEta_ORVETO;
      
    double _Nevt_after_cuts_or;
    double _Nevt_after_cuts_and;
    double _Nevt_after_cuts_orveto;

  };

  // This global object acts as a hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMSTOTEM_2014_1294140);

}

