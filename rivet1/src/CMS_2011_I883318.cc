// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
/// @todo Include more projections as required, e.g. ChargedFinalState, FastJets, ZFinder...

namespace Rivet {


  class CMS_2011_I883318 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CMS_2011_I883318()
      : Analysis("CMS_2011_I883318")
    {        
      setBeams(PROTON, PROTON);
      setNeedsCrossSection(true);
     }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      /// @todo Initialise and register projections here
      UnstableFinalState ufs(-8.0, 8.0, 0.0*GeV);
      addProjection(ufs, "UFS");

      /// @todo Book histograms here, e.g.:
      // _h_XXXX = bookProfile1D(1, 1, 1);
      // _h_YYYY = bookHistogram1D(2, 1, 1);
       _h_dsigdpt = bookHistogram1D(2, 1, 1);
       _h_dsigdy  = bookHistogram1D(3, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      /// @todo Do the event by event analysis here
      // Loop through unstable FS particles and look for B+ meson
       const UnstableFinalState& ufs = applyProjection<UnstableFinalState>(event, "UFS");

      /// @todo Do the event by event analysis here
//      cout << " analysing B+ events " <<endl;

       foreach (const Particle& p, ufs.particles()) {
          const PdgId id= abs(p.pdgId());
//          cout << " PdgId = " << id<< endl;
          if ( id == 521 ) {
              const double ptB= p.momentum().pT();
              const double pty= p.momentum().rapidity();
//              cout << " in loop  PdgId " << id  << endl;
//                  cout << " pt = " << ptB << " y = " << pty << endl;
               if (pty<2.4&&pty>-2.4){
             _h_dsigdpt->fill(ptB, weight);}
               if (ptB>5){
              _h_dsigdy ->fill(pty, weight);
              }
          }
       }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // value for xsection is in picobarn
      double invlumi = crossSection()/1000000/sumOfWeights();
      // apply factor of 2, since we count B+ and B-, but xsection is for B+
      scale(_h_dsigdpt, invlumi/2.0); 
      scale(_h_dsigdy, invlumi/2.0);

    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here


  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D *_h_dsigdpt;
    AIDA::IHistogram1D *_h_dsigdy;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2011_I883318);

}
