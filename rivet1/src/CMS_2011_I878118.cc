// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

#include "HepMC/IO_AsciiParticles.h"
/// @todo Include more projections as required, e.g. ChargedFinalState, FastJets, ZFinder...

namespace Rivet {


  class CMS_2011_I878118 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CMS_2011_I878118()
      : Analysis("CMS_2011_I878118")
    {    }
      double Nevent ;

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      /// @todo Initialise and register projections here

      /// @todo Book histograms here, e.g.:
      // _h_XXXX = bookProfile1D(1, 1, 1);
      // _h_YYYY = bookHistogram1D(2, 1, 1);
      UnstableFinalState ufs(-6.0, 6.0, 0.0*GeV);
      addProjection(ufs, "UFS");
       Nevent = 0;
      _hpsi_prompt_y1= bookHistogram1D(8, 1, 1);
      _hpsi_prompt_y2= bookHistogram1D(9, 1, 1);
      _hpsi_prompt_y3= bookHistogram1D(10, 1, 1);

      _hpsi_non_prompt_y1= bookHistogram1D(11, 1, 1);
      _hpsi_non_prompt_y2= bookHistogram1D(12, 1, 1);
      _hpsi_non_prompt_y3= bookHistogram1D(13, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      /// @todo Do the event by event analysis here

      const UnstableFinalState& ufs = applyProjection<UnstableFinalState>(event, "UFS");
      bool nonPrompt = false;
      ParticleVector jpsi;
//      cout << " new event " << Nevent << endl;
      foreach (const Particle& p, ufs.particles()) {
          const PdgId id= abs(p.pdgId());
          if ( id == 443 ) {
          if ( p.hasAncestor(5) || p.hasAncestor(-5) ) { nonPrompt = true ;}
          double y=fabs(p.momentum().rapidity()) ;
          double pt = p.momentum().pT();
          
          if( y < 1.2 ) {
            if (nonPrompt) _hpsi_non_prompt_y1->fill(pt, weight);
            else if (!nonPrompt) _hpsi_prompt_y1->fill(pt, weight);
          }
          else if ( 1.2 < y and y < 1.6 ) {
            if (nonPrompt) _hpsi_non_prompt_y2->fill(pt, weight);
            else if (!nonPrompt) _hpsi_prompt_y2->fill(pt, weight);
          }
           
          else if ( 1.6 < y and y < 2.4 ) {
            if (nonPrompt) _hpsi_non_prompt_y3->fill(pt, weight);
            else if (!nonPrompt) _hpsi_prompt_y3->fill(pt, weight);
          }
          }
       }
/*
       if (nonPrompt){               
             const HepMC::GenEvent* pevt= &event.genEvent();
          //   print out the hepmc record in a human readbale format on file:IO_AsciiParticles.dat
             ascii_io.write_event (pevt);
       }
*/

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      /// @todo Normalise, scale and otherwise manipulate histograms here

      // scale(_h_YYYY, crossSection()/sumOfWeights()); # norm to cross section
      // normalize(_h_YYYY); # normalize to unity
      double norm = crossSection()/nanobarn/sumOfWeights()*0.0593;
      cout << " xsections for CMS_2011_I878118 " << endl;
      scale(_hpsi_prompt_y1      , norm);
      scale(_hpsi_prompt_y2      , norm);
      scale(_hpsi_prompt_y3      , norm);
      scale(_hpsi_non_prompt_y1      , norm);
      scale(_hpsi_non_prompt_y2      , norm);
      scale(_hpsi_non_prompt_y3      , norm);
    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here

  HepMC::IO_AsciiParticles ascii_io;


  private:

    /// @name Histograms
    //@{
//    AIDA::IProfile1D *_h_XXXX;
//    AIDA::IHistogram1D *_h_YYYY;
     AIDA::IHistogram1D *_hpsi_prompt_y1;
     AIDA::IHistogram1D *_hpsi_prompt_y2;
     AIDA::IHistogram1D *_hpsi_prompt_y3;
     AIDA::IHistogram1D *_hpsi_non_prompt_y1;
     AIDA::IHistogram1D *_hpsi_non_prompt_y2;
     AIDA::IHistogram1D *_hpsi_non_prompt_y3;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2011_I878118);

}
