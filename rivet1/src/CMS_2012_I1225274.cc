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


  class CMS_2012_I1225274 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CMS_2012_I1225274()
      : Analysis("CMS_2012_I1225274")
    {    }
      double nmu_all ;
      double Nevent ;

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      /// @todo Initialise and register projections here
      UnstableFinalState ufs(-6.0, 6.0, 0.0*GeV);
      addProjection(ufs, "UFS");
      // _h_YYYY = bookHistogram1D(2, 1, 1);
      /// @todo Book histograms here, e.g.:

       _hist_pt_1S = bookHistogram1D(2,1,2);
       _hist_pt_2S = bookHistogram1D(3,1,2);
       _hist_pt_3S = bookHistogram1D(4,1,2);
       _hist_pt_1S_r = bookHistogram1D(2,1,2);
       _hist_pt_2S_r = bookHistogram1D(3,1,2);
       _hist_pt_3S_r = bookHistogram1D(4,1,2);
       _hist_pt_1S_ylt2 = bookHistogram1D(29,1,1);
       _hist_pt_2S_ylt2 = bookHistogram1D(30,1,1);
       _hist_pt_3S_ylt2 = bookHistogram1D(31,1,1);
// in bins of rapidity
       _hist_pt_1S_y1 = bookHistogram1D(5,1,2);
       _hist_pt_2S_y1 = bookHistogram1D(6,1,2);
       _hist_pt_3S_y1 = bookHistogram1D(7,1,2);
       _hist_pt_1S_y2 = bookHistogram1D(8,1,2);
       _hist_pt_2S_y2 = bookHistogram1D(9,1,2);
       _hist_pt_3S_y2 = bookHistogram1D(10,1,2);
       _hist_pt_1S_y3 = bookHistogram1D(11,1,2);
       _hist_pt_2S_y3 = bookHistogram1D(12,1,2);
       _hist_pt_3S_y3 = bookHistogram1D(13,1,2);
       _hist_pt_1S_y4 = bookHistogram1D(14,1,2);
       _hist_pt_2S_y4 = bookHistogram1D(15,1,2);
       _hist_pt_3S_y4 = bookHistogram1D(16,1,2);
       _hist_pt_1S_y5 = bookHistogram1D(17,1,2);
       _hist_pt_2S_y5 = bookHistogram1D(18,1,2);
       _hist_pt_3S_y5 = bookHistogram1D(19,1,2);
       _hist_pt_1S_y6 = bookHistogram1D(20,1,2);
       _hist_pt_2S_y6 = bookHistogram1D(21,1,2);
       _hist_pt_3S_y6 = bookHistogram1D(22,1,2);


       _hist_y_1S = bookHistogram1D(23,1,2);
       _hist_y_2S = bookHistogram1D(24,1,2);
       _hist_y_3S = bookHistogram1D(25,1,2);

       _hist_pt_1S_fid = bookHistogram1D(2,1,1);
       _hist_pt_2S_fid = bookHistogram1D(3,1,1);
       _hist_pt_3S_fid = bookHistogram1D(4,1,1);

       _hist_y_1S_fid = bookHistogram1D(23,1,1);
       _hist_y_2S_fid = bookHistogram1D(24,1,1);
       _hist_y_3S_fid = bookHistogram1D(25,1,1);
              
       nmu_all = 0;
       Nevent = 0;

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      /// @todo Do the event by event analysis here
      double nmu = 0 ;
      double nups = 0;
      double nups1 = 0;
      double nups2 = 0;
      double nups3 =0;
      double nup  =0;
      
      // Loop through unstable FS particles and look for Upsilon meson
//      cout << " new event "<< sqrtS() << endl;
      Nevent = Nevent + 1;
      const UnstableFinalState& ufs = applyProjection<UnstableFinalState>(event, "UFS");

      
      int ip =0;
      int nm_553 = 0 ;
      int nm_100553 = 0;
      int nm_200553 = 0;
      for (GenEvent::particle_const_iterator p = event.genEvent().particles_begin();
           p != event.genEvent().particles_end(); ++p) {
           ip++ ;
           int nm_test = 0;
           if ( fabs((*p)->pdg_id()) == 13 ) {
           double pt_muon = (*p)->momentum().perp();
           double eta_muon = (*p)->momentum().eta();
//           cout << " pt_muon "<< pt_muon << endl;
           if(fabs(eta_muon) < 0.8 && pt_muon > 3.75) { nm_test=nm_test + 1  ;}
           if(fabs(eta_muon) > 0.8 && fabs(eta_muon) < 1.6 && pt_muon > 3.5) { nm_test=nm_test + 1 ;}
           if(fabs(eta_muon) > 1.6 && fabs(eta_muon) < 2.4 && pt_muon > 3.0) { nm_test=nm_test + 1 ;}
           if(nm_test == 1 ) {
           if ( (*p)->production_vertex() ) {
               for ( HepMC::GenVertex::particle_iterator mother = (*p)->production_vertex()-> particles_begin(HepMC::parents);
               mother != (*p)->production_vertex()-> particles_end(HepMC::parents); ++mother ) {
//               cout << "\t"; (*mother)->print();
//               cout << "\t" << " mother " << (*mother)->pdg_id() << endl;
               if( (*mother)->pdg_id() == 553 ) { 
                  nm_553=nm_553 + 1 ;
/*                  if(nm_553 == 2) {
                  double pt= (*mother)->momentum().perp();
                  double eta = (*mother)->momentum().eta();
                  cout << " 553: pt " << pt << " eta " << eta << endl;
                  }
*/
               } 
               if( (*mother)->pdg_id() == 100553 ) { 
                   nm_100553= nm_100553 + 1 ;
/*                  if(nm_100553 == 2) {
                  double pt= (*mother)->momentum().perp();
                  double eta = (*mother)->momentum().eta();
                  cout << " 100553: pt " << pt << " eta " << eta << endl;
                  }
*/                  
               } 
               if( (*mother)->pdg_id() == 200553 ) { 
                  nm_200553 = nm_200553 + 1 ;
/*                  if(nm_200553 == 2) {
                  double pt= (*mother)->momentum().perp();
                  double eta = (*mother)->momentum().eta();
                  cout << " 200553: pt " << pt << " eta " << eta << endl;
                  }
*/                  
               } 
               }
               }
           }
           }
      }    



       foreach (const Particle& p, ufs.particles()) {
          const PdgId id= abs(p.pdgId());
//          cout << " PdgId = " << id<< endl;
             if ( id == 553 ) {
//              cout << " in loop  PdgId " << id  << endl;
              nup = nup +1 ; }
          const double pt = p.momentum().pT();
          const double y = fabs(p.momentum().rapidity());
           if( pt > 0. &&  pt < 50. && fabs(y) < 2.4 ) {
             if ( id == 553 ) {
              nups = nups +1 ;
              nups1 = nups1 + 1;
//              cout << " in loop  PdgId " << id  << endl;
//              cout << " pt = " << pt << " y = " << y << endl;
              _hist_pt_1S->fill(pt, weight);
              _hist_pt_1S_r->fill(pt, weight);
              if (y < 2.0)               {_hist_pt_1S_ylt2->fill(pt, weight);}
              if (y < 0.4)               {_hist_pt_1S_y1->fill(pt, weight);}
              else if (y>0.4 && y < 0.8) {_hist_pt_1S_y2->fill(pt, weight);} 
              else if (y>0.8 && y < 1.2) {_hist_pt_1S_y3->fill(pt, weight);} 
              else if (y>1.2 && y < 1.6) {_hist_pt_1S_y4->fill(pt, weight);} 
              else if (y>1.6 && y < 2.0) {_hist_pt_1S_y5->fill(pt, weight);} 
              else if (y>2.0 && y < 2.4) {_hist_pt_1S_y6->fill(pt, weight);} 
              _hist_y_1S ->fill(y, weight);
              if (nm_553 == 2  ) {
                _hist_pt_1S_fid->fill(pt, weight);
                _hist_y_1S_fid ->fill(y, weight);
//                cout << " muons from 553 found " << " pt = "<<pt << " y = " << y<< endl;
              }
             }
             else if (id == 100553 ) {
              nups = nups + 1 ;
              nups2 = nups2 + 1 ;
              _hist_pt_2S->fill(pt, weight);
              _hist_pt_2S_r->fill(pt, weight);
// check whether y cut < 2 or ycut < 1.2 ????              
              if (y < 2.0)               {_hist_pt_2S_ylt2->fill(pt, weight);}
              if (y < 0.4)               {_hist_pt_2S_y1->fill(pt, weight);}
              else if (y>0.4 && y < 0.8) {_hist_pt_2S_y2->fill(pt, weight);} 
              else if (y>0.8 && y < 1.2) {_hist_pt_2S_y3->fill(pt, weight);} 
              else if (y>1.2 && y < 1.6) {_hist_pt_2S_y4->fill(pt, weight);} 
              else if (y>1.6 && y < 2.0) {_hist_pt_2S_y5->fill(pt, weight);} 
              else if (y>2.0 && y < 2.4) {_hist_pt_2S_y6->fill(pt, weight);} 
              _hist_y_2S ->fill(y, weight);
              if (nm_100553 == 2 ) {
                _hist_pt_2S_fid->fill(pt, weight);
                _hist_y_2S_fid ->fill(y, weight);
//                cout << " muons from 100553 found " << " pt = "<<pt << " y = " << y<< endl;
              }
             } 
             else if (id == 200553 ) {
              nups = nups + 1;
              nups3 = nups3 + 1 ;
              _hist_pt_3S->fill(pt, weight);
              _hist_pt_3S_r->fill(pt, weight);
              if (y < 2.0)               {_hist_pt_3S_ylt2->fill(pt, weight);}
              if (y < 0.4)               {_hist_pt_3S_y1->fill(pt, weight);}
              else if (y>0.4 && y < 0.8) {_hist_pt_3S_y2->fill(pt, weight);} 
              else if (y>0.8 && y < 1.2) {_hist_pt_3S_y3->fill(pt, weight);} 
              else if (y>1.2 && y < 1.6) {_hist_pt_3S_y4->fill(pt, weight);} 
              else if (y>1.6 && y < 2.0) {_hist_pt_3S_y5->fill(pt, weight);} 
              else if (y>2.0 && y < 2.4) {_hist_pt_3S_y6->fill(pt, weight);} 
              _hist_y_3S ->fill(y, weight);
              if (nm_200553 == 2  ) {
                _hist_pt_3S_fid->fill(pt, weight);
                _hist_y_3S_fid ->fill(y, weight);
//                cout << " muons from 200553 found " << " pt = "<<pt << " y = " << y<< endl;
              }
             }
                     
          }
       }
       if ( nmu == 2 ){ nmu_all = nmu_all+weight ;}
       if( nups1 > 1 ) { cout << " FATAL nr upsilon1S  " << nups<< " " << nups1<< " " << nups2 << " " <<nups3  << endl; }
       if( nups2 > 1 ) { cout << " FATAL nr upsilon2S  " << nups<< " " << nups1<< " " << nups2 << " " <<nups3  << endl; }
       if( nups3 > 1 ) { cout << " FATAL nr upsilon3S  " << nups<< " " << nups1<< " " << nups2 << " " <<nups3  << endl; }
//       if( nup > 1 ) { cout << " all FATAL nr upsilon > 1 " << nup << endl; }
//       cout << " muons "<< nmu<<" " << nmu_all << endl;
/*       if (nmu == 2 ) {
           cout << " Nevent = "<<Nevent << endl;
          //   print out the hepmc record on screen
          //   event.genEvent().print() ;
          // need to create a pointer for using IO_ascii
          
   	       const HepMC::GenEvent* pevt= &event.genEvent();
          //   print out the hepmc record in a human readbale format on file:IO_AsciiParticles.dat
             ascii_io.write_event (pevt);
       }
*/       
 
    }


    /// Normalise histograms etc., after the run
    void finalize() {
       cout << " di-muons "<<  nmu_all << " fiducial muon efficency  " << nmu_all/sumOfWeights() <<endl;
      double invlumi = crossSection()/nanobarn/sumOfWeights();
      double BR1S = 0.0248 ;
      double BR2S = 0.0193 ;
      double BR3S = 0.0218 ;
      
       scale(_hist_pt_1S, invlumi*BR1S/4.8 ) ;
       scale(_hist_pt_2S, invlumi*BR2S/4.8 ) ;
       scale(_hist_pt_3S, invlumi*BR3S/4.8 ) ;
       scale(_hist_pt_1S_r, invlumi*BR1S/4.8 ) ;
       scale(_hist_pt_2S_r, invlumi*BR2S/4.8 ) ;
       scale(_hist_pt_3S_r, invlumi*BR3S/4.8 ) ;
       scale(_hist_pt_1S_ylt2, invlumi*BR1S/4. ) ;
       scale(_hist_pt_2S_ylt2, invlumi*BR2S/4. ) ;
       scale(_hist_pt_3S_ylt2, invlumi*BR3S/4. ) ;
       scale(_hist_pt_1S_y1, invlumi*BR1S/0.8 ) ;
       scale(_hist_pt_2S_y1, invlumi*BR2S/0.8 ) ;
       scale(_hist_pt_3S_y1, invlumi*BR3S/0.8 ) ;
       scale(_hist_pt_1S_y2, invlumi*BR1S/0.8 ) ;
       scale(_hist_pt_2S_y2, invlumi*BR2S/0.8 ) ;
       scale(_hist_pt_3S_y2, invlumi*BR3S/0.8 ) ;
       scale(_hist_pt_1S_y3, invlumi*BR1S/0.8 ) ;
       scale(_hist_pt_2S_y3, invlumi*BR2S/0.8 ) ;
       scale(_hist_pt_3S_y3, invlumi*BR3S/0.8 ) ;
       scale(_hist_pt_1S_y4, invlumi*BR1S/0.8 ) ;
       scale(_hist_pt_2S_y4, invlumi*BR2S/0.8 ) ;
       scale(_hist_pt_3S_y4, invlumi*BR3S/0.8 ) ;
       scale(_hist_pt_1S_y5, invlumi*BR1S/0.8 ) ;
       scale(_hist_pt_2S_y5, invlumi*BR2S/0.8 ) ;
       scale(_hist_pt_3S_y5, invlumi*BR3S/0.8 ) ;
       scale(_hist_pt_1S_y6, invlumi*BR1S/0.8 ) ;
       scale(_hist_pt_2S_y6, invlumi*BR2S/0.8 ) ;
       scale(_hist_pt_3S_y6, invlumi*BR3S/0.8 ) ;

       scale(_hist_y_1S, invlumi*BR1S/2. ) ;
       scale(_hist_y_2S, invlumi*BR2S/2. ) ;
       scale(_hist_y_3S, invlumi*BR3S/2. ) ;
       
       scale(_hist_pt_1S_fid, invlumi/4.8 ) ;
       scale(_hist_pt_2S_fid, invlumi/4.8 ) ;
       scale(_hist_pt_3S_fid, invlumi/4.8 ) ;

       scale(_hist_y_1S_fid, invlumi/2. ) ;
       scale(_hist_y_2S_fid, invlumi/2. ) ;
       scale(_hist_y_3S_fid, invlumi/2. ) ;


    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here


  private:

    /// @name Histograms
    //@{
       AIDA::IHistogram1D *_hist_pt_1S ;
       AIDA::IHistogram1D *_hist_pt_2S ;
       AIDA::IHistogram1D *_hist_pt_3S ;
       AIDA::IHistogram1D *_hist_pt_1S_r ;
       AIDA::IHistogram1D *_hist_pt_2S_r ;
       AIDA::IHistogram1D *_hist_pt_3S_r ;
       AIDA::IHistogram1D *_hist_pt_1S_ylt2;
       AIDA::IHistogram1D *_hist_pt_2S_ylt2;
       AIDA::IHistogram1D *_hist_pt_3S_ylt2;
       AIDA::IHistogram1D *_hist_pt_1S_y1;
       AIDA::IHistogram1D *_hist_pt_2S_y1;
       AIDA::IHistogram1D *_hist_pt_3S_y1;
       AIDA::IHistogram1D *_hist_pt_1S_y2;
       AIDA::IHistogram1D *_hist_pt_2S_y2;
       AIDA::IHistogram1D *_hist_pt_3S_y2;
       AIDA::IHistogram1D *_hist_pt_1S_y3;
       AIDA::IHistogram1D *_hist_pt_2S_y3;
       AIDA::IHistogram1D *_hist_pt_3S_y3;
       AIDA::IHistogram1D *_hist_pt_1S_y4;
       AIDA::IHistogram1D *_hist_pt_2S_y4;
       AIDA::IHistogram1D *_hist_pt_3S_y4;
       AIDA::IHistogram1D *_hist_pt_1S_y5;
       AIDA::IHistogram1D *_hist_pt_2S_y5;
       AIDA::IHistogram1D *_hist_pt_3S_y5;
       AIDA::IHistogram1D *_hist_pt_1S_y6;
       AIDA::IHistogram1D *_hist_pt_2S_y6;
       AIDA::IHistogram1D *_hist_pt_3S_y6;

       AIDA::IHistogram1D *_hist_y_1S ;
       AIDA::IHistogram1D *_hist_y_2S ;
       AIDA::IHistogram1D *_hist_y_3S ;
       
       AIDA::IHistogram1D *_hist_pt_1S_fid ;
       AIDA::IHistogram1D *_hist_pt_2S_fid ;
       AIDA::IHistogram1D *_hist_pt_3S_fid ;

       AIDA::IHistogram1D *_hist_y_1S_fid ;
       AIDA::IHistogram1D *_hist_y_2S_fid ;
       AIDA::IHistogram1D *_hist_y_3S_fid ;

    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2012_I1225274);

}
