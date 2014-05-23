// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/ParticleBase.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"


namespace Rivet {

  class CMS_EWK_11_015 : public Analysis {
  public:

    /// Constructor
    CMS_EWK_11_015() : Analysis("CMS_EWK_11_015")
    {      
      setNeedsCrossSection(true);
    }  

    /// Add projections and book histograms
    void init() {
               
      _sumW = 0;
      _sumW50 = 0;
      _sumWpT = 0;

      FinalState fs(-2.4, 2.4, 20.0*GeV);
      addProjection(fs, "FS");

      UnstableFinalState ufs(-2, 2, 15.0*GeV);
      addProjection(ufs, "UFS");
	
      ZFinder zfindermu(fs,-2.4, 2.4, 0.0*GeV, MUON, 
				81.0*GeV, 101.0*GeV, 0.1, false, true, 91.2*GeV);

      addProjection(zfindermu, "ZFinderMu");

      ZFinder zfinderel(fs, -2.4, 2.4, 0.0*GeV, ELECTRON, 
				81.0*GeV, 101.0*GeV, 0.1, false, true, 91.2*GeV);

      addProjection(zfinderel, "ZFinderEl");

      /// Histograms in nonboosted region of Z pT
      _h_dR_BB = bookHistogram1D(1, 1, 1);
      _h_dphi_BB = bookHistogram1D(2, 1, 1);
      _h_min_dR_ZB = bookHistogram1D(3, 1, 1);
      _h_A_ZBB = bookHistogram1D(4, 1, 1);
      
      /// Histograms in boosted region of Z pT (pT > 50 GeV)
      _h_dR_BB_boost = bookHistogram1D(5, 1, 1);
      _h_dphi_BB_boost = bookHistogram1D(6, 1, 1);
      _h_min_dR_ZB_boost = bookHistogram1D(7, 1, 1);
      _h_A_ZBB_boost = bookHistogram1D(8, 1, 1);    
      
      _h_min_ZpT = bookHistogram1D(9,1,1); 
            
    }
  

    /// Do the analysis
    void analyze(const Event& e) { 
        
      const double weight = e.weight();

      vector<FourMomentum> Bmom;
      
      const UnstableFinalState& ufs = applyProjection<UnstableFinalState>(e,"UFS");
      const ZFinder& zfindermu = applyProjection<ZFinder>(e, "ZFinderMu");
      const ZFinder& zfinderel = applyProjection<ZFinder>(e, "ZFinderEl");
     
      /// Look for a Z --> mu+ mu- event in the final state
      if( zfindermu.empty() && zfinderel.empty() ) vetoEvent;
      
      const ParticleVector& z = !zfindermu.empty() ? zfindermu.bosons() : zfinderel.bosons();
      const bool is_boosted = ( z[0].momentum().pT() > 50*GeV );           
    
      /// Loop over the unstable particles
      foreach ( const Particle& p, ufs.particles()) {
      
        const PdgId pid = p.pdgId();
       
        /// Look for particles with a bottom quark      
        if( PID::hasBottom(pid) ) {
          
          bool good_B = false;
	  const HepMC::GenParticle& pgen = p.genParticle();
	  HepMC::GenVertex* vgen = pgen.end_vertex();
	  
	  /// Loop over the decay products of each unstable particle.
	  /// Look for a couple of B hadrons.
	  for( HepMC::GenVertex::particles_out_const_iterator it = vgen->particles_out_const_begin(); 
	  				it !=  vgen->particles_out_const_end(); ++it ){
	    
	    /// If the particle produced has a bottom quark do not count it and go to the next loop cycle.	    
	    if( !( PID::hasBottom( (*it)->pdg_id() ) ) ){
	      good_B = true;
	      continue;
	    }
	    
	    else {
	      good_B = false;
	      break;
	    }		
	    					
	  }
	  
	  if( good_B ){
	    Bmom.push_back( p.momentum() );
	  }				
                    	
        }
        else continue;
      }


     
      /// If there are more than two B's in the final state veto the event   
      if( Bmom.size() != 2 ) {
        vetoEvent;
      }
      

      /// Calculate the observables values    
      double dphiBB = abs(Bmom.at(0).phi() - Bmom.at(1).phi());
      double dRBB = deltaR( Bmom.at(0), Bmom.at(1) );

      const FourMomentum& pZ = z[0].momentum();

      const bool closest_B = ( deltaR(pZ, Bmom.at(0)) < deltaR(pZ, Bmom.at(1)) );
      double mindR_ZB = closest_B ? deltaR(pZ, Bmom.at(0)) : deltaR(pZ, Bmom.at(1));
      double maxdR_ZB = closest_B ? deltaR(pZ, Bmom.at(1)) : deltaR(pZ, Bmom.at(0));
      
      double AZBB = ( maxdR_ZB - mindR_ZB ) / ( maxdR_ZB + mindR_ZB );
      
      /// Fill the histograms in the non-boosted region
      _h_dphi_BB->fill(dphiBB,weight);
      _h_dR_BB->fill(dRBB,weight);
      _h_min_dR_ZB->fill(mindR_ZB,weight);     
      _h_A_ZBB->fill(AZBB,weight);
      _sumW+=weight;
      _sumWpT+=weight;
     
      /// Fill the histograms in the boosted region   
      if( is_boosted ){
        _sumW50+=weight;
        _h_dphi_BB_boost->fill(dphiBB,weight);
        _h_dR_BB_boost->fill(dRBB,weight);
        _h_min_dR_ZB_boost->fill(mindR_ZB,weight);          
        _h_A_ZBB_boost->fill(AZBB,weight);
      }


      _h_min_ZpT->fill(0,weight);        
      
      if( pZ.pT() > 40*GeV ){
        _sumWpT+=weight;
        _h_min_ZpT->fill(40,weight);
      }
       if( pZ.pT() > 80*GeV ){
         _sumWpT+=weight; 
         _h_min_ZpT->fill(80,weight);
      }
      if( pZ.pT() > 120*GeV ){
        _sumWpT+=weight;
        _h_min_ZpT->fill(120,weight);
      }          
      
      Bmom.clear();
    }


    /// Finalize
    void finalize() {   
      
      std::cout << "Sum of weights = " << sumOfWeights() << std::endl;
      std::cout << "Total cross section = " << crossSection() << " pb " << std::endl;
      std::cout << "SumW = " << _sumW << std::endl;   
      std::cout << "SumW 50 = " << _sumW50 << std::endl;   
     
      /// Normalize excluding overflow bins 
      normalizeNoOverflows(_h_dR_BB, 0.7*crossSection()*_sumW/sumOfWeights());	/// d01-x01-y01   
      normalizeNoOverflows(_h_dphi_BB, 0.53*crossSection()*_sumW/sumOfWeights());	/// d02-x01-y01
      normalizeNoOverflows(_h_min_dR_ZB, 0.84*crossSection()*_sumW/sumOfWeights());	/// d03-x01-y01
      normalizeNoOverflows(_h_A_ZBB, 0.2*crossSection()*_sumW/sumOfWeights());	/// d04-x01-y01  
      
      normalizeNoOverflows(_h_dR_BB_boost, 0.84*crossSection()*_sumW50/sumOfWeights());	/// d05-x01-y01      
      normalizeNoOverflows(_h_dphi_BB_boost, 0.63*crossSection()*_sumW50/sumOfWeights());	/// d06-x01-y01
      normalizeNoOverflows(_h_min_dR_ZB_boost, 1*crossSection()*_sumW50/sumOfWeights());	/// d07-x01-y01
      normalizeNoOverflows(_h_A_ZBB_boost, 0.25*crossSection()*_sumW50/sumOfWeights());	/// d08-x01-y01
      
      normalizeNoOverflows(_h_min_ZpT, 40*crossSection()*_sumWpT/sumOfWeights());	/// d09-x01-y01      
          
    }
  private:
  
    double _sumW, _sumW50, _sumWpT;
    
    /// Histograms
    AIDA::IHistogram1D * _h_dphi_BB;
    AIDA::IHistogram1D * _h_dR_BB;
    AIDA::IHistogram1D * _h_min_dR_ZB;
    AIDA::IHistogram1D * _h_A_ZBB;
    AIDA::IHistogram1D * _h_dphi_BB_boost;
    AIDA::IHistogram1D * _h_dR_BB_boost;
    AIDA::IHistogram1D * _h_min_dR_ZB_boost;
    AIDA::IHistogram1D * _h_A_ZBB_boost;
    AIDA::IHistogram1D * _h_min_ZpT;    
    
    
    void normalizeNoOverflows(AIDA::IHistogram1D* plot, double integral){
      double factor=1.;
      if (plot->sumBinHeights()>0 && plot->sumAllBinHeights()>0){
        factor = plot->sumAllBinHeights()/plot->sumBinHeights();
      }
      normalize(plot, factor*integral);  
    }
    
  };

  /// This global object acts as a hook for the plugin system
  AnalysisBuilder<CMS_EWK_11_015> plugin_CMS_EWK_11_015;

}
