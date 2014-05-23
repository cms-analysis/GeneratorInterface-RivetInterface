// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  class CMS_2013_I1261026 : public Analysis {

  public:
  
    CMS_2013_I1261026()
      : Analysis("CMS_2013_I1261026")
    {    }


  public:

    void init() {

      // gives the range of eta and min pT for the final state from which I get the jets
        FastJets jetpro (ChargedFinalState(-2.4, 2.4, 0.25*GeV), FastJets::ANTIKT, 0.5);
        addProjection(jetpro, "Jets");

	const ChargedFinalState cfs(-2.4, 2.4, 0.25*GeV);   
        addProjection(cfs, "CFS250");

	//for MinBias trigger
	const ChargedFinalState cfsBSCplus(3.23, 4.65, 500*MeV);
        addProjection(cfsBSCplus, "cfsBSCplus");

	const ChargedFinalState cfsBSCminus(-4.65, -3.23, 500*MeV);
       	addProjection(cfsBSCminus, "cfsBSCminus");	


	_h_AllTrkMeanPt	           = bookProfile1D(1, 1, 1);
	_h_SoftTrkMeanPt           = bookProfile1D(2, 1, 1);
	_h_IntrajetTrkMeanPt       = bookProfile1D(3, 1, 1);
	_h_IntrajetLeaderTrkMeanPt = bookProfile1D(4, 1, 1);
	_h_MeanJetPt		   = bookProfile1D(5, 1, 1);
	_h_JetRate5GeV 	   	   = bookProfile1D(6, 1, 1);
	_h_JetRate30GeV	  	   = bookProfile1D(7, 1, 1);

	for (int ihist = 0; ihist < 5; ++ihist){
	  _h_JetSpectrum[ihist] = bookHisto1D(ihist+8, 1, 1);
  	  _h_JetStruct[ihist] = bookHisto1D(ihist+13, 1, 1);
	  
	  //temp histograms distribution parameters and sequent SEM calculation
  	  _AllTrkSpectrum[ihist]        = Histo1D(200, 0.0, 20.0);
	  _SoftTrkSpectrum[ihist]        = Histo1D(100, 0.0, 15.0);
	  _JetTrkSpectrum[ihist]        = Histo1D(100, 0.0, 20.0);
	  _JetLTrkSpectrum[ihist]        = Histo1D(100, 0.0, 20.0);
	  	  
	}
	
	MultBinCent[0]=20;
	MultBinCent[1]=40;
	MultBinCent[2]=65;
	MultBinCent[3]=95;
	MultBinCent[4]=125;
 
    }

	 
    /// Perform the per-event analysis
     void analyze(const Event& event) {
	const double weight = event.weight();
	
	//MinBias trigger
	const ChargedFinalState& cfsBSCplus = applyProjection<ChargedFinalState>(event, "cfsBSCplus");
	if (cfsBSCplus.empty()) vetoEvent;
	const ChargedFinalState& cfsBSCminus = applyProjection<ChargedFinalState>(event, "cfsBSCminus");
	if (cfsBSCminus.empty()) vetoEvent;
	  

	const ChargedFinalState& cfsp = applyProjection<ChargedFinalState>(event, "CFS250");
	if (cfsp.empty()) vetoEvent;
	 
	const FastJets& jetpro = applyProjection<FastJets>(event, "Jets");
      	const Jets& jets = jetpro.jetsByPt(5.0*GeV);

	const int _mult = cfsp.size();
	 
	int multbin[6] = {10, 30, 50, 80, 110, 140}; 
	for (int ibin = 0; ibin < 5; ++ibin){
	  if(_mult > multbin[ibin] && _mult <= multbin[ibin + 1]){
	  
	    passedEv[ibin]++;
	    EventDecomp(event, _h_JetStruct[ibin], &JetStructNorm[ibin], &_AllTrkSpectrum[ibin],
	    &_SoftTrkSpectrum[ibin], &_JetTrkSpectrum[ibin], &_JetLTrkSpectrum[ibin], weight);

	    for(signed int ijets = 0; ijets < (int)jets.size(); ijets++){
	      if(fabs(jets[ijets].momentum().eta()) < 1.9){
	        _h_JetSpectrum[ibin]->fill(jets[ijets].momentum().pT()/GeV, weight);
	        if(jets[ijets].momentum().pT() > 5*GeV) JetCounter5GeV[ibin] += weight;
	        if(jets[ijets].momentum().pT() > 30*GeV) JetCounter30GeV[ibin] += weight;
	      } 	
	    }	 
	  }
	}

     }

    /// Normalise histograms etc., after the run
    void finalize() {
	
	double SEM;
	
	for (int i = 0; i < 5; ++i){
	   
	  //all trk mean pT vs Nch
	  _h_AllTrkMeanPt->fill(MultBinCent[i], _AllTrkSpectrum[i].mean(), GetMeanError(_AllTrkSpectrum[i]));

	  //soft trk mean pT vs Nch
	  _h_SoftTrkMeanPt->fill(MultBinCent[i], _SoftTrkSpectrum[i].mean(), GetMeanError(_SoftTrkSpectrum[i]));
	 
	  //intrjet trk mean pT vs Nch
	  _h_IntrajetTrkMeanPt->fill(MultBinCent[i], _JetTrkSpectrum[i].mean(), GetMeanError(_JetTrkSpectrum[i]));

	  //intrjet leader trk mean pT vs Nch
	   _h_IntrajetLeaderTrkMeanPt->fill(MultBinCent[i], _JetLTrkSpectrum[i].mean(), GetMeanError(_JetLTrkSpectrum[i]));
	
	  //jet mean pT vs Nch
	   SEM = (_h_JetSpectrum[i] -> stdDev())/(sqrt(_h_JetSpectrum[i] -> sumW())) / _h_JetSpectrum[i]->mean();
	  _h_MeanJetPt->fill(MultBinCent[i], _h_JetSpectrum[i]->mean(), SEM);
	  
	  //jet rates
	  AvJetRate5[i] = JetCounter5GeV[i] / passedEv[i];
	  AvJetRate30[i] = JetCounter30GeV[i] / passedEv[i];
		 
	  if(JetCounter5GeV[i] != 0) SEM = 1 / sqrt(JetCounter5GeV[i]);
	  else SEM=0;
	  _h_JetRate5GeV->fill(MultBinCent[i], AvJetRate5[i], SEM);
		
	  if(JetCounter30GeV[i] != 0) SEM = 1 / sqrt(JetCounter30GeV[i]);
	  else SEM=0;
	  _h_JetRate30GeV->fill(MultBinCent[i], AvJetRate30[i], SEM);
	  
          scale(_h_JetSpectrum[i], 4.0 / JetCounter5GeV[i]);	 
          scale(_h_JetStruct[i], 0.08 / JetStructNorm[i]);
	 }
	 	 
    }

      double GetMeanError(const Histo1D& hist){
	double SEM = hist.stdDev() / sqrt(hist.numEntries()); //Standard error of the mean
	return SEM / hist.mean(); // relative SEM
      }
 
      void EventDecomp(const Event& event, Histo1DPtr JetSruct, double* JetStructNorm, Histo1D* AllTrk, Histo1D* SoftTrk, Histo1D* JetTrk, Histo1D* _JetLTrk, const double weight){

	 struct TrkInJet{double pt; double eta; double phi; double R;}; 
	 TrkInJet JetConstituents[100][100];//1-st index - the number of the jet, 2-nd index - track in the jet
	 TrkInJet JetsEv[100];
	 int j[100];
	 int jCount=0;

	 for(int i = 0; i < 100; i++){
	    j[i]=0;
	    JetsEv[i].pt=0;
            JetsEv[i].eta=0;
            JetsEv[i].phi=0;
            for(int k=0; k < 100; k++){
		JetConstituents[i][k].pt=0;
		JetConstituents[i][k].phi=0;
		JetConstituents[i][k].eta=0;
		JetConstituents[i][k].R=0;
	    }
         }

	 const FastJets& jetpro = applyProjection<FastJets>(event, "Jets");
         const Jets& jets = jetpro.jetsByPt(5.0*GeV);

	//-------Event Decomposiotn--Begining----------------------------------

	 for(signed int ijets = 0; ijets < (int)jets.size(); ijets++){
	   JetsEv[ijets].pt = jets[ijets].momentum().pT();
	   JetsEv[ijets].eta = jets[ijets].momentum().eta();
	   JetsEv[ijets].phi = jets[ijets].momentum().phi();
	   jCount++;
	 }
	  
	 const ChargedFinalState& cfsp = applyProjection<ChargedFinalState>(event, "CFS250");
	 foreach (const Particle& p, cfsp.particles()) {
	 	   
	  AllTrk -> fill(p.momentum().pT()/GeV, weight); 
	  int flag = 0;		
	  for(int i = 0; i < jCount; i++){
	     const double delta_phi = deltaPhi(JetsEv[i].phi, p.momentum().phi());
	     const double delta_eta = JetsEv[i].eta - p.momentum().eta();
	     const double R = sqrt(delta_phi * delta_phi + delta_eta * delta_eta); 
	     if(R <= 0.5){
		flag++;
		JetConstituents[i][j[i]].pt = p.momentum().pT();	
		JetConstituents[i][j[i]].R = R;
		j[i]++;
	     }
	   }
	
	   if(flag == 0) SoftTrk -> fill(p.momentum().pT(), weight);
	   
	 }
	
	 for(int i = 0; i < jCount; i++){
  	   double PtInjetLeader = 0;
	   if(JetsEv[i].eta > -1.9 && JetsEv[i].eta < 1.9){// only fully reconstructed jets for internal jet studies
	     for(int k = 0; k < j[i]; k++){
			
	       JetTrk -> fill(JetConstituents[i][k].pt * weight);
	       JetSruct -> fill(JetConstituents[i][k].R, JetConstituents[i][k].pt/JetsEv[i].pt);
	       *JetStructNorm += JetConstituents[i][k].pt / JetsEv[i].pt;
	       if(PtInjetLeader < JetConstituents[i][k].pt) PtInjetLeader = JetConstituents[i][k].pt;
	     }
             if(PtInjetLeader != 0)_JetLTrk -> fill(PtInjetLeader, weight);
	  }
	}
	
    } 

  private:
	
	double JetStructNorm[5];
	double MultBinCent[5]; 
    	double AvJetRate5[5], AvJetRate30[5];

	// counters
	double JetCounter5GeV[5], JetCounter30GeV[5], passedEv[5];

	Profile1DPtr _h_AllTrkMeanPt;
	Profile1DPtr _h_SoftTrkMeanPt;
	Profile1DPtr _h_IntrajetTrkMeanPt;
	Profile1DPtr _h_IntrajetLeaderTrkMeanPt;
	Profile1DPtr _h_MeanJetPt;
	Profile1DPtr _h_JetRate5GeV;
	Profile1DPtr _h_JetRate30GeV;
	
	Histo1DPtr _h_JetSpectrum[5];	
	Histo1DPtr _h_JetStruct[5];

	//temp histograms
	Histo1D _AllTrkSpectrum[5];
	Histo1D _SoftTrkSpectrum[5];
	Histo1D _JetTrkSpectrum[5];
	Histo1D _JetLTrkSpectrum[5];
	
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2013_I1261026);

}
