// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/VetoedFinalState.hh"


/// @todo Include more projections as required, e.g. ChargedFinalState, FastJets, ZFinder...

namespace Rivet {


  class CMS_2013_I1261026 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CMS_2013_I1261026()
      : Analysis("CMS_2013_I1261026")
    {    }

    //@}


  public:

    void init() {

      // gives the range of eta and min pT for the final state from which I get the jets
        FastJets jetpro (ChargedFinalState(-2.4, 2.4, 0.25*GeV), FastJets::ANTIKT, 0.5);
        addProjection(jetpro, "Jets");

	const ChargedFinalState cfs(-2.4, 2.4, 0.25*GeV);   
        addProjection(cfs, "CFS250");

	//for MinBias trigger---------------------------------------------
	const ChargedFinalState cfsBSCplus(3.23, 4.65, 500*MeV);
        addProjection(cfsBSCplus, "cfsBSCplus");

	const ChargedFinalState cfsBSCminus(-4.65, -3.23, 500*MeV);
       	addProjection(cfsBSCminus, "cfsBSCminus");
	

	//----------------------------------------------------------------------
	_h_AllTrkMeanPt            = bookProfile1D(1, 1, 1);
	_h_SoftTrkMeanPt           = bookProfile1D(2, 1, 1);
	_h_IntrajetTrkMeanPt       = bookProfile1D(3, 1, 1);
	_h_IntrajetLeaderTrkMeanPt = bookProfile1D(4, 1, 1);
	_h_MeanJetPt		   = bookProfile1D(5, 1, 1);
	_h_JetRate5GeV 	   	   = bookProfile1D(6, 1, 1);
	_h_JetRate30GeV	  	   = bookProfile1D(7, 1, 1);

	_h_JetSpectrum1 = bookHistogram1D(8, 1, 1);
	_h_JetSpectrum2 = bookHistogram1D(9, 1, 1);
	_h_JetSpectrum3 = bookHistogram1D(10, 1, 1);
	_h_JetSpectrum4 = bookHistogram1D(11, 1, 1);
	_h_JetSpectrum5 = bookHistogram1D(12, 1, 1);

	
	_h_JetStruct1 = bookHistogram1D(13, 1, 1);
	_h_JetStruct2 = bookHistogram1D(14, 1, 1);
	_h_JetStruct3 = bookHistogram1D(15, 1, 1);
	_h_JetStruct4 = bookHistogram1D(16, 1, 1);
	_h_JetStruct5 = bookHistogram1D(17, 1, 1);
	
	
	//-----------------------------------------------------------------------------------------------
	//histograms that are need for computation of distribution parameter and sequent SEM calculation
	_AllTrkSpectrum1        = bookHistogram1D("AllTrkSpectrum1", 200, 0.0, 20.0);
	_AllTrkSpectrum2        = bookHistogram1D("AllTrkSpectrum2", 200, 0.0, 20.0);
	_AllTrkSpectrum3        = bookHistogram1D("AllTrkSpectrum3", 200, 0.0, 20.0);
	_AllTrkSpectrum4        = bookHistogram1D("AllTrkSpectrum4", 200, 0.0, 20.0);
	_AllTrkSpectrum5        = bookHistogram1D("AllTrkSpectrum5", 200, 0.0, 20.0);

	_SoftTrkSpectrum1        = bookHistogram1D("SoftTrkSpectrum1", 100, 0.0, 15.0);
	_SoftTrkSpectrum2        = bookHistogram1D("SoftTrkSpectrum2", 100, 0.0, 15.0);
	_SoftTrkSpectrum3        = bookHistogram1D("SoftTrkSpectrum3", 100, 0.0, 15.0);
	_SoftTrkSpectrum4        = bookHistogram1D("SoftTrkSpectrum4", 100, 0.0, 15.0);
	_SoftTrkSpectrum5        = bookHistogram1D("SoftTrkSpectrum5", 100, 0.0, 15.0);

	_JetTrkSpectrum1        = bookHistogram1D("JetTrkSpectrum1", 100, 0.0, 20.0);
	_JetTrkSpectrum2        = bookHistogram1D("JetTrkSpectrum2", 100, 0.0, 20.0);
	_JetTrkSpectrum3        = bookHistogram1D("JetTrkSpectrum3", 100, 0.0, 20.0);
	_JetTrkSpectrum4        = bookHistogram1D("JetTrkSpectrum4", 100, 0.0, 20.0);
	_JetTrkSpectrum5        = bookHistogram1D("JetTrkSpectrum5", 100, 0.0, 20.0);

	_JetLTrkSpectrum1        = bookHistogram1D("JetLTrkSpectrum1", 100, 0.0, 20.0);
	_JetLTrkSpectrum2        = bookHistogram1D("JetLTrkSpectrum2", 100, 0.0, 20.0);
	_JetLTrkSpectrum3        = bookHistogram1D("JetLTrkSpectrum3", 100, 0.0, 20.0);
	_JetLTrkSpectrum4        = bookHistogram1D("JetLTrkSpectrum4", 100, 0.0, 20.0);
	_JetLTrkSpectrum5        = bookHistogram1D("JetLTrkSpectrum5", 100, 0.0, 20.0);

	//----------------------------------------------------------------------------------------------------
	MultBinCent[0]=20;
	MultBinCent[1]=40;
	MultBinCent[2]=65;
	MultBinCent[3]=95;
	MultBinCent[4]=125;
    }

	 

	 void EvetnDecompostion(const Event& event, IHistogram1D* JetSruct, double* JetStructNorm, IHistogram1D * AllTrk, IHistogram1D * SoftTrk, IHistogram1D * JetTrk, IHistogram1D * _JetLTrk, const double weight){

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
            for(int k=0; k<100; k++){
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
	 	   
	   AllTrk->fill(p.momentum().pT()/GeV, weight);
      
	   int flag = 0;
		

	  for(int i = 0;i < jCount;i++){
	     const double delta_phi = deltaPhi(JetsEv[i].phi, p.momentum().phi());
	     const double delta_eta = JetsEv[i].eta-p.momentum().eta();
	     const double R = sqrt(delta_phi*delta_phi + delta_eta*delta_eta); 
	     if(R <= 0.5){
		flag++;
		JetConstituents[i][j[i]].pt=p.momentum().pT();	
		JetConstituents[i][j[i]].R=R;
		j[i]++;
	     }
	   }
	
	   if(flag==0) SoftTrk->fill(p.momentum().pT(), weight);
	   
	}
	
	
	for(int i = 0; i < jCount; i++){
		 double PtInjetLeader=0;
		if(JetsEv[i].eta > -1.9 && JetsEv[i].eta < 1.9){// only fully reconstructed jets for internal jet studies
		for(int k=0; k<j[i]; k++){
			
			JetTrk->fill(JetConstituents[i][k].pt*weight);
			JetSruct->fill(JetConstituents[i][k].R, JetConstituents[i][k].pt/JetsEv[i].pt);
			*JetStructNorm += JetConstituents[i][k].pt/JetsEv[i].pt;
			if(PtInjetLeader<JetConstituents[i][k].pt)PtInjetLeader=JetConstituents[i][k].pt;
		}
               if(PtInjetLeader!=0)_JetLTrk->fill(PtInjetLeader, weight);
		}

	}
	

	  //------Event Decomposiotn--End-----------------------------------
	  //---------------------------------------------------------------------

	 }

    /// Perform the per-event analysis
     void analyze(const Event& event) {
      const double weight = event.weight();

      /// @todo Do the event by event analysis here

	  //MinBias trigger-----------------------------------------------------------------------
	  const ChargedFinalState& cfsBSCplus = applyProjection<ChargedFinalState>(event, "cfsBSCplus");
	  if (cfsBSCplus.empty()) vetoEvent;
	  const ChargedFinalState& cfsBSCminus = applyProjection<ChargedFinalState>(event, "cfsBSCminus");
	  if (cfsBSCminus.empty()) vetoEvent;
	  //--------------------------------------------------------------------------------------


	const ChargedFinalState& cfsp = applyProjection<ChargedFinalState>(event, "CFS250");
	if (cfsp.empty()) vetoEvent;
	 
	const FastJets& jetpro = applyProjection<FastJets>(event, "Jets");
      	const Jets& jets = jetpro.jetsByPt(5.0*GeV);

	const int _mult = cfsp.size();

	 


	  if(_mult>10 && _mult<=30){
	  
	  passedEv[0]++;
	  EvetnDecompostion(event, _h_JetStruct1, &JetStructNorm[0], _AllTrkSpectrum1, _SoftTrkSpectrum1, _JetTrkSpectrum1, _JetLTrkSpectrum1, weight);

	  for(signed int ijets = 0; ijets < (int)jets.size(); ijets++){
		  if(std::fabs(jets[ijets].momentum().eta())<1.9){
		  _h_JetSpectrum1->fill(jets[ijets].momentum().pT()/GeV, weight);
		  if(jets[ijets].momentum().pT()>5)JetCounter5GeV[0]+=weight;
		  if(jets[ijets].momentum().pT()>30)JetCounter30GeV[0]+=weight;
		  } 	
	  }	 
	  }


	  if(_mult>30 && _mult<=50){

	  passedEv[1]++;
	  EvetnDecompostion(event, _h_JetStruct2, &JetStructNorm[1], _AllTrkSpectrum2, _SoftTrkSpectrum2, _JetTrkSpectrum2, _JetLTrkSpectrum2, weight);

	  for(signed int ijets = 0; ijets < (int)jets.size(); ijets++){
		  if(std::fabs(jets[ijets].momentum().eta())<1.9){
		  _h_JetSpectrum2->fill(jets[ijets].momentum().pT()/GeV, weight);
		  if(jets[ijets].momentum().pT()>5)JetCounter5GeV[1]+=weight;
		  if(jets[ijets].momentum().pT()>30)JetCounter30GeV[1]+=weight;
		  }
	 }
	  }


	  if(_mult>50 && _mult<=80){

	  passedEv[2]++;
	  EvetnDecompostion(event, _h_JetStruct3, &JetStructNorm[2], _AllTrkSpectrum3, _SoftTrkSpectrum3, _JetTrkSpectrum3, _JetLTrkSpectrum3, weight);

	  for(signed int ijets = 0; ijets < (int)jets.size(); ijets++){
		  if(std::fabs(jets[ijets].momentum().eta())<1.9){
		  _h_JetSpectrum3->fill(jets[ijets].momentum().pT()/GeV, weight);
		  if(jets[ijets].momentum().pT()>5)JetCounter5GeV[2]+=weight;
		  if(jets[ijets].momentum().pT()>30)JetCounter30GeV[2]+=weight;
		  }
	 }
	  }

	  if(_mult>80 && _mult<=110){

	  passedEv[3]++;
	  EvetnDecompostion(event, _h_JetStruct4, &JetStructNorm[3], _AllTrkSpectrum4, _SoftTrkSpectrum4, _JetTrkSpectrum4, _JetLTrkSpectrum4, weight);

	  for(signed int ijets = 0; ijets < (int)jets.size(); ijets++){
		  if(std::fabs(jets[ijets].momentum().eta())<1.9){
		  _h_JetSpectrum4->fill(jets[ijets].momentum().pT()/GeV, weight);
		  if(jets[ijets].momentum().pT()>5)JetCounter5GeV[3]+=weight;
		  if(jets[ijets].momentum().pT()>30)JetCounter30GeV[3]+=weight;
		  } 
	 }
	  }


	  if(_mult>110 && _mult<=140){

	  passedEv[4]++;
	  EvetnDecompostion(event, _h_JetStruct5, &JetStructNorm[4], _AllTrkSpectrum5, _SoftTrkSpectrum5, _JetTrkSpectrum5, _JetLTrkSpectrum5, weight);

	  for(signed int ijets = 0; ijets < (int)jets.size(); ijets++){
		  if(std::fabs(jets[ijets].momentum().eta())<1.9){
		  _h_JetSpectrum5->fill(jets[ijets].momentum().pT()/GeV, weight);
		  if(jets[ijets].momentum().pT()>5)JetCounter5GeV[4]+=weight;
		  if(jets[ijets].momentum().pT()>30)JetCounter30GeV[4]+=weight;
		  } 
	 }
	  }


	}//end of analyze()

	double GetMeanError(IHistogram1D* hist){
	  double hist_rms = hist->rms();
	  double mean_hist= hist->mean();
	  double Nenetries = hist->entries();
	  double SEM = hist_rms/sqrt(Nenetries); //Standard error of the mean
	  SEM = SEM/mean_hist; // relative SEM
	  return SEM;
	}

    /// Normalise histograms etc., after the run
    void finalize() {

	//all trk mean pT vs Nch------------------------------------------------------
	 Mean[0] =  _AllTrkSpectrum1->mean();
	 Mean[1] =  _AllTrkSpectrum2->mean();
	 Mean[2] =  _AllTrkSpectrum3->mean();
	 Mean[3] =  _AllTrkSpectrum4->mean();
	 Mean[4] =  _AllTrkSpectrum5->mean();
	 SEM[0] = GetMeanError(_AllTrkSpectrum1);
	 SEM[1] = GetMeanError(_AllTrkSpectrum2);
	 SEM[2] = GetMeanError(_AllTrkSpectrum3);
	 SEM[3] = GetMeanError(_AllTrkSpectrum4);
	 SEM[4] = GetMeanError(_AllTrkSpectrum5);
	 for(int i=0; i<5; i++){
	   _h_AllTrkMeanPt->fill(MultBinCent[i], Mean[i], SEM[i]);
	 }
	//---------------------------------------------------------------------------

	  //soft trk mean pT vs Nch--------------------------------------------------
	 Mean[0] =  _SoftTrkSpectrum1->mean();
	 Mean[1] =  _SoftTrkSpectrum2->mean();
	 Mean[2] =  _SoftTrkSpectrum3->mean();
	 Mean[3] =  _SoftTrkSpectrum4->mean();
	 Mean[4] =  _SoftTrkSpectrum5->mean();
	 SEM[0] = GetMeanError(_SoftTrkSpectrum1);
	 SEM[1] = GetMeanError(_SoftTrkSpectrum2);
	 SEM[2] = GetMeanError(_SoftTrkSpectrum3);
	 SEM[3] = GetMeanError(_SoftTrkSpectrum4);
	 SEM[4] = GetMeanError(_SoftTrkSpectrum5);
	 for(int i=0; i<5; i++){
	   _h_SoftTrkMeanPt->fill(MultBinCent[i], Mean[i], SEM[i]);
	 }
	//----------------------------------------------------------------



	//intrjet trk mean pT vs Nch--------------------------------------------------
	 Mean[0] =  _JetTrkSpectrum1->mean();
	 Mean[1] =  _JetTrkSpectrum2->mean();
	 Mean[2] =  _JetTrkSpectrum3->mean();
	 Mean[3] =  _JetTrkSpectrum4->mean();
	 Mean[4] =  _JetTrkSpectrum5->mean();
	 SEM[0] = GetMeanError(_JetTrkSpectrum1);
	 SEM[1] = GetMeanError(_JetTrkSpectrum2);
	 SEM[2] = GetMeanError(_JetTrkSpectrum3);
	 SEM[3] = GetMeanError(_JetTrkSpectrum4);
	 SEM[4] = GetMeanError(_JetTrkSpectrum5);
 	 for(int i=0; i<5; i++){
	   _h_IntrajetTrkMeanPt->fill(MultBinCent[i], Mean[i], SEM[i]);
	 }
	//----------------------------------------------------------------


	//intrjet leader trk mean pT vs Nch--------------------------------------------------
	 Mean[0] = _JetLTrkSpectrum1->mean();
	 Mean[1] = _JetLTrkSpectrum2->mean();
	 Mean[2] = _JetLTrkSpectrum3->mean();
	 Mean[3] = _JetLTrkSpectrum4->mean();
	 Mean[4] = _JetLTrkSpectrum5->mean();
	 SEM[0] = GetMeanError(_JetLTrkSpectrum1);
	 SEM[1] = GetMeanError(_JetLTrkSpectrum2);
	 SEM[2] = GetMeanError(_JetLTrkSpectrum3);
	 SEM[3] = GetMeanError(_JetLTrkSpectrum4);
	 SEM[4] = GetMeanError(_JetLTrkSpectrum5);
	 for(int i=0; i<5; i++){
	   _h_IntrajetLeaderTrkMeanPt->fill(MultBinCent[i], Mean[i], SEM[i]);
	 }
	//----------------------------------------------------------------


	//jet mean pT vs Nch--------------------------------------------------
	 Mean[0] = _h_JetSpectrum1->mean();
	 Mean[1] = _h_JetSpectrum2->mean();
	 Mean[2] = _h_JetSpectrum3->mean();
	 Mean[3] = _h_JetSpectrum4->mean();
	 Mean[4] = _h_JetSpectrum5->mean();
	 SEM[0] = GetMeanError(_h_JetSpectrum1);
	 SEM[1] = GetMeanError(_h_JetSpectrum2);
	 SEM[2] = GetMeanError(_h_JetSpectrum3);
	 SEM[3] = GetMeanError(_h_JetSpectrum4);
	 SEM[4] = GetMeanError(_h_JetSpectrum5);
	 for(int i=0; i<5; i++){
	    _h_MeanJetPt->fill(MultBinCent[i], Mean[i], SEM[i]);
	 }
	//----------------------------------------------------------------
	 

	 for(int i=0; i<5; i++){
	   AvJetRate5[i]=JetCounter5GeV[i]/passedEv[i];
	   AvJetRate30[i]=JetCounter30GeV[i]/passedEv[i];
		 
	    if(JetCounter5GeV[i]!=0)SEM[i] =1/sqrt(JetCounter5GeV[i]);
	    else SEM[i]=0;
	    _h_JetRate5GeV->fill(MultBinCent[i], AvJetRate5[i], SEM[i]);
		
	    if(JetCounter30GeV[i]!=0)SEM[i] =1/sqrt(JetCounter30GeV[i]);
	    else SEM[i]=0;
	    _h_JetRate30GeV->fill(MultBinCent[i], AvJetRate30[i], SEM[i]);
	 }
	 
	 scale(_h_JetSpectrum1, 4.0/JetCounter5GeV[0]);
	 scale(_h_JetSpectrum2, 4.0/JetCounter5GeV[1]);
	 scale(_h_JetSpectrum3, 4.0/JetCounter5GeV[2]);
	 scale(_h_JetSpectrum4, 4.0/JetCounter5GeV[3]);
	 scale(_h_JetSpectrum5, 4.0/JetCounter5GeV[4]);
	 
	 scale(_h_JetStruct1, 0.08/JetStructNorm[0]);
	 scale(_h_JetStruct2, 0.08/JetStructNorm[1]);
	 scale(_h_JetStruct3, 0.08/JetStructNorm[2]);
	 scale(_h_JetStruct4, 0.08/JetStructNorm[3]);
	 scale(_h_JetStruct5, 0.08/JetStructNorm[4]);
	          
	 AIDA::IHistogramFactory& temp = histogramFactory();

	 temp.destroy(_AllTrkSpectrum1);
	 temp.destroy(_AllTrkSpectrum2);
	 temp.destroy(_AllTrkSpectrum3);
	 temp.destroy(_AllTrkSpectrum4);
	 temp.destroy(_AllTrkSpectrum5);

	 temp.destroy(_SoftTrkSpectrum1);
	 temp.destroy(_SoftTrkSpectrum2);
	 temp.destroy(_SoftTrkSpectrum3);
	 temp.destroy(_SoftTrkSpectrum4);
	 temp.destroy(_SoftTrkSpectrum5);

	 temp.destroy(_JetTrkSpectrum1);
         temp.destroy(_JetTrkSpectrum2);
	 temp.destroy(_JetTrkSpectrum3);
	 temp.destroy(_JetTrkSpectrum4);
	 temp.destroy(_JetTrkSpectrum5);

         temp.destroy(_JetLTrkSpectrum1);
         temp.destroy(_JetLTrkSpectrum2);
	 temp.destroy(_JetLTrkSpectrum3);
	 temp.destroy(_JetLTrkSpectrum4);
	 temp.destroy(_JetLTrkSpectrum5);
	 
    }


  private:
	
	//double JetStructNorm1,JetStructNorm2,JetStructNorm3,JetStructNorm4,JetStructNorm5;
	double JetStructNorm[5];
	double MultBinCent[5]; 
    	double AvJetRate5[5], AvJetRate30[5];//

	double Mean[5], SEM[5];

	// counters
	double JetCounter5GeV[5], JetCounter30GeV[5], passedEv[5];

	//histrograms--------------------------------------------

	AIDA::IProfile1D * _h_AllTrkMeanPt;
	AIDA::IProfile1D * _h_SoftTrkMeanPt;
	AIDA::IProfile1D * _h_IntrajetTrkMeanPt;
	AIDA::IProfile1D * _h_IntrajetLeaderTrkMeanPt;
	AIDA::IProfile1D * _h_MeanJetPt;
	AIDA::IProfile1D * _h_JetRate5GeV;
	AIDA::IProfile1D * _h_JetRate30GeV;
	
	AIDA::IHistogram1D * _h_JetSpectrum1;
	AIDA::IHistogram1D * _h_JetSpectrum2;
	AIDA::IHistogram1D * _h_JetSpectrum3;
	AIDA::IHistogram1D * _h_JetSpectrum4;
	AIDA::IHistogram1D * _h_JetSpectrum5;
	
	AIDA::IHistogram1D * _h_JetStruct1;
	AIDA::IHistogram1D * _h_JetStruct2;
	AIDA::IHistogram1D * _h_JetStruct3;
	AIDA::IHistogram1D * _h_JetStruct4;
	AIDA::IHistogram1D * _h_JetStruct5;

	//temp hisrograms---------------------------------------
	AIDA::IHistogram1D * _AllTrkSpectrum1;
	AIDA::IHistogram1D * _AllTrkSpectrum2;
	AIDA::IHistogram1D * _AllTrkSpectrum3;
	AIDA::IHistogram1D * _AllTrkSpectrum4;
	AIDA::IHistogram1D * _AllTrkSpectrum5;

	AIDA::IHistogram1D * _SoftTrkSpectrum1;
	AIDA::IHistogram1D * _SoftTrkSpectrum2;
	AIDA::IHistogram1D * _SoftTrkSpectrum3;
	AIDA::IHistogram1D * _SoftTrkSpectrum4;
	AIDA::IHistogram1D * _SoftTrkSpectrum5;

	AIDA::IHistogram1D * _JetTrkSpectrum1;
	AIDA::IHistogram1D * _JetTrkSpectrum2;
	AIDA::IHistogram1D * _JetTrkSpectrum3;
	AIDA::IHistogram1D * _JetTrkSpectrum4;
	AIDA::IHistogram1D * _JetTrkSpectrum5;

	AIDA::IHistogram1D * _JetLTrkSpectrum1;
	AIDA::IHistogram1D * _JetLTrkSpectrum2;
	AIDA::IHistogram1D * _JetLTrkSpectrum3;
	AIDA::IHistogram1D * _JetLTrkSpectrum4;
	AIDA::IHistogram1D * _JetLTrkSpectrum5;	    
	
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2013_I1261026);

}
