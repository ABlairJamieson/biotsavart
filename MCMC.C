
#include "MCMC.h"

  
MCMC::MCMC( 
		 string     aOutDir,
		 UInt_t     aSeed,
		 Int_t      aNPar,
		 Char_t  ** aParNames,
		 Double_t * aParInit,
		 Double_t * aParSigma ){


  std::cout<<"<MCMC::MCMC> NPar="<<aNPar
	   <<" oudir="<<aOutDir.c_str()<<" seed="<<aSeed<<std::endl;
  
  // make copies of inputs
  fOutDir = aOutDir;
  fSeed = aSeed;
  fNPar = aNPar;
  fParNames = new Char_t*[aNPar];
  fParInit = new Double_t[aNPar];
  fParSigma = new Double_t[aNPar];
  fNFixed = 0;
  for (Int_t i=0;i<aNPar;i++){
    std::cout<<"<MCMC::MCMC> iPar="<<i
	     <<" name="<<aParNames[i]
	     <<" Init="<<aParInit[i]
	     <<" Sigma="<<aParSigma[i]
	     <<std::endl;

    fParNames[i] = new Char_t[ strlen(aParNames[i])+1 ];
    strcpy( fParNames[i], aParNames[i] );
    fParInit[i] = aParInit[i];
    fParSigma[i] = aParSigma[i];
    if ( fParSigma[i] <= 0. ) fNFixed++; 
  }

  // make proposal kernels
  Char_t aname[256];
  fPropKernel = new TF1*[aNPar];
  for (Int_t i=0; i<aNPar; i++){
    fPropKernel[i] = NULL;
    if ( fParSigma[i] > 0. ) {
      sprintf(aname,"fstudentt_%03d",i);
      fPropKernel[i] = new TF1(aname,"gaus",-fParSigma[i]*10.,fParSigma[i]*10.);
      fPropKernel[i]->SetParameters(1.,0.,fParSigma[i]/3.);

      std::cout<<"<MCMC::MCMC> Setup Proposal Kernel "<<aname<<" sigma="<<fParSigma[i]/3.<<std::endl;
      std::cout<<"                   a few randoms from this kernel: ";
      for (Int_t ii=0; ii<5; ii++) std::cout<<fPropKernel[i]->GetRandom()<<" ";
      std::cout<<std::endl;
    }
  }

  // Save 2d hists here
  string anOutFile = fOutDir + "/MCMC.root";
  f2d=new TFile(anOutFile.c_str(),"recreate");

  // make TH1Ds
  fNTH1D = fNPar;
  fParvsPar = new TH1D*[fNTH1D];
  fParvsPar1 = new TH1D*[fNTH1D];
  fParvsPar2 = new TH1D*[fNTH1D];
  std::cout<<"<MCMC::MCMC> Number of TH2Ds of Par vs Par is "<<fNTH1D<<std::endl;  

  // setup parameter storage
  fParCurr = new Double_t[fNPar];
  fParProp = new Double_t[fNPar];
  for (Int_t i=0; i<fNPar; i++){
    fParCurr[i] = fParInit[i];
    fParProp[i] = fParInit[i];
  }

  // make the tree
  TTree * ffTree = new TTree("Tmcmc","MCMC results");
  fTree = ffTree;
  // general branches in all versions
  fTree->Branch("Step",&fStep,"Step/I");
  fTree->Branch("LogLPro",&fLogLProp,"LogLPro/D");
  fTree->Branch("LogLAcc",&fLogLCurr,"LogLAcc/D");
  fTree->Branch("ProbAcc",&fProbAcc,"ProbAcc/D");
  fTree->Branch("Rndm",&fRandom,"Rndm/D");
  // branches for proposed parameter values
  Char_t bname1[256];
  Char_t bname2[256];
  for (Int_t i=0; i<fNPar; i++){
    if ( fParSigma[i]>0. ){
      sprintf(bname1,"Pro%s",fParNames[i]);
      sprintf(bname2,"Pro%s/D",fParNames[i]);
      fTree->Branch(bname1,&(fParProp[i]),bname2);
      sprintf(bname1,"Acc%s",fParNames[i]);
      sprintf(bname2,"Acc%s/D",fParNames[i]);
      fTree->Branch(bname1,&(fParCurr[i]),bname2);
    }
  }
  fTree->Print();

  // Initialize current parameters
  fgin = new Double_t[fNPar];
  fflag = 0;
  //fcnLikelihood( fNPar, fgin, fLogLCurr, fParCurr, fflag );
  std::cout<<"Calling Likelihood"<<std::endl;
  MCMCFcn( fNPar, fgin, fLogLCurr, fParCurr, fflag );
  fLogLProp = fLogLCurr;
  
  // Initialize rest to zero or NULL
  fNPoints = 0;
  fChainLength = 0;
  fNIgnore = 0;
  
  return;
}

MCMC::~MCMC(){
  f2d->Close();
  if ( fParNames != NULL ){
    for (Int_t i=0; i<fNPar; i++){
      delete [] fParNames[i];
    }
    delete [] fParNames;
  }
  if (fParInit!=NULL) delete [] fParInit;
  if (fPropKernel!=NULL) {
    for (Int_t i=0; i<fNPar; i++){
      if (fParSigma[i]>0.&&fPropKernel[i]!=NULL) delete fPropKernel[i];
    }
    delete [] fPropKernel;
  }

  // leave TH2D and graphs alone, they went to files
  if (fParSigma!=NULL) delete [] fParSigma;
  if (fParCurr!=NULL) delete [] fParCurr;
  if (fParProp!=NULL) delete [] fParProp;

  fNPar = fNFixed = fChainLength = fNIgnore = fNPoints = fNTH1D = 0;

  return;
}

void MCMC::PrintState(Int_t i){
  std::cout<<"<MCMC::PrintState> Step "<<i<<" / "<<fChainLength<<std::endl;
  std::cout<<"  Parameter    Current   Proposed"<<std::endl;
  for (Int_t ipar=0; ipar<fNPar; ipar++){
    std::cout<<"  "<<fParNames[ipar]<<" "<<fParCurr[ipar]<<" "<<fParProp[ipar]<<std::endl;
  }
  std::cout<<"  logL         Current   Proposed  Probability of accept"<<std::endl;
  std::cout<<"              "<<fLogLCurr<<" "<<fLogLProp<<" "<<fProbAcc<<std::endl;
  std::cout<<"-------------------------------------------"<<std::endl;
  return;
}

void MCMC::RunMCMC(Int_t aChainLength, Int_t aNIgnore, Int_t aPrintFreq ){
  // Make Copy of inputs
  fChainLength = aChainLength;
  fNIgnore = aNIgnore;
  fPrintFreq = aPrintFreq;
  fNPoints = fChainLength - fNIgnore;
  
  std::cout<<"<MCMC::RunMCMC> ChainLength="<<fChainLength<<" Ingnore First "<<fNIgnore<<" steps"<<std::endl;
  
  // Setup storage
  Double_t * accmin  = new Double_t[fNPar];
  Double_t * accmax  = new Double_t[fNPar];

  for (Int_t i=0; i<fNPar; i++){
    accmin[i] =  999999.0;
    accmax[i] = -999999.0;
  }

  // do MCMC
  Double_t aMinLogl=999999.;
  gRandom->SetSeed(fSeed);
  std::cout<<"<MCMC::MCMC> gRandom using Seed "<<gRandom->GetSeed()<<std::endl;
  for (Int_t i=0; i<fChainLength; i++){
    // generate new proposal parameters
    for (Int_t ipar=0; ipar<fNPar; ipar++){
      if ( fParSigma[ipar] > 0. ){
	fParProp[ipar] = fParCurr[ipar] + 
	  fPropKernel[ipar]->GetRandom();
      }	
    }

    // calculate likelihood for proposal parameters
    MCMCFcn( fNPar, fgin, fLogLProp, fParProp, fflag );

    // calculate acceptance probability
    fProbAcc = TMath::Min(1.,TMath::Exp(fLogLCurr-fLogLProp));
    //
    if (i%fPrintFreq==0) PrintState(i);

    fRandom = gRandom->Rndm();
    if ( fRandom <= fProbAcc ){
      // keep this step
      fLogLCurr = fLogLProp;
      if ( fLogLCurr<aMinLogl ) aMinLogl = fLogLCurr;
      for (Int_t ipar=0; ipar<fNPar; ipar++){
	fParCurr[ipar] = fParProp[ipar];
      }
    }

    //
    // Save the current step info in the tree
    fStep = i;
    fTree->Fill();

    for (Int_t ipar=0; ipar<fNPar; ipar++){
      if ( accmin[ipar] > fParCurr[ipar] ) accmin[ipar] = fParCurr[ipar];
      if ( accmax[ipar] < fParCurr[ipar] ) accmax[ipar] = fParCurr[ipar];     
    }

    // autosave every 500 steps
    if ( fStep%500 == 499 ) fTree->AutoSave() ;

    // done this step, repeat
  }
  

  std::cout<<std::endl
	   <<std::endl
	   <<std::endl
	   <<std::endl
	   <<"<MCMC::RunMCMC> Done MCMC, make graphs"
	   <<std::endl
	   <<std::endl
	   <<std::endl
	   <<std::endl
	   <<std::endl
	   <<std::endl
	   <<std::endl
	   <<std::endl
	   <<std::endl
	   <<std::endl
	   <<std::endl;

  // Save the 2d hists 
  f2d->cd();
  // (that way if things explode making likelihood graphs have some results)
  Char_t aname[256];
  Char_t atitle[256];
  std::cout<<"<MCMC::MCMC> Number of TH2Ds of Par vs Par is "<<fNTH1D<<std::endl;
  Int_t aidx=0;
  for (Int_t i=0; i<fNPar; i++){
    fParvsPar[aidx] = NULL;
    fParvsPar1[aidx] = NULL;
    fParvsPar2[aidx] = NULL;
    if ( fParSigma[i] > 0. ){
      if ( aidx >= fNTH1D ){
	std::cout<<"ERROR! aidx >= fNTH1D Exiting!"<<std::endl;
	exit(0);
      }
      sprintf(aname,"h%s%s",fParNames[i],fParNames[i]);
      sprintf(atitle,"%s All",fParNames[i]);
      fParvsPar[aidx] = 
	new TH1D(aname,atitle,
		 100,accmin[i],accmax[i]);
      fParvsPar[aidx]->Sumw2();
      sprintf(aname,"h%s%s1",fParNames[i],fParNames[i]);
      sprintf(atitle,"%s First Half",fParNames[i]);
      fParvsPar1[aidx] = 
	new TH1D(aname,atitle,
		 100,accmin[i],accmax[i]);
      fParvsPar1[aidx]->SetLineColor(2);
      fParvsPar1[aidx]->SetMarkerColor(2);
      fParvsPar1[aidx]->SetMarkerStyle(7);
      fParvsPar1[aidx]->Sumw2();
      
      sprintf(aname,"h%s%s2",fParNames[i],fParNames[i]);
      sprintf(atitle,"%s Second Half",fParNames[i]);
      fParvsPar2[aidx] = 
	new TH1D(aname,atitle,
		 100,accmin[i],accmax[i]);
      fParvsPar2[aidx]->SetLineColor(4);
      fParvsPar2[aidx]->SetMarkerColor(4);
      fParvsPar2[aidx]->SetMarkerStyle(7);
      fParvsPar2[aidx]->Sumw2();
      
    }      
    aidx++;
  }

  // calculate correlation matrix while filling histoz
  // not very efficient but coded quickly!
  Double_t *dsumxi  = new Double_t[ (fNPar-fNFixed)*(fNPar-fNFixed) ];
  Double_t *dsumyi  = new Double_t[ (fNPar-fNFixed)*(fNPar-fNFixed) ];
  Double_t *dmeanxi  = new Double_t[ (fNPar-fNFixed)*(fNPar-fNFixed) ];
  Double_t *dmeanyi  = new Double_t[ (fNPar-fNFixed)*(fNPar-fNFixed) ];
  Double_t *dxiyi   = new Double_t[ (fNPar-fNFixed)*(fNPar-fNFixed) ];
  Double_t *dcovmat = new Double_t[ (fNPar-fNFixed)*(fNPar-fNFixed) ];
  Double_t *dcormat = new Double_t[ (fNPar-fNFixed)*(fNPar-fNFixed) ];
  for (Int_t ii=0; ii<(fNPar-fNFixed)*(fNPar-fNFixed); ii++){
    dsumxi[ii] = 0.;
    dsumyi[ii] = 0.;
    dmeanxi[ii] = 0.;
    dmeanyi[ii] = 0.;
    dxiyi[ii] = 0.;
    dcovmat[ii] = 0.;
    dcormat[ii] = 0.;
  }

  for (Long64_t iev=fNIgnore; iev<fTree->GetEntries(); iev++){
    fTree->GetEntry(iev);
    aidx=0;


    Int_t iipar=-1;
    for (Int_t ipar=0; ipar<fNPar; ipar++){
      if (fParSigma[ipar]>0.) iipar++;
      Int_t jjpar=-1;
      for (Int_t jpar=0; jpar<fNPar; jpar++){      
	if ( fParSigma[ipar]>0. && fParSigma[jpar]>0. ){
	  jjpar++;
	  dsumxi[iipar*(fNPar-fNFixed)+jjpar] += fParCurr[ipar];
	  dsumyi[iipar*(fNPar-fNFixed)+jjpar] += fParCurr[jpar];
	  dxiyi[iipar*(fNPar-fNFixed)+jjpar] += fParCurr[ipar] * fParCurr[jpar];
	}
      }
    }
  }

  aidx=0;
  for (Int_t ipar=0; ipar<fNPar-fNFixed; ipar++){
    for (Int_t jpar=0; jpar<fNPar-fNFixed; jpar++){      
      aidx = ipar*(fNPar-fNFixed)+jpar;
      dmeanxi[aidx] = dsumxi[aidx]/float(fNPoints);
      dmeanyi[aidx] = dsumyi[aidx]/float(fNPoints);
      dcovmat[aidx] = 
	( dxiyi[aidx] - 
	  dmeanxi[aidx]*dsumyi[aidx] - 
	  dmeanyi[aidx]*dsumxi[aidx] )/float(fNPoints)
	+ dmeanxi[aidx]*dmeanyi[aidx];
      aidx++;
    }
  }
  for (Int_t ipar=0; ipar<fNPar-fNFixed; ipar++){
    for (Int_t jpar=0; jpar<fNPar-fNFixed; jpar++){
      dcormat[ipar*(fNPar-fNFixed)+jpar] =
	dcovmat[ipar*(fNPar-fNFixed)+jpar] / 
	sqrt( dcovmat[ipar*(fNPar-fNFixed)+ipar] *
	      dcovmat[jpar*(fNPar-fNFixed)+jpar] );
    }
  }
  
  // Fill histograms
  for (Long64_t iev=fNIgnore; iev<fTree->GetEntries(); iev++){
    fTree->GetEntry(iev);
    aidx=0;
    for (Int_t ipar=0; ipar<fNPar; ipar++){
      if ( fParSigma[ipar]>0. ){
	fParvsPar[aidx]->Fill( fParCurr[ipar] );
	if ( Int_t(iev)-fNIgnore >= fChainLength/2 ){
	  fParvsPar2[aidx]->Fill( fParCurr[ipar] );
	} else {
	  fParvsPar1[aidx]->Fill( fParCurr[ipar] );
	}
      }
      aidx++;
    }
  }
  
  //setup output file
  string anOutFile = fOutDir + "/FitResult.txt";
  FILE *ou = fopen(anOutFile.c_str(),"w");
  if(!ou){
    cout << "<MCMC::WriteFitResult> ERROR: file could not be opened." << endl;
    return;
  } 

  fprintf(ou,"# MCMC Status : Successful\n" );
  fprintf(ou,"# MCMC MinimumLogLikelihood : %f\n", aMinLogl );


  // Make 1d histograms, split into two histograms
  // one for first  half of event chain, one for second half
  std::cout<<"# -----------------------------------------------------------"<<std::endl;
  std::cout<<"# MCMC Fit Results:"<<std::endl;
  std::cout<<"# Parameter   Second-Half +-  First-Half +-  (First-Second)/sqrt(sigmas)"<<std::endl;
  fprintf(ou,"# -----------------------------------------------------------\n");
  fprintf(ou,"# MCMC Fit Results:\n");
  fprintf(ou,"# Parameter   All  +-  Second-Half +-  First-Half +-  (First-Second)/sqrt(sigmas)\n");
  aidx=0;
  for (Int_t ipar=0; ipar<fNPar; ipar++){
    if ( fParSigma[ipar]<=0. ){
      fprintf(ou,"% 3d % 20s % 12.6g % 12.6g  % 12.6g % 12.6g fixed\n",ipar,fParNames[ipar],fParInit[ipar],0.,fParInit[ipar],0.);
    } else {
      // Now write fit results...
      std::cout<<fParNames[ipar]<<"  "
	       <<fParvsPar1[ipar]->GetMean()<<"+-"<<fParvsPar1[ipar]->GetRMS()<<"  "
	       <<fParvsPar2[ipar]->GetMean()<<"+-"<<fParvsPar2[ipar]->GetRMS()<<"  "

	       <<std::endl;
      
      fprintf(ou,"% 3d % 20s % 12.6g % 12.6g  % 12.6g % 12.6g  % 12.6g % 12.6g % 12.6g\n",ipar,fParNames[ipar],
	      fParvsPar[ipar]->GetMean(),fParvsPar[ipar]->GetRMS(),
	      fParvsPar2[ipar]->GetMean(),fParvsPar2[ipar]->GetRMS(),
	      fParvsPar1[ipar]->GetMean(),fParvsPar1[ipar]->GetRMS(),
	      (fParvsPar[ipar]->GetMean()-fParvsPar2[ipar]->GetMean())/sqrt(fParvsPar[ipar]->GetRMS()*fParvsPar[ipar]->GetRMS()+fParvsPar2[ipar]->GetRMS()*fParvsPar2[ipar]->GetRMS()));
      
    }
    aidx++;
  }

  std::cout<<"-----------------------------------------------------------"<<std::endl;


  // write out covariance matrix
  fprintf(ou,"\nMCMC Fit covariance matrix:\n");
  fprintf(ou,"            ");
  for (Int_t ipar=0; ipar<fNPar; ipar++){
    if ( fParSigma[ipar]>0. ) fprintf(ou,"% 12s  ",fParNames[ipar]);
  }
  aidx=0;
  Int_t iipar=-1;
  Int_t jjpar=-1;
  for (Int_t ipar=0; ipar<fNPar; ipar++){
    if ( fParSigma[ipar]>0. ) {
      fprintf(ou,"\n% 20s  ",fParNames[ipar]);
      iipar++;
    }
    jjpar=-1;
    for (Int_t jpar=0; jpar<fNPar; jpar++){
      if ( fParSigma[ipar]>0. && fParSigma[jpar]>0.){
	jjpar++;
	fprintf(ou,"% 12.6g  ",dcovmat[iipar*(fNPar-fNFixed)+jjpar]);
      }
      aidx++;
    }
  }

  // Write correlation matrix
  fprintf(ou,"\n\nMCMC Fit correlation matrix:\n");
  fprintf(ou,"            ");
  for (Int_t ipar=0; ipar<fNPar; ipar++){
    if ( fParSigma[ipar]>0. ) fprintf(ou,"% 12s  ",fParNames[ipar]);
  }
  aidx=0;
  iipar=-1;
  jjpar=-1;
  for (Int_t ipar=0; ipar<fNPar; ipar++){
    if ( fParSigma[ipar]>0. ) {
      fprintf(ou,"\n% 20s  ",fParNames[ipar]);
      iipar++;
    }
    jjpar=-1;
    for (Int_t jpar=0; jpar<fNPar; jpar++){
      if ( fParSigma[ipar]>0. && fParSigma[jpar]>0.){
	jjpar++;
	fprintf(ou,"% 12.6g  ",dcormat[iipar*(fNPar-fNFixed)+jjpar]);
      }
      aidx++;
    }
  }


  f2d->Write();
  
  std::cout<<"  .... cleanup .... "<<std::endl;

  delete [] accmin;
  delete [] accmax;
}
