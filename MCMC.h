// MCMC
//
// This class is used to perform a Markov-Chain Monte Carlo.
//
// Inputs at instance creation:
//   Likelihood function (same format as for Minuit)
//   Parameter Names
//   Initial Parameter Guesses
//   Initial/Prior Parameter Uncertainty/or Fixed?
//   Length of Chain to run
//
// Outputs: (Ignore first few calls based on input)
//   TH1Ds:  each parameter distribution (for non fixed params)
//   TTree:  each proposed and accepted par value and likelihood at each step 
//   
// Blair Jamieson (Sept. 2007)
#ifndef _MCMC_h_
#define _MCMC_h_

//#include "TObject.h"
#include "TRint.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"

#include "MCMCFcn.h"

#include <cstdio>
#include <iostream>
using namespace std;

class MCMC {
  
 public:
  
  // constructor 
  // note that ParSigma is -1 for fixed parameteres
  MCMC( //void (*afcnLikelihood)(Int_t &, Double_t *,Double_t &, Double_t *,Int_t),
	  string aOutDir,
	  UInt_t aSeed,
	  Int_t      aNPar,
	  Char_t  ** aParNames,
	  Double_t * aParInit,
	  Double_t * aParSigma );
  
  // destructor cleanup memory use
  ~MCMC();

  // main work done when this is called to regenerate pdfs
  void RunMCMC(Int_t aChainLength, Int_t aNIgnore = 100, Int_t aPrintFreq=100);

  // print state
  void PrintState(Int_t i);

  // report output using routines similar to tminuit
  // mnpout returns parameter name, value and error
  void mnpout(Int_t ipar, TString &aparname, Double_t &aparval, Double_t &aparerr ){
    aparname = fParNames[ipar];
    if ( fParSigma[ipar] <= 0. ){
      aparval  = fParInit[ipar];
      aparerr  = 0.;
    } else {
      aparval = fParvsPar2[ipar]->GetMean();
      aparerr = fParvsPar2[ipar]->GetRMS();
    }
  }
  // mnerrs returns asymmetric + and - error estimate
  void mnerrs(Int_t ipar, Double_t &apluserr, Double_t &aminerr){
    if ( fParSigma[ipar] <= 0. ){
      apluserr  = 0.;
      aminerr   = 0.;
    } else {
      // use 1 sigma as 65% of area of PDF
      // or 32.5% to + of median and 332.5% to - of median.
      Double_t totarea = 0.325*fParvsPar2[ipar]->Integral("width");
      Int_t ibinc; Int_t i1,i2;
      fParvsPar2[ipar]->GetMaximumBin(ibinc,i1,i2);
      Int_t plusbin=ibinc;
      Int_t minubin=ibinc;
      Double_t plusarea = (fParvsPar2[ipar]->Integral(ibinc,ibinc,"width"))/2.;
      Double_t minusarea = plusarea;
      for (Int_t ibin=ibinc+1; ibin< fParvsPar2[ipar]->GetNbinsX(); ibin++ ){
	if (plusarea>totarea) break;
	plusbin++;
	plusarea += fParvsPar2[ipar]->Integral(ibin,ibin,"width");
      }
      apluserr = 
	fParvsPar2[ipar]->GetBinCenter(plusbin) - 
	fParvsPar2[ipar]->GetBinCenter(ibinc); 
      for (Int_t ibin=ibinc-1; ibin>=1; ibin-- ){
	if (minusarea>totarea) break;
	minubin--;
	minusarea += fParvsPar2[ipar]->Integral(ibin,ibin,"width");
      }
      aminerr = 
	fParvsPar2[ipar]->GetBinCenter(ibinc)-
	fParvsPar2[ipar]->GetBinCenter(minubin); 
    }
  }
 
 private:

  // Keep track of Inputs:
  string fOutDir;       // Output directory name
  UInt_t fSeed;         // Initial random number seed
  Int_t fNPar;          // Number of paramters
  Int_t fNFixed;        // Number of the parameters that are fixed
  Char_t ** fParNames;  //! Parameter names
  Double_t * fParInit;  //! Initial parameter values
  Double_t * fParSigma; //! Prior/proposal kernel (sigma of cauchy distr. about init)
  TF1 **     fPropKernel; //! Function for proposal kernel for this param.
  
  // Chain length values
  Int_t      fChainLength; // number of steps to take
  Int_t      fNIgnore;     // number of initial steps to ignore
  Int_t      fPrintFreq;   // printout frequency

  // Current values
  Int_t      fStep;
  Double_t * fParCurr;  //! Current Parameter values
  Double_t * fParProp;  //! Proposal Parameter values
  Double_t   fLogLCurr; // Current likelihood
  Double_t   fLogLProp; // Likelihood for proposal values
  Double_t   fProbAcc;  // Current probability of accepting step
  Double_t   fRandom;   // Current uniform random value

  // Output plots / storage of step info
  Int_t      fNPoints;   // number of points in plots
  Int_t      fNTH1D;     // number of combos

  TTree   *  fTree;      //
  TH1D    ** fParvsPar;    //! [Parnum1*fNPar+Parnum2] 2d histograms (accepted vals)
  TH1D    ** fParvsPar1;   //! [Parnum1*fNPar+Parnum2] 2d histograms (accepted vals)
  TH1D    ** fParvsPar2;   //! [Parnum1*fNPar+Parnum2] 2d histograms (accepted vals)

  // working variables 
  Double_t * fgin;       //!
  Int_t      fflag;
  TFile    * f2d;        //! TFile in which to put the 2d hists

};

#endif
