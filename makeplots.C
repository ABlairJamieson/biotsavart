// simple script to produce plots from root file output
void makeplots( const string & tag, const string & fname = "hyperkbs2.root" ){
  gROOT->LoadMacro("t2kstyle.h");
  t2kstyle();

  
  TFile * fin = new TFile( fname.c_str(), "read" );

  TH1D* hbmag = (TH1D*)fin->Get("init_Bmag");
  TCanvas * cbmag = new TCanvas;
  hbmag->SetLineWidth(3);
  hbmag->SetStats(1);
  hbmag->Draw();
  hbmag->GetXaxis()->SetRangeUser(0.0, 0.2);
  string fout = tag;
  fout += string("_Bmag.pdf");
  cbmag->Print(fout.c_str());

  TH1D* hbperp = (TH1D*)fin->Get("init_Bperp");
  TCanvas * cbperp = new TCanvas;
  hbperp->SetLineWidth(3);
  hbperp->SetStats(1);
  hbperp->Draw(); 
  hbperp->GetXaxis()->SetRangeUser(0.0, 0.2);
  fout = tag;
  fout += string("_Bperp.pdf");
  cbperp->Print(fout.c_str());

  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);

  TH2D* hphizB = (TH2D*)fin->Get("init_hphiz_side_B");
  TCanvas * cphizB = new TCanvas;
  hphizB->SetMaximum(0.1);
  hphizB->Draw("colz");
  fout = tag;
  fout += string("_phizB.pdf");
  cphizB->Print(fout.c_str());

  
  TH2D* hxytopB = (TH2D*)fin->Get("init_hxy_top_B");
  TCanvas * cxytopB = new TCanvas;
  hxytopB->SetMaximum(0.1);
  hxytopB->Draw("colz");
  fout = tag;
  fout += string("_xytopB.pdf");
  cxytopB->Print(fout.c_str());

}
