#include <TStyle.h>
#include <TROOT.h>
#include <TColor.h>
#ifndef _t2kstyle_h_
#define _t2kstyle_h_

// T2K style definition
// Adopted from BaBar collaboration
// Add the following lines to the start of your rootlogon.C file
//TStyle *t2kStyle= new TStyle("T2K","T2K approved plots style");
void t2kstyle(){  
// use plain black on white colors
gStyle->SetFrameBorderMode(0);
gStyle->SetCanvasBorderMode(0);
gStyle->SetPadBorderMode(0);
 
gStyle->SetPadColor(0);
gStyle->SetCanvasColor(0);
gStyle->SetStatColor(0);
 
if (0) gStyle->SetFillColor(kWhite);
gStyle->SetLegendBorderSize(1); 
 
// set the paper & margin sizes
gStyle->SetPaperSize(20,26);
gStyle->SetPadTopMargin(0.05);
gStyle->SetPadRightMargin(0.05);
gStyle->SetPadBottomMargin(0.16);
gStyle->SetPadLeftMargin(0.2);

// use large Times-Roman fonts
gStyle->SetTextFont(132);
gStyle->SetTextSize(0.08);
gStyle->SetLabelFont(132,"x");
gStyle->SetLabelFont(132,"y");
gStyle->SetLabelFont(132,"z");
gStyle->SetLabelSize(0.05,"x");
gStyle->SetTitleSize(0.06,"x");
gStyle->SetLabelSize(0.05,"y");
gStyle->SetTitleSize(0.06,"y");
gStyle->SetLabelSize(0.05,"z");
gStyle->SetTitleSize(0.06,"z");
gStyle->SetLabelFont(132,"t");
gStyle->SetTitleFont(132,"x");
gStyle->SetTitleFont(132,"y");
gStyle->SetTitleFont(132,"z");
gStyle->SetTitleFont(132,"t"); 
gStyle->SetTitleFillColor(0);
gStyle->SetTitleX(0.25);
gStyle->SetTitleFontSize(0.08);
gStyle->SetTitleFont(132,"pad");
 
// use bold lines and markers
gStyle->SetMarkerStyle(20);
gStyle->SetHistLineWidth(2);
gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

// get rid of X error bars and y error bar caps
gStyle->SetErrorX(0.001);



// do not display any of the standard histogram decorations
gStyle->SetOptTitle(0);
gStyle->SetOptStat(0);
gStyle->SetOptFit(0);

// put tick marks on top and RHS of plots
gStyle->SetPadTickX(1);
gStyle->SetPadTickY(1);
  

// Add a greyscale palette for 2D plots
 if (0){
   int ncol=50;
   double dcol = 1./float(ncol);
   double gray = 1;
   TColor **theCols = new TColor*[ncol];
   for (int i=0; i<ncol; i++) {
     theCols[i] = new TColor(999-i, 0.0, 0.7, 0.7);
   }
   for (int j = 0; j < ncol; j++) {
     theCols[j]->SetRGB(gray,gray,gray);
     gray -= dcol;
   }
   int ColJul[100];
   for  (int i=0; i<100; i++) ColJul[i]=999-i;
   gStyle->SetPalette(ncol,ColJul);
   
 } else {
// Define a nicer color palette (red->blue)
// Uncomment these lines for a color palette (default is B&W)
//gStyle->SetPalette(1,0);  // use the nice red->blue palette
 const Int_t NRGBs = 5;
 const Int_t NCont = 255;

 Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
 Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
 Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
 Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
 TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
 gStyle->SetNumberContours(NCont); 

 }
 //gROOT->SetStyle("T2K");
}
// End of definition of t2kStyle
#endif


