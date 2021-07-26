/*
 * hyperkbs2.C
 * Program to calculate magnetic field of Winnipeg PTF
 *
 * Author: Blair Jamieson (Mar. 2018, Jul. 2021)
 *
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include "TStyle.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLine.h"
#include "t2kstyle.h"
#include "TGraph2D.h"
#include "TArrow.h"

using std::cout;
using std::endl;

const double pi = std::acos(-1.0);
const double mu0_over_4pi = 1.0e-7;
const double T_to_mgauss = 1.0e4 * 1.0e3; 

// Earth's field in this problem is (303, 0, -366  ) mgauss
// where +z is vertical direction, +x is direction of B-field in horizontal plane
const TVector3 BEarth{ 0.303e-4, 0.0, -0.366e-4 }; // Tesla
//const TVector3 BEarth{ 0.0, 0.0, -0.366e-4 };
//const TVector3 BEarth{ 0.0, 0.0, 0.0 };

// HyperK dimensions
//const double hk_cavern_H = 60.0; //m
//const double hk_cavern_R = 74.0/2.0; //m

// Use root TVector3 to represent vectors

// represent wire element by vector of length three
// all lenghts in meters!
struct WireElement{
  WireElement();
  // constructors take vector dl and location rp
  WireElement( TVector3 adl, TVector3 arp ) : dl( adl ), rp( arp ){;}

  // get back info about element
  TVector3 get_dl() const { return dl; }
  TVector3 get_rp() const { return rp; }
 
  // method to get dl cross ( r-r' )
  void get_dl_x_r( const TVector3 & ar, TVector3 & dlxr ) const;

  // method to get dl cross ( r-r' ) / |r-r'|^3
  void get_dl_x_r3( const TVector3 & ar, TVector3 & dlxr ) const;
 
private:
  TVector3 dl; // Direction
  TVector3 rp; // Location of current element  
};

typedef std::vector< WireElement > Wire;

// method to get dl cross ( r-r' )
void WireElement::get_dl_x_r( const TVector3 & ar, TVector3 & dlxr ) const{
  TVector3 rr =ar - rp;
  dlxr = dl.Cross( rr );
}

// method to get dl cross ( r-r' ) / |r-r'|^3
void WireElement::get_dl_x_r3( const TVector3 & ar, TVector3 & dlxr ) const{
    TVector3 rr =ar - rp;
    dlxr = dl.Cross( rr );
    dlxr *= 1.0/( pow( rr.Mag(), 3 ) );
}


// Rectangular coil
struct RectangularCoil{
  
RectangularCoil() :
  fN(100), fL(1.0), fW(1.0), fDirL(TVector3(1.,0.,0.)), fDirW(TVector3(0.,1.,0.)),
    fCenter(TVector3(0.,0.,1.)) { InitWire(); }
  
RectangularCoil( int aN, float aL, float aW, const TVector3& aDirL, const TVector3& aDirW, const TVector3& aCenter ) :
  fN(aN), fL(aL), fW(aW), fDirL(aDirL), fDirW(aDirW), fCenter(aCenter) { InitWire(); }

  void set_coil( int aN, float aL, float aW, const TVector3& aDirL, const TVector3& aDirW, const TVector3& aCenter ) {
       fN = aN;       fL = aL;          fW = aW;
    fDirW = aDirW; fDirL = aDirL;  fCenter = aCenter;
    InitWire();
  }
  
  Wire get_wire() const { return fWire; }
  
private:
  void InitWire();
  int   fN;         // number of wire segments per side.
  float fL;         // coil length (meters)
  float fW;         // coil width (meters)
  TVector3 fDirL;   // direction unit vector along length of coil
  TVector3 fDirW;   // direction unit vector along width of coil
  TVector3 fCenter; // position of coil center (meters)
  Wire  fWire;
};


void RectangularCoil::InitWire(){
  fWire.clear();

  // step along L then along W then -L finally -W
  float dL = fL/fN;
  float dW = fW/fN;
  // going in +L
  for ( int i=0 ; i<fN; ++i ){
    TVector3 loc = fCenter + // center of coil
      -fW/2 * fDirW + // below center by 1/2 W direction
      -fL/2 * fDirL + dL * (i+0.5) * fDirL;  // stepping along L
    TVector3 dir = fDirL;
    fWire.push_back( WireElement( dir, loc ) );
  }
  // going in +W
  for ( int i=0 ; i<fN; ++i ){
    TVector3 loc = fCenter + // center of coil
      fL/2 * fDirL + // right of center by 1/2 W direction
      -fW/2 * fDirW + dW * (i+0.5) * fDirW;  // stepping along W
    TVector3 dir = fDirW;
    fWire.push_back( WireElement( dir, loc ) );
  }
  // going in -L
  for ( int i=0 ; i<fN; ++i ){
    TVector3 loc = fCenter + // center of coil
      fW/2 * fDirW + // above center by 1/2 W direction
      fL/2 * fDirL - dL * (i+0.5) * fDirL;  // stepping along -L
    TVector3 dir = fDirL;
    fWire.push_back( WireElement( dir, loc ) );
  }
  // going in -W
  for ( int i=0 ; i<fN; ++i ){
    TVector3 loc = fCenter + // center of coil
      -fL/2 * fDirL + // left of center by 1/2 W direction
      fW/2 * fDirW - dW * (i+0.5) * fDirW;  // stepping along -W
    TVector3 dir = fDirW;
    fWire.push_back( WireElement( dir, loc ) );
  }
}

  


// class to calculate magnetic field at a point.
// initialize class with:
//   1) current in amps
//   2) std::vector< WireElement > to represent locations and directions of current
// methods then exist to get field calculated at:
//   1) a point (TVector3)
//   2) a vector of points ( std::vector< TVector3 > )
//   3) as a TH1D along a line
//   4) as a TH2D on one of planes
struct BiotSavart{
  BiotSavart( double aI, Wire ww ) :
    I(aI), w( ww ) { ; }

  // field at point
  TVector3 get_B( const TVector3 &p ) const;

  // field at vector of points
  std::vector< TVector3 > get_B( const std::vector< TVector3 > & p ) const;

  // Build and fill a 1D histogram of |B| along line
  // ceneter at point loc and going in direction dir (dir must be normalized)
  TH1D* get_B( const TVector3 &  loc, const TVector3 & dir,
	       const int nbinsx, const double xmin, const double xmax ) const;

  // Build and fill 2D histogram of B in plane
  // Define plane by a center location and two (presumed to be orthogonal) vectors as TVector3
  // assumes xdir and ydir are normalized vectors
  // Define region of histogram as coordinates relative to center location
  TH2D* get_B( const TVector3 & c, const TVector3 & xdir, const TVector3 & ydir,
	       const int nbinsx, const double xmin, const double xmax,
	       const int nbinsy, const double ymin, const double ymax ) const;

  const Wire & get_Wire() const { return w; }
  double get_I() const { return I; }
  
private:
  double I;
  Wire w;
};

// field at point
TVector3 BiotSavart::get_B( const TVector3& p ) const{
  TVector3 B(0.,0.,0.);
  TVector3 dB(0.,0.,0.);
  for ( auto elem : w ){
    elem.get_dl_x_r3( p, dB );
    B += dB;
  }
  B *= mu0_over_4pi * I;
  return B;
}

// field at vector of points
std::vector< TVector3 > BiotSavart::get_B( const std::vector< TVector3 > & vp ) const{
  std::vector< TVector3 > B;
  for ( const auto p : vp ){
    B.push_back( get_B( p ) );
  }
  return B;
}

// Build and fill a 1D histogram of B along line
// center at point loc and going in direction dir (dir must be normalized)
TH1D* BiotSavart::get_B( const TVector3 &  loc, const TVector3 & dir,
			 const int nbinsx, const double xmin, const double xmax ) const {
  static int count=0;
  count++;
  std::string hname{"get_B1D"};
  hname += std::to_string( count ); 
  TH1D* h = new TH1D( hname.c_str(), " ; loc (m) ; B ( Tesla )", nbinsx, xmin, xmax );
  for ( int i=1; i<=nbinsx; ++i){
    double bc = h->GetBinCenter( i );
    TVector3 p = loc + bc*dir;
    TVector3 B = get_B( p );
    h->SetBinContent( i, B.Mag() );
  }
  return h;
}

// Build and fill 2D histogram of B in plane
// Define plane by a center location and two (presumed to be orthogonal) vectors as TVector3
// assumes xdir and ydir are normalized vectors
// Define region of histogram as coordinates relative to center location
TH2D* BiotSavart::get_B( const TVector3 & c, const TVector3 & xdir, const TVector3 & ydir,
			 const int nbinsx, const double xmin, const double xmax,
			 const int nbinsy, const double ymin, const double ymax ) const {
  static int count=0;
  count++;
  std::string hname{"get_B2D"};
  hname += std::to_string( count ); 
  TH2D* h = new TH2D( hname.c_str(), " ; loc (m) ; loc (m); B ( Tesla )",
		      nbinsx, xmin, xmax, nbinsy, ymin, ymax );
  TAxis * xax = h->GetXaxis();
  TAxis * yax = h->GetYaxis();
  for ( int i=1; i<=nbinsx; ++i){
    for ( int j=1; j<=nbinsy; ++j){
      double x = xax->GetBinCenter( i );
      double y = yax->GetBinCenter( j );
      TVector3 p = c + x * xdir + y * ydir;
      TVector3 B = get_B( p );
      h->SetBinContent( i, j, B.Mag() );
    }
  }
  return h;
}

// coil types is vector of strings "XY", "XZ", or "YZ"
void plot_coils( const std::vector< RectangularCoil >& Coils, const std::vector< std::string > & coil_types ) {
  std::vector<double> x,y,z;

  for ( unsigned i=0 ; i<coil_types.size(); ++i ){
    std::ostringstream osname;
    osname << "hcoil" << i;
    std::ostringstream ostitle;
    ostitle << coil_types[i] <<" Coil; "<< coil_types[i][0] <<" (m); "<< coil_types[i][1] <<" m";
    TH2D* harrow = new TH2D( osname.str().c_str(), ostitle.str().c_str(), 100, -2., 2., 100, -2., 2. );

    Wire wr = Coils[i].get_wire();
    for (unsigned j=1; j<wr.size(); ++j ){
      float x0=0., y0=0., x1=0., y1=0.;
      switch ( coil_types[i][0] ) {
      case 'X':
	x0 = wr[j-1].get_rp().X();
	x1 = wr[j].get_rp().X();
	break;
      case 'Y':
	x0 = wr[j-1].get_rp().Y();
	x1 = wr[j].get_rp().Y();
	break;
      default:
	x0 = wr[j-1].get_rp().Z();
	x1 = wr[j].get_rp().Z();
	break;
      }
      switch ( coil_types[i][1] ) {
      case 'X':
	y0 = wr[j-1].get_rp().X();
	y1 = wr[j].get_rp().X();
	break;
      case 'Y':
	y0 = wr[j-1].get_rp().Y();
	y1 = wr[j].get_rp().Y();
	break;
      default:
	y0 = wr[j-1].get_rp().Z();
	y1 = wr[j].get_rp().Z();
	break;
      }
      TArrow * ta = new TArrow(x0,y0,x1,y1);
      harrow->GetListOfFunctions()->Add( (TObject*)ta );

    }
    
  }
  
  for ( const RectangularCoil& rec_coil : Coils ){
    Wire wr = rec_coil.get_wire();
    

    
    for ( const WireElement& we : wr ){
      TVector3 loc = we.get_rp();
      x.push_back( loc.X() );
      y.push_back( loc.Y() );
      z.push_back( loc.Z() );
    }
  }
  TGraph2D * tg2d = new TGraph2D( x.size(), &x[0], &y[0], &z[0] );
  tg2d->SetName("tg_coilloc");
  tg2d->SetTitle("Coil Locations; X(m); Y(m); Z(m)");

  tg2d->Write();
}

void biotsavart(){
  t2kstyle();
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  
  // Output file
  TFile * fout = new TFile("biotsavart.root","recreate");

  // Make a list of coils
  std::vector< RectangularCoil > Coils;

  // coil dimensions
  const float xcoil_len = 1.310; // m
  const float ycoil_len = 1.260; // m
  const float zcoil_len = 1.360; // m

  // define coils
  //                                     L          W          dirL                dirW                center
  // xcoils (in YZ plane shifted in x)
  Coils.push_back( RectangularCoil( 100, xcoil_len, xcoil_len, TVector3(0.,1.,0.), TVector3(0.,0.,1.), TVector3(-ycoil_len/2,0.,0.) ) );
  Coils.push_back( RectangularCoil( 100, xcoil_len, xcoil_len, TVector3(0.,1.,0.), TVector3(0.,0.,1.), TVector3( ycoil_len/2,0.,0.) ) );
  // ycoils (in XZ plane shifted in y)
  Coils.push_back( RectangularCoil( 100, ycoil_len, ycoil_len, TVector3(1.,0.,0.), TVector3(0.,0.,1.), TVector3(0.,-xcoil_len/2,0.) ) );
  Coils.push_back( RectangularCoil( 100, ycoil_len, ycoil_len, TVector3(1.,0.,0.), TVector3(0.,0.,1.), TVector3(0., xcoil_len/2,0.) ) );
  // zcoils (in XY plane shifted in z)
  Coils.push_back( RectangularCoil( 100, zcoil_len, zcoil_len, TVector3(1.,0.,0.), TVector3(0.,1.,0.), TVector3(0.,0.,-zcoil_len/2) ) );
  Coils.push_back( RectangularCoil( 100, zcoil_len, zcoil_len, TVector3(1.,0.,0.), TVector3(0.,1.,0.), TVector3(0.,0., zcoil_len/2) ) );

  std::vector<std::string> coil_types = { "YZ", "YZ", "XZ", "XZ", "XY", "XY" };
  
  // plot coil locations
  plot_coils( Coils, coil_types );


  // coil currents
  std::vector< float > Currents{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

  // Biot-Savart calculators
  std::vector< BiotSavart > BSvec;
  for ( unsigned i=0; i<Currents.size(); ++i ){
    BSvec.push_back(  BiotSavart( Currents[i], Coils[i].get_wire() ) );
  }

  // generate field maps in XY, YZ, XZ planes
  TH2D * hXY = new TH2D("hXY","Btot (milli-gauss); X (m); Y (m)", 100, -zcoil_len/2, zcoil_len/2, 100, -zcoil_len/2, zcoil_len/2 );
  TH2D * hXZ = new TH2D("hXZ","Btot (milli-gauss); X (m); Z (m)", 100, -ycoil_len/2, ycoil_len/2, 100, -ycoil_len/2, ycoil_len/2 );
  TH2D * hYZ = new TH2D("hYZ","Btot (milli-gauss); Y (m); Z (m)", 100, -xcoil_len/2, xcoil_len/2, 100, -xcoil_len/2, xcoil_len/2 );

  {
    std::cout<<"Calculating for XY plane"<<std::endl;
    float z = 0.0;
    //XY
    for (unsigned ix=1; ix<=hXY->GetNbinsX(); ++ix ){
      float x = hXY->GetXaxis()->GetBinCenter( ix );
      for (unsigned iy=1; iy<=hXY->GetNbinsY(); ++iy ){
	float y = hXY->GetYaxis()->GetBinCenter( iy );
	TVector3 Btot(0.,0.,0.);
	for ( const BiotSavart& b : BSvec ){
	  Btot += b.get_B( TVector3( x, y, z ) );
	}
	hXY->SetBinContent( ix, iy, Btot.Mag() * T_to_mgauss );
      }
    }
  }
  
  {
    std::cout<<"Calculating for XZ plane"<<std::endl;
    float y = 0.0;
    //XZ
    for (unsigned ix=1; ix<=hXZ->GetNbinsX(); ++ix ){
      float x = hXZ->GetXaxis()->GetBinCenter( ix );
      for (unsigned iy=1; iy<=hXZ->GetNbinsY(); ++iy ){
	float z = hXZ->GetYaxis()->GetBinCenter( iy );
	TVector3 Btot(0.,0.,0.);
	for ( const BiotSavart& b : BSvec ){
	  Btot += b.get_B( TVector3( x, y, z ) );
	}
	hXZ->SetBinContent( ix, iy, Btot.Mag() * T_to_mgauss );
      }
    }
  }  
  
  {
    std::cout<<"Calculating for YZ plane"<<std::endl;
    float x = 0.0;
    //YZ
    for (unsigned ix=1; ix<=hYZ->GetNbinsX(); ++ix ){
      float y = hYZ->GetXaxis()->GetBinCenter( ix );
      for (unsigned iy=1; iy<=hYZ->GetNbinsY(); ++iy ){
	float z = hYZ->GetYaxis()->GetBinCenter( iy );
	TVector3 Btot(0.,0.,0.);
	for ( const BiotSavart& b : BSvec ){
	  Btot += b.get_B( TVector3( x, y, z ) );
	}
	hYZ->SetBinContent( ix, iy, Btot.Mag() * T_to_mgauss );
      }
    }
  }
  
  
  hXY->SetMaximum( 1000.0 );
  hXZ->SetMaximum( 1000.0 );
  hYZ->SetMaximum( 1000.0 );
  
  fout->Write();
  fout->Close();
}


int main(){

  biotsavart();
  

  return 0;
}
  

