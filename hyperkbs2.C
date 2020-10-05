/*
 * hyperkbs2.C
 * Program to calculate magnetic field around each PMT in hyperk.
 *
 * Parameterise positions of coils.  Using saddle coil for horizontal component
 * and some form of solenoid with more coils at top and bottom to even out
 * magnetic field.
 *
 * +z is vertical direction, +x is direction of B-field in horizontal plane
 * 
 * Author: Blair Jamieson (Mar. 2018)
 *
 * This version to calculate field from Kyla's calculation
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
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "MCMC.h"

using std::cout;
using std::endl;


// hk design report values:
//const std::vector< double > VertBodyZ = {
//  2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0,
//  30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0};

/// horizontally oriented coil Z locations (at outermost radius of cavern)
/// producing vertical magnetic Field
//const std::vector< double > VertBodyZ = {1.60851,3.21572,4.82029,6.42086,8.01594,9.60397,11.1832,12.7516,14.3071,15.847,17.3682,18.8671,20.3393,21.7792,23.1796,24.5311,25.8213,27.0326,28.1373,29.1031,29.7559};

// Values from Kyla:
const std::vector< double > VertBodyZ = {0.85398,1.70777,2.56117,3.41399,4.26601,5.11704,5.96687,6.81528,7.66204,8.50691,9.34965,10.19,11.0277,11.8624,12.6938,13.5215,14.3453,15.1647,15.9792,16.7884,17.5917,18.3886,19.1785,19.9605,20.7339,21.4976,22.2507,22.9921,23.7199,24.4329,25.1292,25.8058,26.4604,27.0899,27.6875,28.2517,28.7711,29.2429,29.5925,29.8746};

// Values after 10000 step MCMC
/*const std::vector< double > VertBodyZ = {
  1.2946153,
  2.4046319,
  2.8513065, 
  5.0431714, 
  5.541175,  
  5.6609456, 
  6.1858474, 
  6.5182402, 
  6.6520096, 
  7.5331453, 
  8.4744004,
  9.3058129,
  10.12976, 
  11.754781,
  12.122738,
  13.687004,
  13.758501,
  15.948217,
  16.190154,
  16.692468,
  17.375495,
  19.297744,
  19.483985,
  20.091876,
  20.44237, 
  21.232053,
  22.89082, 
  23.552511,
  24.631542,
  24.786075,
  25.006081,
  26.54446, 
  26.701242,
  27.765176,
  28.443082,
  28.540903,
  28.599587,
  28.827276,
  29.002775,
  29.142188 };
*/

/// horizontally oriented coil R locations (at top and bottom of cavern)
//const std::vector< double > VertEndR = {35.1422,34.3698,33.0804,31.3654,29.1741,26.4126,22.8925,18.1936,10.9021};


// Values from Kyla:
//const std::vector< double > VertEndR = {35.2685,34.9645,34.5561,33.9537,33.2445,32.4154,31.4608,30.3741,29.1437,27.7537,26.1811,24.3924,22.3379,19.9374,17.0471,13.3546,7.79596};
// Values after 10000 step MCMC
const std::vector< double > VertEndR = {
  36.313228, 
  34.93007 , 
  34.767463, 
  34.513694, 
  33.720215,
  33.327414, 
  31.442605,
  31.005992,
  29.34141 , 
  28.035282, 
  26.409661, 
  24.4938  , 
  22.659889, 
  20.028578, 
  17.125543, 
  13.522629,
  7.9255305 };

/// number of equally spaced vertically oriented coil X locations
/// producing horizontal magnetic Field
const unsigned HorzN = 18;

const int npar = 3 + VertBodyZ.size() + VertEndR.size();


const double pi = std::acos(-1.0);
const double mu0_over_4pi = 1.0e-7;

// Earth's field in this problem is (303, 0, -366  ) mgauss
// where +z is vertical direction, +x is direction of B-field in horizontal plane
const TVector3 BEarth{ 0.303e-4, 0.0, -0.366e-4 }; // Tesla
//const TVector3 BEarth{ 0.0, 0.0, -0.366e-4 };
//const TVector3 BEarth{ 0.0, 0.0, 0.0 };

// HyperK dimensions
const double hk_cavern_H = 60.0; //m
const double hk_cavern_R = 74.0/2.0; //m

const double hyperk_H = 54.800 / 2.0;// + 0.3; //m (half-height of HK from center to center-of PMT support structure)
const double hyperk_R = 70.800 / 2.0;// + 0.3; // m (center radius of PMT support structure)
const double hyperk_Rface = 70.800 / 2.0; // m (center radius of PMT cathode edges)
//const double hpyerk_circumference = hyperk_R * 2 * pi; // meters
const double hyperk_pmtsep = 0.703; // meters

//const unsigned npmts = 44000;

/// hk_pmtlocs
/// Function to calculate center location of support structure for each PMT
/// Returns vectors of TVector3: one for top, one for bottom, and one for sides
void hk_pmtlocs( std::vector< TVector3 >& toplocs,
		 std::vector< TVector3 >& botlocs,
		 std::vector< TVector3 >& sidelocs ){
  static bool firstcall=true;
  // calculate positions of PMTs?
  // want positions of center of support structure for that PMT.

  // number of pmts on top and bottom (each)
  // build up locations in plane
  for ( double fx =-hyperk_Rface + hyperk_pmtsep; fx < hyperk_Rface-hyperk_pmtsep/2; fx += hyperk_pmtsep ){ 
    for ( double fy =-hyperk_Rface + hyperk_pmtsep; fy < hyperk_Rface-hyperk_pmtsep/2; fy += hyperk_pmtsep ){ 
      double fr = std::sqrt( fx*fx + fy*fy );
      if ( fr < hyperk_Rface ){
	toplocs.push_back( TVector3{ fx, fy, hyperk_H } );
	botlocs.push_back( TVector3{ fx, fy, -hyperk_H } );
      }
    }
  }
  
  // number on each layer on sides:
  unsigned nperlayer = 2 * pi * hyperk_Rface / hyperk_pmtsep;

  for ( double zlayer = -(hyperk_H-0.3)+hyperk_pmtsep; zlayer < hyperk_H-0.3-hyperk_pmtsep/2; zlayer+=hyperk_pmtsep ){
    for ( unsigned ith = 0; ith < nperlayer; ++ith ){
      double theta = 2.0*pi*ith/double(nperlayer);
      sidelocs.push_back( TVector3{ hyperk_R*std::cos( theta ), hyperk_R*std::sin( theta ), zlayer } );
    }
  }

  // Final config
  std::cout<<"Final configuration:"
	   <<" ntop ="<<toplocs.size()
	   <<" nbot ="<<botlocs.size()
	   <<" nsides ="<<sidelocs.size()
	   <<" Total="<<toplocs.size()+botlocs.size()+sidelocs.size()
	   <<std::endl;

  // plot locations
  std::vector<double> xloc;
  std::vector<double> yloc;
  std::vector<double> zloc;
  for ( unsigned i=0; i<toplocs.size(); ++i) {
    xloc.push_back( toplocs[i].X() );
    yloc.push_back( toplocs[i].Y() );
    zloc.push_back( toplocs[i].Z() );
  }
  for ( unsigned i=0; i<botlocs.size(); ++i) {
    xloc.push_back( botlocs[i].X() );
    yloc.push_back( botlocs[i].Y() );
    zloc.push_back( botlocs[i].Z() );
  }
  for ( unsigned i=0; i<sidelocs.size(); ++i) {
    xloc.push_back( sidelocs[i].X() );
    yloc.push_back( sidelocs[i].Y() );
    zloc.push_back( sidelocs[i].Z() );
  }

  if (firstcall==true){
    firstcall=false;
    TGraph2D * tg2d = new TGraph2D( xloc.size(), &xloc[0], &yloc[0], &zloc[0] );
    tg2d->SetName("tg_pmtsupportlocs");
    tg2d->SetTitle("HK PMT Support Locations; X (m); Y(m); Z(m)");
    tg2d->Write();
  }

}



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

/// HKEndCapCoils
/// Vertical coils at top and bottom to try to keep field more uniform
/// Coils are rings in xy-plane with with different radius.
struct HKEndCapCoils{
  HKEndCapCoils(); // default 7 rings spaced at 2m, starting at 2m radius, all 1A
  HKEndCapCoils( double current, unsigned nrings );
  HKEndCapCoils( double current, const std::vector<double>& locs );//radius locations of coils
  
  Wire get_wire() const { return fw; }
  double get_current() const { return fI; }
  double  get_penalty() const { return fPenalty; }
  friend void plot_coils( const HKEndCapCoils& hk );
protected:
  void Init();
  double fI;
  unsigned fnrings;
  double fPenalty;
  Wire fw;
};

HKEndCapCoils::HKEndCapCoils(){
  fnrings=7;
  fI=1.0;
  fPenalty = 0.0;
  Init();
}


void plot_coils( const HKEndCapCoils& hk ) {
  std::vector<double> x,y,z;
  for ( const WireElement& we : hk.fw ){
    TVector3 loc = we.get_rp();
    x.push_back( loc.X() );
    y.push_back( loc.Y() );
    z.push_back( loc.Z() );
  }
  TGraph2D * tg2d = new TGraph2D( x.size(), &x[0], &y[0], &z[0] );
  tg2d->SetName("tg_vertcoilloc");
  tg2d->SetTitle("Veritcal Coil Locations; X(m); Y(m); Z(m)");
  tg2d->Write();
}


HKEndCapCoils::HKEndCapCoils( double current, unsigned nrings) :
  fI(current), fnrings(nrings), fPenalty(0.0) {
  Init();
}

void HKEndCapCoils::Init(){
  const int nsteps = 360; // break into nsteps around each coil
  const double dphi = 2.0*pi/nsteps; // radians to step by

  // Add end cap coils
  double curz = hk_cavern_H/2.0;
  for ( double coilR = 2.0; coilR<hk_cavern_R-1.0; coilR += 2.0 ){
    double dl = 2.0 * pi * coilR / nsteps; // size of each step in meters
    for ( int sign = -1; sign<=2; sign += 2 ){
      curz = sign * fabs( curz );
      for ( int step=0; step<nsteps; ++step ){
	double phi = step*dphi;
	double curx = coilR * std::cos( phi );
	double cury = coilR * std::sin( phi );
	TVector3 stepdir( -std::sin(phi)*dl, std::cos( phi )*dl, 0.0 );
	TVector3 steploc( curx, cury, curz );
	fw.push_back( WireElement( stepdir, steploc ) );
      }
    }    
  }
}

HKEndCapCoils::HKEndCapCoils( double current, const std::vector<double>& locs ) :
  fI( current ), fnrings(locs.size()), fPenalty(0.0) {
  const int nsteps = 360; // break into nsteps around each coil
  const double dphi = 2.0*pi/nsteps; // radians to step by

  for ( unsigned i=0; i<locs.size()-1; ++i ){
    if ( locs[i+1] > locs[i] ) {
      fPenalty += (locs[i]-locs[i+1])*(locs[i]-locs[i+1])*1000000.0;
    }
  }    

  // Add end cap coils
  double curz = hk_cavern_H/2.0;
  for ( std::vector<double>::const_iterator it = locs.begin(); it != locs.end(); ++it ){
    double coilR = *it;
    double dl = 2.0 * pi * coilR / nsteps; // size of each step in meters
    for ( int sign = -1; sign<=2; sign += 2 ){
      curz = sign * fabs( curz );
      for ( int step=0; step<nsteps; ++step ){
	double phi = step*dphi;
	double curx = coilR * std::cos( phi );
	double cury = coilR * std::sin( phi );
	TVector3 stepdir( -std::sin(phi)*dl, std::cos( phi )*dl, 0.0 );
	TVector3 steploc( curx, cury, curz );
	fw.push_back( WireElement( stepdir, steploc ) );
      }
    }    
  }  
}


/// HKVerticalCoils
/// Vertical coils to compensate vertical portion of HK magnetic field
/// Coils are rings in xy-plane at different spacings along z-axis.
/// Assumed to be symmetric about z=0.
/// Specify spacings by std::vector<double> that gives number of coils, and their spacing
struct HKVerticalCoils {
  HKVerticalCoils(); // default all equal spacing coils and 30 coils total, all with 1A current
  HKVerticalCoils( double current, const std::vector<double>& top_half_locs );

  Wire get_wire() const { return fw; }
  double get_current() const { return fI; }
  double get_penalty() const;  
  friend void plot_coils( const HKVerticalCoils& hk );
  friend std::ostream& operator<<( std::ostream& os, const HKVerticalCoils& hk );
protected:
  void Init();
  Wire fw;
  double fI;
  std::vector<double> locs;
};

double HKVerticalCoils::get_penalty() const {
  double penalty=0.0;
  for( unsigned i=0; i<locs.size(); ++i ) {
    if ( locs[i] < 0.0 ) penalty= locs[i]*locs[i]*1000000.0;
    if ( i+1 < locs.size() && locs[i+1] < locs[i] ) {
      penalty += (locs[i] - locs[i+1] )*(locs[i] - locs[i+1] )*1000000.0;
    }
    if ( fabs( locs[i] ) > hk_cavern_H/2.0 ) {
      penalty += (locs[i] - hk_cavern_H/2.0)*(locs[i] - hk_cavern_H/2.0)*1000000.0;
    }
  }
  return penalty;
}


std::ostream& operator<<( std::ostream& os, const HKVerticalCoils& hk ){
  os << "HKVerticalCoils Nlocs="<<hk.locs.size()<<" I="<<hk.fI<<" Nwire-elemets="<<hk.fw.size()<<std::endl;

  TVector3 btot(0.,0.,0.);
  for (unsigned i=0; i<hk.fw.size(); ++i){
    WireElement we = hk.fw[i];
    TVector3 dir = we.get_dl();
    TVector3 loc = we.get_rp();
    TVector3 belem(0.,0.,0.);
    we.get_dl_x_r3( TVector3{0.,0.,0.}, belem );
    belem *= mu0_over_4pi * hk.fI;
    btot+=belem;
    os << std::setprecision( 4 );
    os << "  (x,y,z) = ("
       << std::setw(8) << loc.X() <<", "
       << std::setw(8) << loc.Y()<<", "
       << std::setw(8)<<loc.Z()<<")"
       << "  (dx,dy,dz) = ("
       << std::setw(8) << dir.X() <<", "
       << std::setw(8) << dir.Y()<<", "
       << std::setw(8) << dir.Z()<<")"
       <<"    (dBx,dBy,dBz)_0 = ("
       << std::setw(8) << belem.X() <<", "
       << std::setw(8) << belem.Y()<<", "
       << std::setw(8) << belem.Z()<<")"	 
       <<std::endl;
  }
  os<<" (Bx,By,Bz)_0 = ("
    << std::setw(8) << btot.X() <<", "
    << std::setw(8) << btot.Y()<<", "
    << std::setw(8) << btot.Z()<<")"	 
    <<std::endl;
  

  return os;
}

void plot_coils( const HKVerticalCoils& hk ) {
  std::vector<double> x,y,z;
  for ( const WireElement& we : hk.fw ){
    TVector3 loc = we.get_rp();
    x.push_back( loc.X() );
    y.push_back( loc.Y() );
    z.push_back( loc.Z() );
  }
  TGraph2D * tg2d = new TGraph2D( x.size(), &x[0], &y[0], &z[0] );
  tg2d->SetName("tg_vertcoilloc");
  tg2d->SetTitle("Veritcal Coil Locations; X(m); Y(m); Z(m)");
  tg2d->Write();
}



HKVerticalCoils::HKVerticalCoils(){
  // equally spaced 30 coils in total:
  double coilsep = hk_cavern_H/30.0;
  for ( unsigned i=1; i<16; ++i ) locs.push_back( coilsep*i );
  fI = 1.0;//Amp
  Init();
}

HKVerticalCoils::HKVerticalCoils( double current, const std::vector<double>& top_half_locs ) :
  fI(current), locs( top_half_locs ) {
  Init();
}

void HKVerticalCoils::Init(){
  const int nsteps = 360; // break into nsteps around each coil
  const double dphi = 2.0*pi/nsteps; // radians to step by

  double curz = 0.0;
  double dl = 2.0 * pi * hk_cavern_R / nsteps; // size of each step in meters
  for ( unsigned icoil = 0; icoil <= locs.size(); ++icoil ){
    if ( icoil == locs.size() ) curz = 0.0;  // add coil at zero
    else curz = locs[icoil];
    for ( int sign = -1; sign<=2; sign += 2 ){
      if ( icoil == locs.size() && sign != -1 ) continue; // add coil at zero
      curz = sign * fabs( curz );
      for ( int step=0; step<nsteps; ++step ){
	double phi = step*dphi;
	double curx = hk_cavern_R * std::cos( phi );
	double cury = hk_cavern_R * std::sin( phi );
	TVector3 stepdir( -std::sin(phi)*dl, std::cos( phi )*dl, 0.0 );
	TVector3 steploc( curx, cury, curz );
	fw.push_back( WireElement( stepdir, steploc ) );
      }
    }
  }
}

/// HKHorizontalCoils
/// Horizontal coils compensate horizontal component of magnetic field.
/// Based on sine phi coil
/// See spacing of coils from : C.P. Bidinosti, Journal of Magnetic Resonance 177 (2005) 31.
/// Pick number of coils in quadrant, rest is fixed.
/// Default is 15 coils in each quadrant and 1Amp.
struct HKHorizontalCoils{
  HKHorizontalCoils();
  HKHorizontalCoils(  double current, unsigned ncoils );

  unsigned get_ncoils() const { return fN; }
  double   get_current() const { return fI; }
  Wire     get_wire() const { return fw; }

  friend void plot_coils( const HKHorizontalCoils& hk);
  friend std::ostream& operator<<( std::ostream& os, const HKHorizontalCoils& hk );
private:
  void Init();
  Wire fw;
  double fI;
  unsigned fN;
};

std::ostream& operator<<( std::ostream& os, const HKHorizontalCoils& hk ){
  os << "HKHorizontalCoils Ncoils="<<hk.fN<<" I="<<hk.fI<<" Nwire-elemets="<<hk.fw.size()<<std::endl;
  TVector3 btot(0.,0.,0.);
  for (unsigned i=0; i<hk.fw.size(); ++i){
    WireElement we = hk.fw[i];
    TVector3 dir = we.get_dl();
    TVector3 loc = we.get_rp();
    TVector3 belem(0.,0.,0.);
    we.get_dl_x_r3( TVector3{0.,0.,0.}, belem );
    belem *= mu0_over_4pi * hk.fI;
    btot+=belem;
    os << std::setprecision( 4 );
    os << "  (x,y,z) = ("
       << std::setw(8) << loc.X() <<", "
       << std::setw(8) << loc.Y()<<", "
       << std::setw(8) << loc.Z()<<")"
       << "  (dx,dy,dz) = ("
       << std::setw(8) << dir.X() <<", "
       << std::setw(8) << dir.Y()<<", "
       << std::setw(8) << dir.Z()<<")"
       <<"    (dBx,dBy,dBz)_0 = ("
       << std::setw(8) << belem.X() <<", "
       << std::setw(8) << belem.Y()<<", "
       << std::setw(8) << belem.Z()<<")"	 
       <<std::endl;
  }
  os<<" (Bx,By,Bz)_0 = ("
    << std::setw(8) << btot.X() <<", "
    << std::setw(8) << btot.Y()<<", "
    << std::setw(8) << btot.Z()<<")"	 
    <<std::endl;
  

  return os;
}

void plot_coils( const HKHorizontalCoils& hk ) {
  std::vector<double> x,y,z;
  for ( const WireElement& we : hk.fw ){
    TVector3 loc = we.get_rp();
    x.push_back( loc.X() );
    y.push_back( loc.Y() );
    z.push_back( loc.Z() );
  }
  TGraph2D * tg2d = new TGraph2D( x.size(), &x[0], &y[0], &z[0] );
  tg2d->SetName("tg_horizcoilloc");
  tg2d->SetTitle("Horizontal Coil Locations; X(m); Y(m); Z(m)");
  tg2d->Write();
}

HKHorizontalCoils::HKHorizontalCoils(){
  fI = 1.0; //Amp
  fN = 15;
  Init();
}

HKHorizontalCoils::HKHorizontalCoils(  double current, unsigned ncoils ):
  fI( current ), fN( ncoils ) {
  Init();
}

void HKHorizontalCoils::Init(){
  bool use_straight_edge = true;
  const int nsteps = 120;// number of steps along each vertical wire
  for (unsigned iwire=0; iwire<fN; ++iwire){ // number of pairs of loops
    for ( int xsign=-1; xsign<2; xsign += 2 ){ // one loop at each side in x
      for ( int ysign=-1; ysign<2; ysign += 2 ){ // each loop has two halves
	double dz = hk_cavern_H/nsteps; // reset dz for each wire.
	double phij= std::acos( 1.0 - (2*iwire+1)/double(2*fN) );
	//std::cout<<"iwire="<<iwire<<" xsign="<<xsign<<" ysign="<<ysign
	//	 <<" phiij="<<phij<<"rad"<<std::endl;
	double curx = xsign * hk_cavern_R * std::cos( phij );
	double cury = ysign * hk_cavern_R * std::sin( phij );
	TVector3 dir{ 0., 0., -ysign*dz };
	//std::cout<<" dir="<<-ysign*dz<<std::endl;
	double curz;
	for ( int istep=0; istep<nsteps; ++istep ){ // vertical part 
	  curz = ysign*(hk_cavern_H/2 - dz*(float(istep)+0.5));
	  TVector3 loc{ curx, cury, curz };
	  //std::cout<<"HKHorizontal(1): add dir="<<dir.X()<<","<<dir.Y()<<","<<dir.Z()
	  //	   <<" loc="<<loc.X()<<","<<loc.Y()<<","<<loc.Z()<<std::endl;
	  fw.push_back( WireElement( dir, loc ) );
	}
	if ( use_straight_edge ){ // across top
	  int nstepsy = 2 * fabs(cury)/dz;
	  double dy = 2 * fabs(cury) / nstepsy; // retune dz due to int roundoff
	  curz = -ysign * hk_cavern_H / 2.0;
	  TVector3 dir2 = TVector3{ 0., -ysign*dy, 0. };
	  //std::cout<<" dir="<<-ysign*dz<<std::endl;
	  for ( int istep=0; istep<nstepsy; ++istep) {// horizontal part
	    double yloc = cury - ysign * dy * (float(istep)+0.5);
	    TVector3 loc{ curx, yloc, curz };
	    //std::cout<<"HKHorizontal(2): add dir="<<dir2.X()<<","<<dir2.Y()<<","<<dir2.Z()
	    //	   <<" loc="<<loc.X()<<","<<loc.Y()<<","<<loc.Z()<<std::endl;
	    fw.push_back( WireElement( dir2, loc ) );
	  } 
	} else { // around edge of top
	  int nstepsy = 2 * fabs(cury)/dz;
	  double dy = 2 * fabs(cury) / nstepsy; // retune dz due to int roundoff
	  curz = -ysign * hk_cavern_H / 2.0;
	  //std::cout<<" dir="<<-ysign*dz<<std::endl;
	  double prevx = curx;
	  for ( int istep=0; istep<nstepsy; ++istep) {// horizontal part
	    double yloc = cury - ysign * dy * (float(istep)+0.5);
	    double xloc = std::sqrt( hk_cavern_R*hk_cavern_R - yloc*yloc );
	    if ( curx < 0 ) xloc = -xloc;
	    TVector3 dir2 = TVector3{ xloc-prevx, -ysign*dy, 0. };
	    TVector3 loc{ xloc, yloc, curz };
	    //std::cout<<"HKHorizontal(2): add dir="<<dir2.X()<<","<<dir2.Y()<<","<<dir2.Z()
	    //	   <<" loc="<<loc.X()<<","<<loc.Y()<<","<<loc.Z()<<std::endl;
	    fw.push_back( WireElement( dir2, loc ) );
	    prevx = xloc;
	  } 
	}
      }
    }
  }
}




/// The chi2 function to minimize.
double BFieldChi2( const double * pars ){
  static unsigned ncalls=0;
  static std::vector< TVector3 > pmtlocs;
  static std::vector< TVector3 > toplocs, botlocs;
  if ( ncalls==0 ){
    hk_pmtlocs( toplocs, botlocs, pmtlocs );
  }

  double Ih = pars[0]; // current for horizontal compensation
  double Iv = pars[1]; // current for vertical compensation
  double Iend = pars[2]; // current in end cap coils
  for ( unsigned ipar=0; ipar<npar; ++ipar ){
    std::cout<<"par["<<ipar<<"]="<<pars[ipar]<<"  ";
  }
  std::cout<<std::endl;
  std::vector<double> vertlocs, endlocs;
  for (unsigned i=0; i<VertBodyZ.size(); ++i ){
    vertlocs.push_back( pars[3+i] );
  }
  for (unsigned i=0; i<VertEndR.size(); ++i ){
    endlocs.push_back( pars[3+VertBodyZ.size()+i] );
  }
  
  HKHorizontalCoils horzcoils( Ih, HorzN );
  HKVerticalCoils   vertcoils( Iv, vertlocs );
  HKEndCapCoils     endcoils( Iend, endlocs );

  BiotSavart horzbs( Ih, horzcoils.get_wire() );
  BiotSavart vertbs( Iv, vertcoils.get_wire() );
  BiotSavart endbs( Iend, endcoils.get_wire() );
  
  double chi2=0.0;

  // Put a very strong constraint requiring the coil spacings to add
  // up to the height of hk
  double spacing_penalty = vertcoils.get_penalty();
  chi2 += spacing_penalty;
  spacing_penalty = endcoils.get_penalty();
  chi2 += spacing_penalty;
  
  // Now add contribution from each PMT location
  TVector3 down( 0., 0., -1.0 );
  for ( const TVector3& loc : toplocs ){
    TVector3 btot = horzbs.get_B( loc ) + vertbs.get_B( loc ) + endbs.get_B( loc ) + BEarth;
    //double delta_chi2 = pow( btot * down, 2 ) * 4.0e10; // only component perp. to PMT (down)
    TVector3 bperp = btot.Cross( down );
    double delta_chi2 = pow( bperp.Mag(), 2 ) * 4.0e10; // only component perp. to PMT (down)
    chi2 += delta_chi2;
  }
  TVector3 up( 0., 0., 1.0 );
  for ( const TVector3& loc : toplocs ){
    TVector3 btot = horzbs.get_B( loc ) + vertbs.get_B( loc ) + endbs.get_B( loc ) + BEarth;
    //double delta_chi2 = pow( btot * up, 2 ) * 4.0e10; // only component perp. to PMT (up)
    TVector3 bperp = btot.Cross( up );
    double delta_chi2 = pow( bperp.Mag(), 2 ) * 4.0e10; // only component perp. to PMT (up)
    chi2 += delta_chi2;
  }
  for ( const TVector3& loc : toplocs ){
    TVector3 btot = horzbs.get_B( loc ) + vertbs.get_B( loc ) + endbs.get_B( loc ) + BEarth;
    TVector3 side = loc.Unit();
    TVector3 bperp = btot.Cross( side );
    //double delta_chi2 = pow( btot * side.Unit(), 2 ) * 4.0e10; // only component perp. to PMT (side))
    double delta_chi2 = pow( bperp.Mag(), 2 ) * 4.0e10; // only component perp. to PMT (side))
    chi2 += delta_chi2;
  }

  std::cout<<"BFieldChi2: "
	   <<" spacing_penalty="<< spacing_penalty
	   <<" chi2="<<std::setprecision(8)<<chi2
	   <<" ncalls="<<++ncalls<<std::endl;
  
  return chi2;
}

void MCMCFcn( int &npar, double *gin, double &chi2, double *par, int iflag ){
  chi2 = BFieldChi2( par );
}

  
/// make_histograms
/// Function to build a few simple histograms for 
/// Converts all magnetic fields to gauss, and adds earth's field
void make_histograms( const std::vector<double>&  pars, const std::string & tag ){
  // pars 0 to 2 are currents Ih, Iv, Iend
  // par 3 to 2+VertBodyZ.size() are horizontally oriented coil locations in Z
  // par 3+VertBodyZ.size() to 2+VertBodyZ.size()+VertEndR.size() are horizontal end cap coil locations in R
  bool verbose = true;

  
  double Ih = pars[0]; // current for horizontal compensation
  double Iv = pars[1]; // current for vertical compensation
  double Iend = pars[2]; // current in end cap coils
  for ( unsigned ipar=0; ipar<npar; ++ipar ){
    std::cout<<"par["<<ipar<<"]="<<pars[ipar]<<"  ";
  }
  std::cout<<std::endl;
  std::vector<double> vertlocs, endlocs;
  for (unsigned i=0; i<VertBodyZ.size(); ++i ){
    vertlocs.push_back( pars[3+i] );
  }
  for (unsigned i=0; i<VertEndR.size(); ++i ){
    endlocs.push_back( pars[3+VertBodyZ.size()+i] );
  }

  if (verbose) std::cout<<"Building coils"<<std::endl;
  HKHorizontalCoils horzcoils( Ih, HorzN );
  HKVerticalCoils   vertcoils( Iv, vertlocs );
  HKEndCapCoils     endcoils( Iend, endlocs );


  if (verbose) std::cout<<"Building BiotSavart calculation"<<std::endl;
  BiotSavart horzbs( Ih, horzcoils.get_wire() );
  BiotSavart vertbs( Iv, vertcoils.get_wire() );
  BiotSavart endbs( Iend, endcoils.get_wire() );

  if (verbose) std::cout<<"Plotting coil locations"<<std::endl;
  plot_coils( horzcoils );
  plot_coils( vertcoils );
  plot_coils( endcoils );

  if (verbose) std::cout<<"Getting PMT locations"<<std::endl;
  std::vector< TVector3 > sidelocs;
  std::vector< TVector3 > toplocs, botlocs;
  hk_pmtlocs( toplocs, botlocs, sidelocs );

  if (verbose) std::cout<<"Building 1D Bfield histograms"<<std::endl;
  std::ostringstream os;
  os.str(""); os.clear(); os << tag << "_Bperp";
  TH1D* hBpmt_perp = new TH1D( os.str().c_str(), " ;B_{perp} (gauss); PMTs/bin", 200, 0., 1.0 );
  os.str(""); os.clear(); os << tag << "_Bmag";
  TH1D* hBpmt = new TH1D( os.str().c_str(), " ;B_{mag} (gauss); PMTs/bin", 200, 0., 1.0 );
  os.str(""); os.clear(); os << tag << "_Bx";
  TH1D* hBx = new TH1D( os.str().c_str(), " ;B_{x} (gauss); PMTs/bin", 200, -1.0, 1.0 );
  os.str(""); os.clear(); os << tag << "_By";
  TH1D* hBy = new TH1D( os.str().c_str(), " ;B_{y} (gauss); PMTs/bin", 200, -1.0, 1.0 );
  os.str(""); os.clear(); os << tag << "_Bz";
  TH1D* hBz = new TH1D( os.str().c_str(), " ;B_{z} (gauss); PMTs/bin", 200, -1.0, 1.0 );
  
  // main histogram to make is of the magnetic field at each PMT
  // then plot Bfield accross face of PMTs with highest Bfield
  double bmax=0.0;
  TVector3 locbmax{0.,0.,0.};
  if (verbose) std::cout<<"Filling 1D Bfiled for toplocs"<<std::endl;
  TVector3 down( 0., 0., -1.0 );
  for ( unsigned ipmt=0; ipmt<toplocs.size(); ++ipmt ){
    TVector3 loc =  toplocs[ipmt];
    TVector3 btot = horzbs.get_B( loc ) + vertbs.get_B( loc ) + endbs.get_B( loc );
    //std::cout<<"ipmt "<<ipmt<<" btot="<<btot.X()<<","<<btot.Y()<<","<<btot.Z()
    //	     <<" bearth="<<BEarth.X()<<","<<BEarth.Y()<<","<<BEarth.Z()<<std::endl;
    btot += BEarth;
    btot *= 1.0e4;//convert to gauss
    hBpmt->Fill( btot.Mag() );
    hBx->Fill( btot.X() );
    hBy->Fill( btot.Y() );
    hBz->Fill( btot.Z() );
    if ( btot.Mag() > bmax ){
      locbmax = loc;
      bmax = btot.Mag();
    }
    TVector3 bperp = btot.Cross( down );
    hBpmt_perp->Fill( bperp.Mag() );
  }
  if (verbose) std::cout<<"Filling 1D Bfiled for botlocs"<<std::endl;
  TVector3 up( 0., 0., 1.0 );
  for ( unsigned ipmt=0; ipmt<botlocs.size(); ++ipmt ){
    TVector3 loc =  botlocs[ipmt];
    TVector3 btot = horzbs.get_B( loc ) + vertbs.get_B( loc ) + endbs.get_B( loc );
    //std::cout<<"ipmt "<<ipmt<<" btot="<<btot.X()<<","<<btot.Y()<<","<<btot.Z()
    //	     <<" bearth="<<BEarth.X()<<","<<BEarth.Y()<<","<<BEarth.Z()<<std::endl;
    btot += BEarth;
    btot *= 1.0e4;//convert to gauss
    hBpmt->Fill( btot.Mag() );
    hBx->Fill( btot.X() );
    hBy->Fill( btot.Y() );
    hBz->Fill( btot.Z() );
    if ( btot.Mag() > bmax ){
      locbmax = loc;
      bmax = btot.Mag();
    }
    TVector3 bperp = btot.Cross( up );
    hBpmt_perp->Fill( bperp.Mag() );
  }
  if (verbose) std::cout<<"Filling 1D Bfiled for sidelocs"<<std::endl;
  for ( unsigned ipmt=0; ipmt<sidelocs.size(); ++ipmt ){
    TVector3 loc =  sidelocs[ipmt];
    TVector3 btot = horzbs.get_B( loc ) + vertbs.get_B( loc ) + endbs.get_B( loc );
    //std::cout<<"ipmt "<<ipmt<<" btot="<<btot.X()<<","<<btot.Y()<<","<<btot.Z()
    //	     <<" bearth="<<BEarth.X()<<","<<BEarth.Y()<<","<<BEarth.Z()<<std::endl;
    btot += BEarth;
    btot *= 1.0e4;//convert to gauss
    hBpmt->Fill( btot.Mag() );
    hBx->Fill( btot.X() );
    hBy->Fill( btot.Y() );
    hBz->Fill( btot.Z() );
    if ( btot.Mag() > bmax ){
      locbmax = loc;
      bmax = btot.Mag();
    }
    TVector3 side = loc.Unit();
    TVector3 bperp = btot.Cross( side );
    hBpmt_perp->Fill( bperp.Mag() );
  }


  // make histograms of magnetic field around detecto
  if (verbose) std::cout<<"Making 2D Bfield plots in XY at Z=0 "<<std::endl;
  os.str(""); os.clear(); os << tag << "_hxyBx";
  TH2D* hxyBx = new TH2D(os.str().c_str()," Bx (gauss); x (m); y (m) ",100,-hk_cavern_R, hk_cavern_R,100,-hk_cavern_R,hk_cavern_R);
  os.str(""); os.clear(); os << tag << "_hxyBy";
  TH2D* hxyBy = new TH2D(os.str().c_str()," By (gauss); x (m); y (m) ",100,-hk_cavern_R, hk_cavern_R,100,-hk_cavern_R,hk_cavern_R);
  os.str(""); os.clear(); os << tag << "_hxyBz";
  TH2D* hxyBz = new TH2D(os.str().c_str()," Bz (gauss); x (m); y (m) ",100,-hk_cavern_R, hk_cavern_R,100,-hk_cavern_R,hk_cavern_R);
  os.str(""); os.clear(); os << tag << "_hxyB";
  TH2D* hxyB = new TH2D(os.str().c_str()," B (gauss); x (m); y (m) ",100,-hk_cavern_R, hk_cavern_R,100,-hk_cavern_R,hk_cavern_R);

  for ( int ix=1; ix<=100; ++ix ){
    for ( int iy=1; iy<=100; ++iy ){
      double x = hxyBx->GetXaxis()->GetBinCenter( ix );
      double y = hxyBx->GetYaxis()->GetBinCenter( iy );
      TVector3 loc = TVector3{ x, y, 0.0 } ;
      TVector3 btot = horzbs.get_B( loc ) + vertbs.get_B( loc ) + endbs.get_B( loc ) + BEarth;
      btot *= 1.0e4;//convert to gauss
      hxyBx->SetBinContent( ix, iy, btot.X() );
      hxyBy->SetBinContent( ix, iy, btot.Y() );
      hxyBz->SetBinContent( ix, iy, btot.Z() );
      hxyB->SetBinContent( ix, iy, btot.Mag() );
    }
  }

  if (verbose) std::cout<<"Making 2D Bfield plots in YZ at X=0 "<<std::endl;
  os.str(""); os.clear(); os << tag << "_hyzBx";
  TH2D* hyzBx = new TH2D(os.str().c_str()," Bx (gauss); y (m); z (m) ",100,-hk_cavern_R,hk_cavern_R,100,-hk_cavern_H/2, hk_cavern_H/2 );
  os.str(""); os.clear(); os << tag << "_hyzBy";
  TH2D* hyzBy = new TH2D(os.str().c_str()," By (gauss); y (m); z (m) ",100,-hk_cavern_R,hk_cavern_R,100,-hk_cavern_H/2, hk_cavern_H/2 );
  os.str(""); os.clear(); os << tag << "_hyzBz";
  TH2D* hyzBz = new TH2D(os.str().c_str()," Bz (gauss); y (m); z (m) ",100,-hk_cavern_R,hk_cavern_R,100,-hk_cavern_H/2, hk_cavern_H/2 );
  os.str(""); os.clear(); os << tag << "_hyzB";
  TH2D* hyzB = new TH2D(os.str().c_str()," B (gauss); y (m); z (m) ",100,-hk_cavern_R,hk_cavern_R,100,-hk_cavern_H/2, hk_cavern_H/2 );

  for ( int iy=1; iy<=100; ++iy ){
    for ( int iz=1; iz<=100; ++iz ){
      double y = hyzBx->GetXaxis()->GetBinCenter( iy );
      double z = hyzBx->GetYaxis()->GetBinCenter( iz );
      TVector3 loc = TVector3{ 0.0, y, z } ;
      TVector3 btot = horzbs.get_B( loc ) + vertbs.get_B( loc ) + endbs.get_B( loc ) + BEarth;
      btot *= 1.0e4;//convert to gauss

      hyzBx->SetBinContent( iy, iz, btot.X() );
      hyzBy->SetBinContent( iy, iz, btot.Y() );
      hyzBz->SetBinContent( iy, iz, btot.Z() );
      hyzB->SetBinContent( iy, iz, btot.Mag() );
    }
  }

  hxyBx->SetMinimum(-0.6); hxyBx->SetMaximum(0.6);
  hxyBy->SetMinimum(-0.6); hxyBy->SetMaximum(0.6);
  hxyBz->SetMinimum(-0.6); hxyBz->SetMaximum(0.6);
  hxyB->SetMinimum(0.0);  hxyB->SetMaximum(1.0);

  hyzBx->SetMinimum(-0.6); hyzBx->SetMaximum(0.6);
  hyzBy->SetMinimum(-0.6); hyzBy->SetMaximum(0.6);
  hyzBz->SetMinimum(-0.6); hyzBz->SetMaximum(0.6);
  hyzB->SetMinimum(0.0); hyzBz->SetMaximum(1.0);


  if (verbose) std::cout<<"Making 2D Bfield plots in XY at Z=top PMT locations "<<std::endl;
  os.str(""); os.clear(); os << tag << "_hxy_top_B";
  TH2D* hxytopB = new TH2D(os.str().c_str()," B (gauss); x (m); y (m) ",100,-hyperk_R,hyperk_R,100,-hyperk_R,hyperk_R);
  os.str(""); os.clear(); os << tag << "_hxy_bot_B";
  TH2D* hxybotB = new TH2D(os.str().c_str()," B (gauss); x (m); y (m) ",100,-hyperk_R,hyperk_R,100,-hyperk_R,hyperk_R);

  os.str(""); os.clear(); os << tag << "_hxy_top_Bperp";
  TH2D* hxytopBp = new TH2D(os.str().c_str()," Bperp (gauss); x (m); y (m) ",100,-hyperk_R,hyperk_R,100,-hyperk_R,hyperk_R);
  os.str(""); os.clear(); os << tag << "_hxy_bot_Bperp";
  TH2D* hxybotBp = new TH2D(os.str().c_str()," Bperp (gauss); x (m); y (m) ",100,-hyperk_R,hyperk_R,100,-hyperk_R,hyperk_R);

  os.str(""); os.clear(); os << tag << "_hxy_top_Bx";
  TH2D* hxytopBx = new TH2D(os.str().c_str()," Bx (gauss); x (m); y (m) ",100,-hyperk_R,hyperk_R,100,-hyperk_R,hyperk_R);
  os.str(""); os.clear(); os << tag << "_hxy_bot_Bx";
  TH2D* hxybotBx = new TH2D(os.str().c_str()," Bx (gauss); x (m); y (m) ",100,-hyperk_R,hyperk_R,100,-hyperk_R,hyperk_R);

  os.str(""); os.clear(); os << tag << "_hxy_top_By";
  TH2D* hxytopBy = new TH2D(os.str().c_str()," By (gauss); x (m); y (m) ",100,-hyperk_R,hyperk_R,100,-hyperk_R,hyperk_R);
  os.str(""); os.clear(); os << tag << "_hxy_bot_By";
  TH2D* hxybotBy = new TH2D(os.str().c_str()," By (gauss); x (m); y (m) ",100,-hyperk_R,hyperk_R,100,-hyperk_R,hyperk_R);

  os.str(""); os.clear(); os << tag << "_hxy_top_Bz";
  TH2D* hxytopBz = new TH2D(os.str().c_str()," B (gauss); x (m); y (m) ",100,-hyperk_R,hyperk_R,100,-hyperk_R,hyperk_R);
  os.str(""); os.clear(); os << tag << "_hxy_bot_Bz";
  TH2D* hxybotBz = new TH2D(os.str().c_str()," Bz (gauss); x (m); y (m) ",100,-hyperk_R,hyperk_R,100,-hyperk_R,hyperk_R);


  for ( int ix=1; ix<=100; ++ix ){
    for ( int iy=1; iy<=100; ++iy ){
      {
      double x = hxytopB->GetXaxis()->GetBinCenter( ix );
      double y = hxytopB->GetYaxis()->GetBinCenter( iy );
      TVector3 loc = TVector3{x, y, hyperk_H } ;
      TVector3 btot = horzbs.get_B( loc ) + vertbs.get_B( loc ) + endbs.get_B( loc ) + BEarth;
      btot *= 1.0e4;//convert to gauss
      hxytopB->SetBinContent( ix, iy, btot.Mag() );
      hxytopBx->SetBinContent( ix, iy, btot.X() );
      hxytopBy->SetBinContent( ix, iy, btot.Y() );
      hxytopBz->SetBinContent( ix, iy, btot.Z() );
      TVector3 bperp = btot.Cross( up );
      hxytopBp->SetBinContent( ix, iy, bperp.Mag() );
      
      }
      {
      double x = hxytopB->GetXaxis()->GetBinCenter( ix );
      double y = hxytopB->GetYaxis()->GetBinCenter( iy );
      TVector3 loc = TVector3{x, y, -hyperk_H } ;
      TVector3 btot = horzbs.get_B( loc ) + vertbs.get_B( loc ) + endbs.get_B( loc ) + BEarth;
      btot *= 1.0e4;//convert to gauss
      hxybotB->SetBinContent( ix, iy, btot.Mag() );
      hxybotBx->SetBinContent( ix, iy, btot.X() );
      hxybotBy->SetBinContent( ix, iy, btot.Y() );
      hxybotBz->SetBinContent( ix, iy, btot.Z() );
      TVector3 bperp = btot.Cross( up );
      hxybotBp->SetBinContent( ix, iy, bperp.Mag() );

      }
    }
  }

  if (verbose) std::cout<<"Making 2D Bfield plots in phi at PMT locations "<<std::endl;
  os.str(""); os.clear(); os << tag << "_hphiz_side_B";
  TH2D* hphizsideB = new TH2D(os.str().c_str()," B (gauss); #phi (rad); z (m) ",360,0.,2.0*pi, 100, -hyperk_H, hyperk_H);
  os.str(""); os.clear(); os << tag << "_hphiz_side_Bx";
  TH2D* hphizsideBx = new TH2D(os.str().c_str()," Bx (gauss); #phi (rad); z (m) ",360,0.,2.0*pi, 100, -hyperk_H, hyperk_H);
  os.str(""); os.clear(); os << tag << "_hphiz_side_By";
  TH2D* hphizsideBy = new TH2D(os.str().c_str()," By (gauss); #phi (rad); z (m) ",360,0.,2.0*pi, 100, -hyperk_H, hyperk_H);
  os.str(""); os.clear(); os << tag << "_hphiz_side_Bz";
  TH2D* hphizsideBz = new TH2D(os.str().c_str()," Bz (gauss); #phi (rad); z (m) ",360,0.,2.0*pi, 100, -hyperk_H, hyperk_H);
  os.str(""); os.clear(); os << tag << "_hphiz_side_Bperp";
  TH2D* hphizsideBperp = new TH2D(os.str().c_str()," Bperp (gauss); #phi (rad); z (m) ",360,0.,2.0*pi, 100, -hyperk_H, hyperk_H);

  for ( int iphi=1; iphi<=360; ++iphi ){
    for ( int iz=1; iz<=100; ++iz ){
      double phi = hphizsideB->GetXaxis()->GetBinCenter( iphi );
      double x = hyperk_R*std::cos( phi );
      double y = hyperk_R*std::sin( phi );      
      double z = hphizsideB->GetYaxis()->GetBinCenter( iz );
      TVector3 loc = TVector3{x, y, z} ;
      TVector3 btot = horzbs.get_B( loc ) + vertbs.get_B( loc ) + endbs.get_B( loc ) + BEarth;
      btot *= 1.0e4;//convert to gauss
      hphizsideB->SetBinContent( iphi, iz, btot.Mag() );
      hphizsideBx->SetBinContent( iphi, iz, btot.X() );
      hphizsideBy->SetBinContent( iphi, iz, btot.Y() );
      hphizsideBz->SetBinContent( iphi, iz, btot.Z() );

      TVector3 side = loc.Unit();
      TVector3 bperp = btot.Cross( side );
      hphizsideBperp->SetBinContent( iphi, iz, bperp.Mag() );
    }
  }
  hphizsideB->SetMaximum(0.5);
  hphizsideBx->SetMinimum(-0.5);
  hphizsideBx->SetMaximum(0.5);
  hphizsideBy->SetMinimum(-0.5);
  hphizsideBy->SetMaximum(0.5);
  hphizsideBz->SetMinimum(-0.5);
  hphizsideBz->SetMaximum(0.5);
  hphizsideBz->SetMinimum(0.0);
  hphizsideBz->SetMaximum(0.5);

  TCanvas * c=new TCanvas(tag.c_str(),tag.c_str(), 1000, 600 );
  c->Divide( 3, 2 );
  c->cd(1); hxyBx->Draw("colz");
  c->cd(2); hxyBy->Draw("colz");
  c->cd(3); hxyBz->Draw("colz");
  c->cd(4); hyzBx->Draw("colz");
  c->cd(5); hyzBy->Draw("colz");
  c->cd(6); hyzBz->Draw("colz");
  c->Write();
}


/// scan_chi2
/// Function to scan one or two of the parameters
/// Returns the minimum chi2 point via the parameters passed
/// minimum chi2 value returned
double scan_chi2( const int ipar1, double &parval1, const int nbins1, const double parmin1, const double parmax1, const std::string& parname1,
		  const int ipar2, double &parval2, const int nbins2, const double parmin2, const double parmax2, const std::string& parname2,
		  std::vector<double>& pars){

  std::cout<<"scan_chi2 starting"<<endl;
  std::cout<<"    par1="<<ipar1<<" with nbins="<<nbins1<<" from "<<parmin1<<" to "<<parmax1<<std::endl;
  std::cout<<"    par1="<<ipar2<<" with nbins="<<nbins2<<" from "<<parmin2<<" to "<<parmax2<<std::endl;
  double minchi2=9.e99;

  std::ostringstream os;
  os << "#Chi^{2}; "<< parname1 << "; "<< parname2;
  TH2D* hscanchi2 = new TH2D("hscanchi2",os.str().c_str(), nbins1, parmin1, parmax1, nbins2, parmin2, parmax2 ); 

  double step_p1 = (parmax1-parmin1)/double(nbins1);
  double step_p2 = (parmax2-parmin2)/double(nbins2);
  int ncoilset = 0;
  if (ipar1 <= npar-3 ) ++ncoilset;
  if (ipar2 <= npar-3 ) ++ncoilset;

  for ( double p1 = parmin1 + 0.5*step_p1; p1<parmax1; p1+=step_p1 ){
    for ( double p2 = parmin2 + 0.5*step_p2; p2<parmax2; p2+=step_p2 ){
      pars[ipar1] = p1;
      pars[ipar2] = p2;
      double chi2 = BFieldChi2( &pars[0] );
      if ( chi2 < minchi2 ){
	minchi2 = chi2;
	parval1 = p1;
	parval2 = p2;
      }
      hscanchi2->Fill( p1, p2, chi2 );
    }
  }
  
  pars[ipar1] = parval1;
  pars[ipar2] = parval2;
  make_histograms( pars, "minchi2scan" );

  return minchi2;
}




std::vector<double> get_nominal_coilspacing(){
  std::vector< double > pars;
  pars.push_back( 66.482051 );//66.582051;<-from MCMC );//  66.5 );   //56.5 ); // Ih (saddle)
  pars.push_back( 30.082835 );//  30.3 );  //29.6 ); // 44.9 ); //57.8 );//44.9 ); //Iv (solenoid)
  pars.push_back( 28.251163 );//   28.6 ); //28.6 ) ; //44.9 ); //53.3 );//44.9 ); //Iend (end)
  for (unsigned i=0; i<VertBodyZ.size(); ++i){
    pars.push_back( VertBodyZ[i] );
  }
  for (unsigned i=0; i<VertEndR.size(); ++i){
    pars.push_back( VertEndR[i] );
  }
  return pars;
}

std::vector<double> get_kylaverticalonly_coilspacing(){
  std::vector< double > pars;
  pars.push_back( 66.5 );//  66.5 );   //56.5 ); // Ih (saddle)
  pars.push_back( 23.8 );//  30.3 );  //29.6 ); // 44.9 ); //57.8 );//44.9 ); //Iv (solenoid)
  pars.push_back( 23.8 );//   28.6 ); //28.6 ) ; //44.9 ); //53.3 );//44.9 ); //Iend (end)
  for (unsigned i=0; i<VertBodyZ.size(); ++i){
    pars.push_back( VertBodyZ[i] );
  }
  for (unsigned i=0; i<VertEndR.size(); ++i){
    pars.push_back( VertEndR[i] );
  }
  return pars;
}


std::vector<double> get_kylaverticalonlysmall_coilspacing(){
  std::vector< double > pars;
  pars.push_back( 0.0 );//  66.5 );   //56.5 ); // Ih (saddle)
  pars.push_back( 44.9 );//  30.3 );  //29.6 ); // 44.9 ); //57.8 );//44.9 ); //Iv (solenoid)
  pars.push_back( 44.9 );//   28.6 ); //28.6 ) ; //44.9 ); //53.3 );//44.9 ); //Iend (end)
  for (unsigned i=0; i<VertBodyZ.size(); ++i){
    pars.push_back( VertBodyZ[i] );
  }
  for (unsigned i=0; i<VertEndR.size(); ++i){
    pars.push_back( VertEndR[i] );
  }
  return pars;
}




std::vector<double> get_hkverticalonly_coilspacing(){
  std::vector< double > pars;
  pars.push_back( 0.0 );//  66.5 );   //56.5 ); // Ih (saddle)
  pars.push_back( 67.0 );//73.82 );//  30.3 );  //29.6 ); // 44.9 ); //57.8 );//44.9 ); //Iv (solenoid)
  pars.push_back( 0.0 );//   28.6 ); //28.6 ) ; //44.9 ); //53.3 );//44.9 ); //Iend (end)
  for (unsigned i=0; i<VertBodyZ.size(); ++i){
    pars.push_back( VertBodyZ[i] );
  }

  for (unsigned i=0; i<VertEndR.size(); ++i){
    pars.push_back( VertEndR[i] );
  }
  return pars;
}



void do_mcmc_fit() {

  std::vector< double > sigmas;
  std::vector< double > pars = get_nominal_coilspacing();
  std::vector< std::string > parnames;
  parnames.push_back( "Ih" );
  parnames.push_back( "Iv" );
  parnames.push_back( "Iend" );
  sigmas.push_back( 0.2 );
  sigmas.push_back( 0.2 );
  sigmas.push_back( 0.2 );
  for (unsigned i=0; i<VertBodyZ.size(); ++i ){
    std::ostringstream os;
    os << "SideZ_" << i;
    parnames.push_back( os.str() );
    sigmas.push_back( 0.2 );
  }
  for (unsigned i=0; i<VertEndR.size(); ++i ){
    std::ostringstream os;
    os << "EndR_" << i;
    parnames.push_back( os.str() );
    sigmas.push_back( 0.2 );
  }

  std::vector< char * > parcharnames;
  for ( unsigned i=0; i< parnames.size(); ++i ){
    parcharnames.push_back( const_cast<char*>(parnames[i].c_str()) );
  }

  MCMC * mchain = new MCMC( ".", 1, pars.size(), &parcharnames[0], &pars[0], &sigmas[0] );

  mchain->RunMCMC( 10000, 8000, 1 );

  


}


void do_chi2_fit(){
  // Setup minimizer
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");

  min->SetMaxFunctionCalls(100000);
  min->SetMaxIterations(10000);
  min->SetTolerance(0.001);

  ROOT::Math::Functor fcn(&BFieldChi2,npar);
  std::vector< double > step;
  std::vector< double > pars = get_nominal_coilspacing();
  std::vector< std::string > parnames;
  parnames.push_back( "Ih" );
  parnames.push_back( "Iv" );
  parnames.push_back( "Iend" );
  step.push_back( 0.5 );
  step.push_back( 0.5 );
  step.push_back( 0.5 );
  for (unsigned i=0; i<VertBodyZ.size(); ++i ){
    std::ostringstream os;
    os << "SideZ_" << i;
    parnames.push_back( os.str() );
    step.push_back( 0.5 );
  }
  for (unsigned i=0; i<VertEndR.size(); ++i ){
    std::ostringstream os;
    os << "EndR_" << i;
    parnames.push_back( os.str() );
    step.push_back( 0.5 );
  }
  min->SetFunction( fcn );

  for (int i=0; i<npar; ++i ){
    min->SetVariable( i, parnames[i], pars[i], step[i] );
  }

  min->Minimize();

  min->PrintResults();
  
  // make plots for the minimum
  make_histograms( pars, "min" );
}

void check_geometry_and_field( std::vector< double > pars ){
  make_histograms( pars, "init" );
}

void biotsavart(){
  t2kstyle();
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  
  // Output file
  TFile * fout = new TFile("hyperkbs2.root","recreate");

  //do_mcmc_fit();
  //do_chi2_fit();
  double p1=67., p2 =2.;
  std::vector< double > pars = get_nominal_coilspacing();//get_hkverticalonly_coilspacing();

  //std::vector< double > pars = get_hkverticalonly_coilspacing();

  //std::vector< double > pars = get_kylaverticalonlysmall_coilspacing();

  double chi2 = BFieldChi2( &pars[0] );
  std::cout<<"Chi2 = "<<chi2<<std::endl;
  
  
  
  //double minchi2 = scan_chi2( 1, p1, 10, 63.2, 64.0, "Iv",
  // 			      3, p2, 1, 1.9, 2.1, "z1", pars ); 
  //pars[1] = p1;
  //pars[2] = p2;
  //std::cout<<"MinChi2="<<minchi2<<" for Iv="<<p1<<std::endl;
  check_geometry_and_field( pars );
  
  //double minchi2 = scan_chi2( 1, p1, 9, 52.0, 60.0, "Ih",
  // 			      2, p2, 1, 53.2, 53.4, "Iend", pars );
  //pars[1] = p1;
  //pars[2] = p2;
 

  
  //double minchi2 = scan_chi2( 1, p1, 1, 29.5, 29.7, "Iv",
  //			      2, p2, 10, 20.0, 40.0, "Iend", pars );
//pars[1] = p1;
// pars[2] = p2;
//std::cout<<"MinChi2="<<minchi2<<" for Iv="<<p1<<" Ih="<<p2<<std::endl;
  //check_geometry_and_field( pars );
  
  // plot for "Nominal"
  //std::vector< double > pars = get_nominal_coilspacing();
  //  double minchi2 = scan_chi2( npar-2, p1, 1, 59.0, 60.0, "Ih",
  // 			      npar-1, p2, 1, 75.0, 76.0, "Iv", pars );


  //std::vector< double > pars = get_equal_spaced_coils(); //Ih=56.5     Iv=91.2   Iend=36.5
  //double minchi2 = scan_chi2( npar-2, p1, 5, 80.0, 96.0, "Iv",
  //			      npar-1, p2, 1, 36.4, 36.6, "Iend", pars );

  
  // plot for "Current best"
  //  double minchi2 = scan_chi2( 14, p1, 5, 83.0, 84.0, "Ih",
  //			      15, p2, 17, 80.0, 96.0, "Iv" );
  
  fout->Write();
  fout->Close();
}


int main(){

  biotsavart();
  

  return 0;
}
  

