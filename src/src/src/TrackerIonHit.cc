
#include "TrackerIonHit.hh"

G4Allocator<TrackerIonHit> TrackerIonHitAllocator;


TrackerIonHit::TrackerIonHit() { flag=0;}



TrackerIonHit::~TrackerIonHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerIonHit::TrackerIonHit(const TrackerIonHit& right)
  : G4VHit()
{
  trackID   = right.trackID;
  particleID= right.particleID;
  pos       = right.pos;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const TrackerIonHit& TrackerIonHit::operator=(const TrackerIonHit& right)
{
  trackID   = right.trackID;
  particleID= right.particleID;
  pos       = right.pos;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int TrackerIonHit::operator==(const TrackerIonHit& right) const
{
  return (this==&right) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerIonHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(pos);
    circle.SetScreenSize(4);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,1.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerIonHit::Print()
{
  std::ios init(NULL);
  init.copyfmt(G4cout);

  G4cout << std::setw(2)<< std::right << trackID <<" "<<flag<<" "
	 << std::setw(15)<< particleID<<" "
	 << std::setprecision(4)<<std::setw(8)<<std::fixed
	 << KE/GeV<<" "<<std::setprecision(6)<<std::setw(8)
	 << beta<< " "
	 <<std::setprecision(4)<<std::setw(10)<<std::right
	 <<theta*1000.<<" "<<std::setw(12)<<std::right
	 <<phi*1000.<<" "<<std::setw(10)<<std::right	 
	 <<pos.getX()<<" "<<std::setw(10)<<std::right
	 <<pos.getY()<<" "<<std::setw(10)<<std::right
	 <<pos.getZ()<<" "<<std::setw(12)<<std::right
	 <<time*1000.<<" "<<std::setw(12)<<std::right
	 <<labtime*1000.<<" "<<std::setw(12)<<std::right
	 <<globaltime*1000.<<" "<<std::setw(10)<<std::right
	 <<Edep/MeV<<" "<<std::setw(10)<<std::right
	 <<length/mm
	 << G4endl;

  G4cout.copyfmt( init );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

