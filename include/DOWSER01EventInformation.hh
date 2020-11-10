#ifndef DOWSER01EventInformation_hh
#define DOWSER01EventInformation_hh 1

#include "globals.hh"
#include "G4VUserEventInformation.hh"

class DOWSER01EventInformation : public G4VUserEventInformation
{
public:
  DOWSER01EventInformation()
  :theta0(0.),theta1(0.), energyN(0.),positionX0(0.), positionY0(0.),positionX1(0.), positionY1(0.)  {}
  ~DOWSER01EventInformation() {}
  
  // Mandatory function
  void Print() const {
    G4cout
    <<"Theta0: "<<'\t'<<theta0<<'\n'
    <<"Theta1: "<<'\t'<<theta1<<'\n'
    <<"EnergyN: "<<'\t'<<energyN<<'\n'
    <<"PositionX0: "<<'\t'<<positionX0<<'\n'
    <<"PositionY0: "<<'\t'<<positionY0<<'\n'
    <<"PositionX1: "<<'\t'<<positionX1<<'\n'
    <<"PositionY1: "<<'\t'<<positionY1<<'\n'
    <<G4endl;
  }
    
  
private:
  G4double theta0,theta1, energyN, positionX0, positionY0, positionX1, positionY1;
  
public:
  // Set
  inline void SetTheta0(G4double t) {theta0 = t;}
  inline void SetTheta1(G4double t) {theta1 = t;}
  inline void SetEnergyN(G4double t) {energyN = t;}
  inline void SetPositionX0(G4double t) {positionX0 = t;}
  inline void SetPositionY0(G4double t) {positionY0 = t;}
    inline void SetPositionX1(G4double t) {positionX1 = t;}
    inline void SetPositionY1(G4double t) {positionY1 = t;}
  
  // Get
  inline G4double GetTheta0() const {return theta0;}
  inline G4double GetTheta1() const {return theta1;}
  inline G4double GetEnergyN() const {return energyN;}
  inline G4double GetPositionX0() const {return positionX0;}
  inline G4double GetPositionY0() const {return positionY0;}
    inline G4double GetPositionX1() const {return positionX1;}
    inline G4double GetPositionY1() const {return positionY1;}

};

#endif
