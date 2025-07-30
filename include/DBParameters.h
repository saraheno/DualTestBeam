#pragma once

class DBParameters {
 public:
  // way too slow if track all photons for now
  // so randomly delete photons after creation according to this fraction
  //   dialScint=1.0, dialCer=1.0 to keep all photons
  double m_dialCherC=10./800000.;
  double m_dialScintC=100./20000000.;
  double m_dialCherO= 100./800000.;
  double m_dialScintO=1./20000000.;
  float m_betarel=1/1.544;
  int m_printlimitSCE=10;
  int m_MAXEVENTSCE=10;

  static DBParameters* instance() {
    static DBParameters* theBDP = new DBParameters();
    return theBDP;
  }
  
 private:
  DBParameters() {}
};
