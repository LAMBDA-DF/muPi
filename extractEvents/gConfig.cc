#include "tinyxml2.h"
#include "gConfig.h"
#include "globalConstants.h"
#include <iostream>
#include <sstream>


using namespace std;

gConfig::gConfig(): fStackSize(128), fDefaultSigma(-1), fDefaultCal(-1), fSeedThr(-1), fAddThr(-1), fAddRad(1), fSkirtSize(0)
{}

float gConfig::getExtSigma(const int ext){
  
  std::map<int, float>::const_iterator it = fExt2Sigma.find( ext );
  if ( it == fExt2Sigma.end() ) {
    return fDefaultSigma;
  }
  return it->second;
}

float gConfig::getExtCal(const int ext){
  
  std::map<int, float>::const_iterator it = fExt2Cal.find( ext );
  if ( it == fExt2Cal.end() ) {
    return fDefaultCal;
  }
  return it->second;
}

void gConfig::printVariables(){
  cout << endl;
  cout << "===== Configuration parameters =====\n";
  cout << "sigma:\n";
  cout << "\tdefault: " << fDefaultSigma << endl;
  typedef std::map<int, float>::iterator it_type;
  for(it_type it = fExt2Sigma.begin(); it != fExt2Sigma.end(); it++) {
    cout << "\text" << it->first << ": " << it->second << endl;
  }
  
  cout << "thr:\n";
  cout << "\tseed:       " << fSeedThr << endl;
  cout << "\tadd:        " << fAddThr << endl;
  cout << "\tadd radius: " << fAddRad << endl;
  cout << "\tskirtSize:  " << fSkirtSize << endl;
  
  cout << "cal:\n";
  cout << "\tdefault: " << fDefaultCal << endl;
  for(it_type it = fExt2Cal.begin(); it != fExt2Cal.end(); it++) {
    cout << "\text" << it->first << ": " << it->second << endl;
  }
  
  cout << "extra:\n";
  cout << "\tsaveTracks:    " << fSaveTracks << endl;
  cout << "\tsaveTrackCuts: " << fTracksCuts << endl;  
  
  cout << "system:\n";
  cout << "\tstackSize:    " << fStackSize << endl;
  cout << "\thitMaxSize:   " << fHitMaxSize << endl;
  cout << "====================================\n";
}
  
bool gConfig::readConfFile(const char* confFileName = "extractConfig.xml"){
  
  tinyxml2::XMLDocument doc;
  if(doc.LoadFile( confFileName ) != 0){
    cerr << red << "\nCan't read config file! Will not continue.\n\n" << normal;
    return false;
  }
//   doc.Print();
  
  /* Sigmas */
  if( doc.FirstChildElement("sigma") == 0 ){
    cerr << "Missing \'sigma\' in config file.\n\n";
    return false;
  }
  
  if(doc.FirstChildElement("sigma")->QueryFloatAttribute("default", &fDefaultSigma) !=0){
    cerr << "Missing \'sigma:default\' in config file.\n\n";
    return false;
  }
  
  for(int i=0;i<100;++i){
    ostringstream ossExtSig;
    ossExtSig << "ext" << i;
    float extSig = -1;
    if(doc.FirstChildElement("sigma")->QueryFloatAttribute(ossExtSig.str().c_str(), &extSig) == 0){
      fExt2Sigma[i] = extSig;
    }
  }

  /* Calibration */
  if( doc.FirstChildElement("calibration") == 0 ){
    cerr << "Missing \'calibration\' in config file.\n\n";
    return false;
  }
  
  if(doc.FirstChildElement("calibration")->QueryFloatAttribute("default", &fDefaultCal) !=0){
    cerr << "Missing \'calibration:default\' in config file.\n\n";
    return false;
  }
  
  for(int i=0;i<100;++i){
    ostringstream ossExtSig;
    ossExtSig << "ext" << i;
    float extCal = -1;
    if(doc.FirstChildElement("calibration")->QueryFloatAttribute(ossExtSig.str().c_str(), &extCal) == 0){
      fExt2Cal[i] = extCal;
    }
  }
  
  /* Thresholds */
  if( doc.FirstChildElement("thr") == 0 ){
    cerr << "Missing \'thr\' in config file.\n\n";
    return false;
  }
  
  if(doc.FirstChildElement("thr")->QueryFloatAttribute("seed", &fSeedThr) != 0){
    cerr << "Missing \'thr:fSeedThr\' in config file.\n\n";
    return false;
  }
  
  if(doc.FirstChildElement("thr")->QueryFloatAttribute("add", &fAddThr) != 0){
    cerr << "Missing \'thr:fAddThr\' in config file.\n\n";
    return false;
  }

  if(doc.FirstChildElement("thr")->QueryIntAttribute("addRad", &fAddRad) != 0){
    cerr << "Missing \'thr:fAddRad\' in config file.\n\n";
    return false;
  }

  if(doc.FirstChildElement("thr")->QueryIntAttribute("skirtSize", &fSkirtSize) != 0){
    cerr << "Missing \'thr:fSkirtSize\' in config file.\n\n";
    return false;
  }
  
  /* Extra */
  if( doc.FirstChildElement("extra") == 0 ){
    cerr << "Missing \'extra\' in config file.\n\n";
    return false;
  }
  
  if(doc.FirstChildElement("extra")->QueryBoolAttribute("saveTracks", &fSaveTracks) != 0){
    fSaveTracks = false;
  }
  
  fTracksCuts = doc.FirstChildElement("extra")->Attribute("saveTrackCuts");
  
  /* System */
  const int kDefaultStackSize = 128;
  if( doc.FirstChildElement("systemConfig") == 0 ){
    cerr << "Missing \'systemConfig\' in config file.\n\n";
    return 1;
  }
  
  if( doc.FirstChildElement("systemConfig")->QueryIntAttribute("stackSize", &fStackSize) != 0 ){
    fStackSize = kDefaultStackSize*8;
  }
  
  if( doc.FirstChildElement("systemConfig")->QueryIntAttribute("hitMaxSize", &fHitMaxSize) != 0){
    fHitMaxSize = kDefaultStackSize*20;
  }
  
  /* Variables that will be saved in the NTuple */
  for(int n=0;n<gNBaseTNtupleVars;++n){
    if(n>0) fNTupleVars += ":";
    fNTupleVars += gBaseTNtupleVars[n];
  }
  for(int l=0;l<=fSkirtSize;++l){
    ostringstream newVars;
    for(int n=0;n<gNExtraTNtupleVars;++n){
      newVars << ":" << gExtraTNtupleVars[n] << l;
    }
    fNTupleVars += newVars.str();
  }
  
  return true;
}
