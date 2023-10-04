#include <sstream>
#include <string>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <functional>
#include <cmath>
#include <iostream>
using namespace std;

#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TKey.h"
#include "TDirectory.h"

#include "TLine.h"
#include "TError.h"
#include "TNtuple.h"
#include "TVector.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TGaxis.h"
#include "TApplication.h"
#include <TGClient.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TChain.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TROOT.h>
#include "TStyle.h"
#include <TRootEmbeddedCanvas.h>
#include <RQ_OBJECT.h>


#include "listBox.h"

void printHelp(const char* exeName){
  cout << "\nThis program plots the hits stored in the root file produced by extract.exe\n";
  cout << "\nUse:\n  " << exeName <<" <input root file>\n\n";
  return;
}

bool fileExist(const char *fileName){
  ifstream in(fileName,ios::in);
  if(in.fail()){
    in.close();
    return false;
  }
  in.close();
  return true;
}


int main(int argc,char *argv[]){
  
  TApplication myapp("myapp", 0, 0);
  
  if(argc!=2){
    printHelp(argv[0]);
    return 1;
  }
  
  if( !fileExist(argv[1]) ){
    cerr << "\nError: file doesn't exist!\n\n";
    return 1;
  }
  
  TFile inF(argv[1]);
  TTree *hitSumm = (TTree*)inF.Get("hitSumm");
//   hitSumm->SetDirectory(0);
  
//   TChain *hitSumm = new TChain("hitSumm");
//   hitSumm->Add(argv[1]);
  
  gErrorIgnoreLevel = 1001;
  gStyle->SetOptStat(0);
  
  listBox lb(gClient->GetRoot(), 100,100, hitSumm);
  lb.Connect("CloseWindow()","TApplication",&myapp,"Terminate()");
  myapp.Run();

  return 0;
}
