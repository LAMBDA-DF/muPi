
#ifndef MYLISTBOX
#define MYLISTBOX


#include <string>
#include <vector>
#include <TApplication.h>
#include <TGClient.h>
#include <TGButton.h>
#include <TGListBox.h>
#include <TList.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TRootEmbeddedCanvas.h>

//#include "plotEvtWin.h"

class TTree;
class TH2F;
class TGTextEntry;
class TH1F;
class TArrow;
class TPaveText;
class TGTextView;
class TGFrame;

struct listEntry_t
{
    std::string name;
    int index;
    listEntry_t(std::string name,int index):name(name),index(index){};
};

class listBox : public TGMainFrame {
    
  public:
    listBox(const TGWindow *p, UInt_t w, UInt_t h, TTree *hitSumm);
    virtual ~listBox();
    void DoExit();
    void HandleButtons();
    void PrintSelected();
    void DoubleClicked(TGFrame *f, Int_t btn);
    TCanvas *GetEvtCanvas(){ return fCanvas->GetCanvas(); };
    TCanvas *GetG1Canvas(){ return fG1->GetCanvas(); };
    TCanvas *GetG2Canvas(){ return fG2->GetCanvas(); };
  
    void buildEventList(const char *selectionCuts);
    void changePageEventList(int i);
    void drawHistogram();
    void drawItem(int i);
    void drawEventArrow(int i);
    void displayCurrValues(int i);
    void displayXY(int i);
    void nothingToPlot(const char *line1="No events",const char *line2="nothing to plot");
    void showHelpWindow();
    void selectionChanged(TGFrame *f);
    void plotHit(const Int_t n, const Int_t* x, const Int_t* y, const Float_t* ePix);
    void applyCuts();
    void exec2event(Int_t event, Int_t x, Int_t y, TObject *selected);
    
  private:
    TRootEmbeddedCanvas *fCanvas;
    TRootEmbeddedCanvas *fG1;
    TRootEmbeddedCanvas *fG2;
    TRootEmbeddedCanvas *fG3;
    TGListBox           *fListBox;
    TGCheckButton       *fCheckMulti;
    TGCheckButton       *fCheckLogHisto;
    TList               *fSelected;   
    TCanvas             *fCan;
    TTree               *fHitSumm;
    TTree               *fSelHitSumm;
    TH2F                *fHistTrack;
    TH1F                *fHistSel;
    TH1F                *fHistAll;
    TGTextEntry         *fSelCuts;
    TGTextEntry         *fHistFormula;
    TFile               *fCacheFile;
    TArrow              *fEventArrow;
    
    Float_t             fHistSelMinY;
    Float_t             fHistSelMaxY;
    
    TPaveText           *fNothingToPlot;
    TPaveText           *fNoFormula;
    
    TGFrame             *fHelpWindow;

    std::vector<listEntry_t> *fSelHitList;
    int                 fLastEntryInListBox;
    int                 fFirstEntryInListBox;

    ClassDef(listBox, 0)
};

#endif
