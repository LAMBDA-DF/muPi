#include <TApplication.h>
#include <TGClient.h>
#include <TGButton.h>
#include <TGListBox.h>
#include <TGTextEntry.h>
#include <TList.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <TH2F.h>
#include <TTree.h>
#include <TFile.h>
#include <TLine.h>
#include <TLatex.h>
#include <TArrow.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TGTextView.h>
#include <TGFrame.h>
#include <TGCanvas.h>
#include <TGFrame.h>
#include <TGTextEdit.h>
#include "TGScrollBar.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>

#include "listBox.h"

using namespace std;

const int kMaxNPixels = 50000;
const int kMinRadPix  = 10;
const int kMaxEntriesInBox = 2000;
const int kVerbosity = 3;

void listBox::drawHistogram(){
  
  if(fHistSel!=0){
    delete fHistSel;
    fHistSel=0;
  }
  
  string auxFormulaString(fHistFormula->GetText());
  if(auxFormulaString.find_first_not_of(" /n/r/t")==string::npos){
    fG2->GetCanvas()->Clear();
    fG2->GetCanvas()->Update();
    return;
  }
  
  fSelHitSumm->SetBranchStatus("*",1);
  fG2->GetCanvas()->cd();
  fSelHitSumm->Draw(fHistFormula->GetText(),"","goff");
  fHistSel = (TH1F*) (fSelHitSumm->GetHistogram()->Clone("fHistSel"));
  fSelHitSumm->SetBranchStatus("*",0);
  
  fHistSel->SetLineColor(kBlue);
  fHistSel->Draw();
  
  fHitSumm->SetBranchStatus("*",1);
  fHitSumm->SetLineColor(kBlack);
  fHitSumm->Draw(fHistFormula->GetText(),"","same");
  fHitSumm->SetBranchStatus("*",0);
  
  fHistSelMinY = fHistSel->GetMinimum()+1;
  fHistSelMaxY = fHistSel->GetMaximum()*(1+gStyle->GetHistTopMargin());

  drawEventArrow(fListBox->GetSelected());
  fG2->GetCanvas()->Update();  
}

void listBox::buildEventList(const char *selectionCuts=""){
  
  fCacheFile->cd();
    
  if(fSelHitSumm!=0){
    delete fSelHitSumm;
    fSelHitSumm=0;
  }
  fListBox->RemoveAll();
  fSelHitList->clear();
  
  fHitSumm->SetBranchStatus("*",1);
  
  fSelHitSumm = fHitSumm->CopyTree(selectionCuts);
  fSelHitSumm->SetName("fSelHitSumm");
  
  fSelHitSumm->SetBranchStatus("*",0);
  if (kVerbosity>5) cout << "buildEventList\n";
  Int_t n;
  fSelHitSumm->SetBranchStatus("nSavedPix", 1);
  fSelHitSumm->SetBranchAddress("nSavedPix", &n);
  
  long nEntries=fSelHitSumm->GetEntries();
  long firstEntry = -1;
  long listEntryNumber = 1;
  for(long i = 0; i<nEntries; ++i){
    fSelHitSumm->GetEntry(i);
    if(n>0){
      if(firstEntry<0) firstEntry=i;
      ostringstream name;
      name << listEntryNumber;
      if (listEntryNumber<kMaxEntriesInBox)
      {
        fListBox->AddEntry(name.str().c_str(), i);
      }
      fSelHitList->push_back(listEntry_t(name.str(), i));
      ++listEntryNumber;
    }
  }
  fLastEntryInListBox = fListBox->GetNumberOfEntries();
  if(listEntryNumber>=kMaxEntriesInBox) fListBox->AddEntry("next page", -100);
  fListBox->Layout();
  
  if (kVerbosity>0) cout << fListBox->GetNumberOfEntries() << " events in the list\n"; 
 
  if(nEntries>0){
    fListBox->Select(firstEntry);
    drawItem(firstEntry);
    drawHistogram();
  }
  else{
    nothingToPlot();
  }
  
}


void listBox::plotHit(const Int_t n, const Int_t* x, const Int_t* y, const Float_t* ePix){
  
  if(fHistTrack!=0){
    delete fHistTrack;
  }
    
  const Int_t xMin = *(min_element( x, &(x[n-1])));
  const Int_t xMax = *(max_element( x, &(x[n-1])));
  const Int_t xRad = (xMax-xMin)/2 +1;
  const Int_t xCen = xRad+xMin;
  
  const Int_t yMin = *(min_element( y, &(y[n-1])));
  const Int_t yMax = *(max_element( y, &(y[n-1])));
  const Int_t yRad = (yMax-yMin)/2 +1;
  const Int_t yCen = yRad+yMin;
  
  const Int_t nRad  = max(max(xRad+1,yRad+1), kMinRadPix);
  const Int_t nBins = nRad*2;
  
  fHistTrack = new TH2F("hPix","",nBins,xCen-nRad,xCen+nRad, nBins,yCen-nRad,yCen+nRad);
  fHistTrack->SetStats(0);
  for(Int_t i=0; i<n; ++i){
    fHistTrack->Fill(x[i],y[i],ePix[i]);
  }
  
  fHistTrack->Draw("colz");
}

void listBox::drawEventArrow(int i){
  
  if (kVerbosity>5) cout << "listBox::drawEventArrow 1\n";
  string auxFormulaString(fHistFormula->GetText());
  if(auxFormulaString.find_first_not_of(" /n/r/t")==string::npos){
    cout << "No formula, no arrow!\n";
    return;
  }
  if (kVerbosity>5) cout << "listBox::drawEventArrow 2\n";

  fG2->GetCanvas()->cd();
  fSelHitSumm->ResetBranchAddresses();
  fSelHitSumm->SetBranchStatus("*",0);  
  fSelHitSumm->SetBranchStatus("*",1);
  ostringstream ossHistCuts;
  ossHistCuts << "Entry$ == " << i;
  if (kVerbosity>5) cout << "listBox::drawEventArrow 3\n";
  fSelHitSumm->Draw(fHistFormula->GetText(), ossHistCuts.str().c_str(), "goff");
  if (kVerbosity>5) cout << "listBox::drawEventArrow 4\n";
  Float_t xEvt = fSelHitSumm->GetHistogram()->GetMean();
  if (kVerbosity>5) cout << "listBox::drawEventArrow 5\n";
  if(fEventArrow!=0){
    delete fEventArrow;
    fEventArrow=0;
  }
  if (kVerbosity>5) cout << "listBox::drawEventArrow 6\n";
  fEventArrow = new TArrow(xEvt,fHistSelMaxY,xEvt,fHistSelMinY,0.05,"|>");
  fEventArrow->SetLineColor(kRed);
  fEventArrow->SetFillColor(kRed);
  fEventArrow->Draw();
  if (kVerbosity>5) cout << "listBox::drawEventArrow 7\n";
  fSelHitSumm->SetBranchStatus("*",0);
  if (kVerbosity>5) cout << "listBox::drawEventArrow 8\n";
}

void listBox::drawItem(int i){
  if (i<0)
  {
    changePageEventList(i);
    return;
  }

  if (kVerbosity>1) cout << "listBox::drawItem(" << i << ")\n";
  
  fSelHitSumm->SetBranchStatus("*",0);
  if (kVerbosity>5) cout << "listBox::drawItem 2\n";
  Int_t n;
  fSelHitSumm->SetBranchStatus("nSavedPix", 1);
  fSelHitSumm->SetBranchAddress("nSavedPix", &n);
  fSelHitSumm->GetEntry(i);
  if(n > kMaxNPixels){
    cout << "\nERROR: too many pixels!\n\n";
    nothingToPlot("ERROR","too many pixels");
    return;
  }
  
  Int_t x[kMaxNPixels];
  fSelHitSumm->SetBranchStatus("xPix", 1);
  fSelHitSumm->SetBranchAddress("xPix", &x);
  
  Int_t y[kMaxNPixels];
  fSelHitSumm->SetBranchStatus("yPix", 1);
  fSelHitSumm->SetBranchAddress("yPix", &y);
  
  Float_t ePix[kMaxNPixels];
  fSelHitSumm->SetBranchStatus("ePix", 1);
  fSelHitSumm->SetBranchAddress("ePix", &ePix);
  if (kVerbosity>5) cout << "listBox::drawItem 3\n";
  fCanvas->GetCanvas()->cd();
  // fCanvas->GetCanvas()->SetLogz();
  fCanvas->GetCanvas()->SetRightMargin(0.12);

  fSelHitSumm->GetEntry(i);
  plotHit(n,x,y,ePix);
  fCanvas->GetCanvas()->Update();
  if (kVerbosity>5) cout << "listBox::drawItem 4\n";
  string auxFormulaString(fHistFormula->GetText());
  if(auxFormulaString.find_first_not_of(" /n/r/t")!=string::npos){
    fG2->GetCanvas()->cd();
    if (kVerbosity>5) cout << "listBox::drawItem 5\n";
    if(fCheckLogHisto->GetState()>0)
      fG2->GetCanvas()->SetLogy(1);
    else
      fG2->GetCanvas()->SetLogy(0);
    
    displayCurrValues(i);
    displayXY(i);
    if (kVerbosity>5) cout << "listBox::drawItem 6\n";
    drawEventArrow(i);
    if (kVerbosity>5)cout << "listBox::drawItem 7\n";
  }
  if (kVerbosity>5) cout << "listBox::drawItem 8\n";
  fG2->GetCanvas()->Update();
}

void listBox::DoExit()
{
   Printf("Slot DoExit()");
   gApplication->Terminate(0);
}

void listBox::changePageEventList(int code){
  fListBox->RemoveAll();
  if (kVerbosity>2) cout << "changePageEventList: " << code << endl;
  const unsigned int listSize = fSelHitList->size();
  int selectEntry = 0;

  if(code==-101){ //Prev page
    fFirstEntryInListBox = fFirstEntryInListBox-kMaxEntriesInBox > 0 ?       fFirstEntryInListBox-kMaxEntriesInBox     : 0;
    fLastEntryInListBox  = fLastEntryInListBox- kMaxEntriesInBox   < listSize? fLastEntryInListBox-kMaxEntriesInBox      : listSize;
    selectEntry          = fSelHitList->at(fLastEntryInListBox).index;
  }
  if(code==-100){ //Next page
    fFirstEntryInListBox = fFirstEntryInListBox+kMaxEntriesInBox > 0 ?       fFirstEntryInListBox+kMaxEntriesInBox   : 0;
    fLastEntryInListBox  = fLastEntryInListBox+ kMaxEntriesInBox   < listSize? fLastEntryInListBox+kMaxEntriesInBox      : listSize;
    selectEntry          = fSelHitList->at(fFirstEntryInListBox).index;
  }
  
  if(fFirstEntryInListBox>0) fListBox->AddEntry("prev page", -101);
  for (int i = fFirstEntryInListBox; i <= fLastEntryInListBox; ++i)
  {   
    fListBox->AddEntry(fSelHitList->at(i).name.c_str(), fSelHitList->at(i).index);
  }
  if(fLastEntryInListBox<listSize) fListBox->AddEntry("next page", -100);
  fListBox->Layout();
  drawItem(selectEntry);
  if (kVerbosity>2) cout << "changePageEventList: Select " << selectEntry << endl;
  fListBox->Layout();
  fListBox->Select(selectEntry);
  if (kVerbosity>2) cout << "changePageEventList: done Select " << endl;
  fListBox->Layout();
}

void listBox::nothingToPlot(const char* line1,const char* line2){
  

  fNothingToPlot->Clear();
  fNothingToPlot->SetTextSize(0.08);
  fNothingToPlot->AddText(" ");
  fNothingToPlot->AddText(line1);
  fNothingToPlot->AddText(line2);
  fNothingToPlot->AddText(" ");

  fG1->GetCanvas()->cd();
  fNothingToPlot->Draw();
  fG1->GetCanvas()->Update();
  
  fG2->GetCanvas()->cd();
  fNothingToPlot->Draw();
  fG2->GetCanvas()->Update();

  fG3->GetCanvas()->cd();
  fNothingToPlot->Draw();
  fG3->GetCanvas()->Update();
  
  fCanvas->GetCanvas()->cd();
  fNothingToPlot->Draw();
  fCanvas->GetCanvas()->Update();
}

void listBox::displayXY(int evN){

  float sizeX = 4200;
  float sizeY = 2100;

  if (kVerbosity>5) cout << "listBox::displayXY 1\n";
  
  fSelHitSumm->SetBranchStatus("*",0);
  if (kVerbosity>5) cout << "listBox::displayXY 2\n";
  const Int_t nIntVars = 1;
  const char* intVarsName[nIntVars] = {"ohdu"};
  Int_t intVars[nIntVars];
  for(int i=0;i<nIntVars;++i){
    fSelHitSumm->SetBranchStatus(intVarsName[i],1);
    fSelHitSumm->SetBranchAddress(intVarsName[i], &(intVars[i]));
  }
  const Int_t nFloatVars = 2;
  const char* floatVarsName[nFloatVars] = { "xBary0", "yBary0" };
  Float_t floatVars[nFloatVars];
  for(int i=0;i<nFloatVars;++i){
    fSelHitSumm->SetBranchStatus(floatVarsName[i],1);
    fSelHitSumm->SetBranchAddress(floatVarsName[i], &(floatVars[i]));
  }

  fG3->GetCanvas()->cd();
  fG3->GetCanvas()->Clear();

  std::vector<TMarker*> vMarkEvts(fSelHitList->size());
  for(int h=0;h<fSelHitList->size();++h){
    int evnt = fSelHitList->at(h).index;
    fSelHitSumm->GetEntry(evnt);
    vMarkEvts[h] = new TMarker(floatVars[0]/sizeX, floatVars[1]/sizeY, 6);
    vMarkEvts[h]->SetMarkerColor(40+intVars[0]);
    vMarkEvts[h]->Draw();
  }
  fSelHitSumm->GetEntry(evN);
  

  TMarker *evt = new TMarker(floatVars[0]/sizeX, floatVars[1]/sizeY, 24);
  evt->SetMarkerColor(kRed);

  evt->Draw();
  fG3->GetCanvas()->Update();
  fSelHitSumm->SetBranchStatus("*",0);

}


void listBox::displayCurrValues(int evN){
  if (kVerbosity>5) cout << "listBox::displayCurrValues 1\n";
  
  fSelHitSumm->SetBranchStatus("*",0);
  if (kVerbosity>5) cout << "listBox::displayCurrValues 2\n";
  const Int_t nIntVars = 4;
  const char* intVarsName[nIntVars] = { "runID", "ohdu", "nSat", "nSavedPix" };
  Int_t intVars[nIntVars];
  for(int i=0;i<nIntVars;++i){
    fSelHitSumm->SetBranchStatus(intVarsName[i],1);
    fSelHitSumm->SetBranchAddress(intVarsName[i], &(intVars[i]));
  }
  if (kVerbosity>5) cout << "listBox::displayCurrValues 3\n";
  const Int_t nFloatVars = 4;
  const char* floatVarsName[nFloatVars] = { "E0", "n0", "xVar0", "yVar0" };
  Float_t floatVars[nFloatVars];
  for(int i=0;i<nFloatVars;++i){
    fSelHitSumm->SetBranchStatus(floatVarsName[i],1);
    fSelHitSumm->SetBranchAddress(floatVarsName[i], &(floatVars[i]));
  }
  
  fSelHitSumm->GetEntry(evN);
  
  fG1->GetCanvas()->cd();
  fG1->GetCanvas()->Clear();

  TLatex Tl;		
  float col1  = 0.05;
  float col2  = 0.55;
  float yTex  = 0.1;
  float stepTex = min(1./(nIntVars+1),1./(nFloatVars+1));
  Tl.SetTextSize(stepTex*0.95);
  if (kVerbosity>5) cout << "listBox::displayCurrValues 4\n";
  for(int i=0;i<nIntVars;++i){
    ostringstream lineOSS;     
    lineOSS << " " << intVarsName[i] << ": " << std::setprecision(3) << intVars[i];
  	Tl.DrawLatex(col1,yTex+i*stepTex,lineOSS.str().c_str());
  }
  if (kVerbosity>5) cout << "listBox::displayCurrValues 5\n";
  for(int i=0;i<nFloatVars;++i){
    ostringstream lineOSS;
    lineOSS << " " << floatVarsName[i] << ": " << std::setprecision(3) << floatVars[i];
  	Tl.DrawLatex(col2,yTex+i*stepTex,lineOSS.str().c_str());
  }
  if (kVerbosity>5) cout << "listBox::displayCurrValues 6\n";
  fG1->GetCanvas()->Update();
  fSelHitSumm->SetBranchStatus("*",0);
}

void listBox::exec2event(Int_t event, Int_t x, Int_t y, TObject *selected)
{
  TCanvas *c = (TCanvas *) gTQSender;
  
  if (event == kButton1Down) // mouse click
    printf("Canvas %s: event=%d, x=%d, y=%d\n", c->GetName(), event, x, y);
}

void listBox::selectionChanged(TGFrame *f){
  TGTextLBEntry *lbe = (TGTextLBEntry *)f;
  int id = lbe->EntryId();
  drawItem(id);
}

void listBox::applyCuts(){
  if (kVerbosity>2) cout << fSelCuts->GetText() << "\n";
  
  buildEventList(fSelCuts->GetText());
}



//=================================================
class Editor {

private:
   TGTransientFrame *fMain;   // main frame of this widget
   TGTextView       *fEdit;   // text edit widget
   TGTextButton     *fOK;     // OK button
   TGLayoutHints    *fL1;     // layout of TGTextView
   TGLayoutHints    *fL2;     // layout of OK button

public:
   Editor(const TGWindow *main, UInt_t w, UInt_t h);
   virtual ~Editor();

   void   LoadFile(const char *file);
   void   LoadBuffer(const char *buffer);
   void   AddBuffer(const char *buffer);

   TGTextView *GetEditor() const { return fEdit; }

   void   SetTitle();
   void   Popup();

   // slots
   void   CloseWindow();
   void   DoOK();
   void   DoOpen();
   void   DoSave();
   void   DoClose();
};

Editor::Editor(const TGWindow *main, UInt_t w, UInt_t h)
{
   // Create an editor in a dialog.

   fMain = new TGTransientFrame(gClient->GetRoot(), main, w, h);
//    fMain->Connect("CloseWindow()", "Editor", this, "CloseWindow()");
//    fMain->DontCallClose(); // to avoid double deletions.

   // use hierarchical cleaning
   fMain->SetCleanup(kDeepCleanup);

   fEdit = new TGTextView(fMain, w, h, kSunkenFrame | kDoubleBorder);
   fL1 = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 3, 3, 3, 3);
   fMain->AddFrame(fEdit, fL1);
//    fEdit->Connect("Opened()", "Editor", this, "DoOpen()");
//    fEdit->Connect("Saved()",  "Editor", this, "DoSave()");
//    fEdit->Connect("Closed()", "Editor", this, "DoClose()");

   // set selected text colors
   Pixel_t pxl;
   gClient->GetColorByName("#3399ff", pxl);
   fEdit->SetSelectBack(pxl);
   fEdit->SetSelectFore(TGFrame::GetWhitePixel());
   fEdit->LoadFile("help.txt");

//    fOK = new TGTextButton(fMain, "  &OK  ");
//    fOK->Connect("Clicked()", "Editor", fMain, "CloseWindow()");
//    fL2 = new TGLayoutHints(kLHintsBottom | kLHintsCenterX, 0, 0, 5, 5);
//    fMain->AddFrame(fOK, fL2);

   SetTitle();

   fMain->MapSubwindows();

   fMain->Resize();

   // editor covers right half of parent window
   fMain->CenterOnParent(kTRUE, TGTransientFrame::kRight);
}

Editor::~Editor()
{
   // Delete editor dialog.

   fMain->DeleteWindow();  // deletes fMain
}

void Editor::SetTitle()
{
   // Set title in editor window.

   TGText *txt = GetEditor()->GetText();
   Bool_t untitled = !strlen(txt->GetFileName()) ? kTRUE : kFALSE;

   char title[256];
   if (untitled)
      sprintf(title, "Help window");
   else
      sprintf(title, "Help window - %s", txt->GetFileName());

   fMain->SetWindowName(title);
   fMain->SetIconName(title);
}

void Editor::Popup()
{
   // Show editor.

   fMain->MapWindow();
}

void Editor::LoadBuffer(const char *buffer)
{
   // Load a text buffer in the editor.

   fEdit->LoadBuffer(buffer);
}

void Editor::LoadFile(const char *file)
{
   // Load a file in the editor.

   fEdit->LoadFile(file);
}

void Editor::AddBuffer(const  char *buffer)
{
   // Add text to the editor.

   TGText txt;
   txt.LoadBuffer(buffer);
   fEdit->AddText(&txt);
}

void Editor::CloseWindow()
{
   // Called when closed via window manager action.

   delete this;
}

void Editor::DoOK()
{
   // Handle ok button.

   CloseWindow();
}

void Editor::DoOpen()
{
   SetTitle();
}

void Editor::DoSave()
{
   SetTitle();
}

void Editor::DoClose()
{
   // Handle close button.

   CloseWindow();
}
//=================================================


void listBox::showHelpWindow(){
  if (kVerbosity>5) cout << "showHelpWindow\n";
  
   Editor *ed = new Editor(this, 700, 400);
//    ed->LoadBuffer(editortxt1);
//    ed->AddBuffer(editortxt2);
//    ed->AddBuffer(editortxt3);
//    ed->AddBuffer(editortxt4);
//    ed->AddBuffer(editortxt5);
//    ed->AddBuffer(editortxt6);
   ed->Popup();
  
}

listBox::listBox(const TGWindow *p, UInt_t w, UInt_t h, TTree *hitSumm) :
   TGMainFrame(p, w, h), fSelHitSumm(0), fHistTrack(0), fHistSel(0), fHistSelMinY(0),fHistSelMaxY(0), fEventArrow(0), fHelpWindow(0), fSelHitList(0), fLastEntryInListBox(0),fFirstEntryInListBox(0)
{
  fCacheFile = new TFile("tempCache.root","RECREATE");
  
  fHitSumm=hitSumm;
  
  gStyle->SetOptFit();
  
  //Init "nothing to plot" screen
  fNothingToPlot = new TPaveText(.05,.05,.95,.95);
  
  //Init "no formula" screen
  fNoFormula = fNothingToPlot ;/*new TPaveText(.05,.05,.95,.95);
  fNoFormula->SetTextSize(0.09);
  fNoFormula->AddText("No formula");*/
  
  // Create a horizontal frame containing the selection cuts stuff
  TGHorizontalFrame *hFrameCuts = new TGHorizontalFrame(this, 150, 20, kFixedWidth);
  fSelCuts = new TGTextEntry(hFrameCuts,"flag==0");
  hFrameCuts->AddFrame(fSelCuts, new TGLayoutHints(kLHintsExpandX));
  
  TGTextButton *applyCutsButton = new TGTextButton(hFrameCuts, "Apply selection cuts");
  applyCutsButton->SetToolTipText("Click here to apply the selection cuts you made");
  applyCutsButton->Connect("Pressed()", "listBox", this, "applyCuts()");
  hFrameCuts->AddFrame(applyCutsButton, new TGLayoutHints());
  fSelCuts->Connect("ReturnPressed()", "listBox", this, "applyCuts()");
  
  
  TGTextButton *showHelpButton = new TGTextButton(hFrameCuts, "Show help");
  showHelpButton->SetToolTipText("Click here to show the help window");
  showHelpButton->Connect("Pressed()", "listBox", this, "showHelpWindow()");
  hFrameCuts->AddFrame(showHelpButton, new TGLayoutHints());
  
  
  
  AddFrame(hFrameCuts, new TGLayoutHints(kLHintsExpandX, 2, 2, 5, 1));

  
  // Create main frame (the one that contains the three columns:
  // the event list, the event display and the two extra information canvases
  TGHorizontalFrame *hframe1 = new TGHorizontalFrame(this);
  
  fCanvas  = new TRootEmbeddedCanvas("ec", hframe1, 720, 600);
  fCanvas->SetEditDisabled();
  
  fListBox = new TGListBox(hframe1, 89);
  
  fSelected = new TList;

  fListBox->Resize(100,600);
  fListBox->GetContainer()->Connect("DoubleClicked(TGFrame*, Int_t)", "listBox", this, "DoubleClicked(TGFrame*, Int_t)");
  
  fListBox->Connect("Selected(Int_t)",    "listBox", this, "drawItem(int)");
  
  fListBox->GetContainer()->Connect("CurrentChanged(TGFrame*)", "listBox", this, "selectionChanged(TGFrame*)");
  
  hframe1->AddFrame(fListBox, new TGLayoutHints(kLHintsExpandX));
  hframe1->AddFrame(fCanvas, new TGLayoutHints(kLHintsNormal));
  
  
  // Rightmost column
  TGVerticalFrame *vframe1 = new TGVerticalFrame(hframe1);
  TGVerticalFrame *hframe2 = new TGVerticalFrame(vframe1);
  fG1  = new TRootEmbeddedCanvas("g1", hframe2, 300, 100);
  fG3  = new TRootEmbeddedCanvas("g3", hframe2, 300, 200);
  hframe2->AddFrame(fG1);
  hframe2->AddFrame(fG3);
  vframe1->AddFrame(hframe2);
  
  //Bottom pannel in leftmost column
  TGVerticalFrame *vframe1_Low = new TGVerticalFrame(vframe1);
  // Formula stuff
  TGHorizontalFrame *hframeHistoFormula = new TGHorizontalFrame(vframe1_Low);
  fHistFormula = new TGTextEntry(vframe1_Low,"E0");
  hframeHistoFormula->AddFrame(fHistFormula, new TGLayoutHints(kLHintsExpandX));
  TGTextButton *applyFormulaButton = new TGTextButton(hframeHistoFormula, "Update formula");
  applyFormulaButton->SetToolTipText("Click here to update the histogram formula");
  applyFormulaButton->Connect("Pressed()", "listBox", this, "drawHistogram()");
  fHistFormula->Connect("ReturnPressed()", "listBox", this, "drawHistogram()");
  hframeHistoFormula->AddFrame(applyFormulaButton, new TGLayoutHints());
  vframe1_Low->AddFrame(hframeHistoFormula, new TGLayoutHints(kLHintsExpandX));
  
  // Plot scale and fit stuff
  TGHorizontalFrame *hframeScaleAndFit = new TGHorizontalFrame(vframe1_Low);
  fCheckLogHisto = new TGCheckButton(hframeScaleAndFit, "Log scale", 2);
  TGTextButton *fitButton = new TGTextButton(hframeScaleAndFit, "Fit panel",33);
  hframeScaleAndFit->AddFrame(fCheckLogHisto);
  hframeScaleAndFit->AddFrame(fitButton);
  vframe1_Low->AddFrame(hframeScaleAndFit);
  fG2  = new TRootEmbeddedCanvas("g2", vframe1_Low, 300, 100);
  fCheckLogHisto->Connect("Clicked()", "listBox", this, "HandleButtons()");
  fitButton->Connect("Pressed()", "listBox", this, "HandleButtons()");
  vframe1_Low->AddFrame(fG2, new TGLayoutHints(kLHintsExpandY));
  
  vframe1->AddFrame(vframe1_Low, new TGLayoutHints(kLHintsExpandY));
  hframe1->AddFrame(vframe1, new TGLayoutHints(kLHintsExpandY));
  
  //Add main frame
  AddFrame(hframe1);
  
  fCheckLogHisto->SetOn();
  
  
//   fCheckMulti = new TGCheckButton(this, "&Mutliple selection", 10);
//   AddFrame(fCheckMulti, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5));
//   fCheckMulti->Connect("Clicked()", "listBox", this, "HandleButtons()");

  
  // Create a horizontal frame containing button(s)
  TGHorizontalFrame *hframe = new TGHorizontalFrame(this, 150, 20, kFixedWidth);
  TGTextButton *exit = new TGTextButton(hframe, "&Exit ");
  exit->Connect("Pressed()", "listBox", this, "DoExit()");
  hframe->AddFrame(exit, new TGLayoutHints(kLHintsExpandX, 5, 5, 3, 4));
  AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 2, 2, 5, 1));

  // Set a name to the main frame   
  SetWindowName("List Box");
  MapSubwindows();

  // Initialize the layout algorithm via Resize()
  Resize(GetDefaultSize());
  
  // Map main frame
  MapWindow();

  //Create internal selected events list
  fSelHitList = new std::vector<listEntry_t>();
  
  // Build event list
  buildEventList(fSelCuts->GetText());
}


void listBox::DoubleClicked(TGFrame *f, Int_t btn)
{
  if (btn < 4) {
    if (kVerbosity>5) cout << "DOUBLE CLICK" << endl;
    Printf("Selected entries is: %d\n", fListBox->GetSelected());
    //drawItem(fListBox->GetSelected()-1);
  }
}

listBox::~listBox()
{
   // Clean up main frame...
   Cleanup();
   if (fSelected) {
      fSelected->Delete();
      delete fSelected;
   }
}

void listBox::HandleButtons()
{
  // Handle check button.
  Int_t id;
  TGButton *btn = (TGButton *) gTQSender;
  id = btn->WidgetId();

  printf("HandleButton: id = %d\n", id);

  if (id == 10)  
    fListBox->SetMultipleSelections(fCheckMulti->GetState());

  if (id == 2){ //Log scale checkbox
    drawHistogram();
    fG2->GetCanvas()->cd();
    if(fCheckLogHisto->GetState()>0)
      fG2->GetCanvas()->SetLogy(1);
    else
      fG2->GetCanvas()->SetLogy(0);

    fG2->GetCanvas()->Update();
  }
  
  if (id == 33){
    if(fHistSel!=0)
      fHistSel->FitPanel();
  }
}


void listBox::PrintSelected()
{
   // Writes selected entries in TList if multiselection.

   fSelected->Clear();

   if (fListBox->GetMultipleSelections()) {
      Printf("Selected entries are:\n");
      fListBox->GetSelectedEntries(fSelected);
      fSelected->ls();
   } else {
      Printf("Selected entries is: %d\n", fListBox->GetSelected());
      //fCan = new TCanvas("fCan");
      //fCan->Draw();
      drawItem(fListBox->GetSelected()-1);
   }
}


