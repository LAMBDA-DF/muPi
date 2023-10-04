#include <string.h>
#include <stdio.h>
#include "fitsio.h"

#include <iostream>
#include <sstream>

#include <sys/time.h>
#include <time.h>
#include <inttypes.h>
#include <fstream>
#include <unistd.h>
#include <getopt.h>    /* for getopt_long; standard getopt is in unistd.h */
#include <vector>
#include <algorithm>
#include <ctime>
#include <climits>
#include <cmath>
#include <iomanip>

#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
#endif

#include "globalConstants.h"

using namespace std;

int deleteFile(const char *fileName){
  cout << yellow;
  cout << "Will overwrite: " << fileName << endl << endl;
  cout << normal;
  return unlink(fileName);
}

bool fileExist(const char *fileName){
  ifstream in(fileName,ios::in);
  
  if(in.fail()){
    //cout <<"\nError reading file: " << fileName <<"\nThe file doesn't exist!\n\n";
    in.close();
    return false;
  }
  
  in.close();
  return true;
}

/*========================================================
  ASCII progress bar
==========================================================*/
void showProgress(unsigned int currEvent, unsigned int nEvent) {

  const int nProgWidth=50;

  if ( currEvent != 0 ) {
    for ( int i=0;i<nProgWidth+8;i++)
      cout << "\b";
  }

  double percent = (double) currEvent/ (double) nEvent;
  int nBars = (int) ( percent*nProgWidth );

  cout << " |";
  for ( int i=0;i<nBars-1;i++)
    cout << "=";
  if ( nBars>0 )
    cout << ">";
  for ( int i=nBars;i<nProgWidth;i++)
    cout << " ";
  cout << "| " << setw(3) << (int) (percent*100.) << "%";
  cout << flush;

}

bool readListFile(const string &inListFile, vector<string> &inFileList){
  
  if(gVerbosity){
    cout << bold << "Reading input list file.\n" << normal;
  }
  
  ifstream in(inListFile.c_str());
  
  if(!in.is_open()){
    cerr << "\nCan't open file: " << inListFile << endl << endl;
    return false;
  }
  
  while(!in.eof()){
    char line[kMaxLine];
    in.getline(line,kMaxLine);
    string lineS(line);
    
    size_t found = lineS.find_first_not_of(" \t\v\r\n#");
    if(found == string::npos) continue;
    
    unsigned int i=0;
    istringstream eventISS(&(lineS[found]));
    string aux;
    const unsigned int nCol = 1;
    while(eventISS >> aux){
      inFileList.push_back(aux);
      if(i>=nCol){
        i=-1;
        cout << aux << endl;
        return false;
      }
      ++i;
    }
    if(i!=nCol){
      cout << "Error: there is more than one column in the list file!\n\n";
      return false;
    }
    
  }
  
  in.close();
  
  return true;
}

void parseTime(const char *card, tm &t)
{  
  string timeS(card); 
  for(unsigned int i=0;i<timeS.size();++i) if(timeS[i]<48 || timeS[i]>57) timeS[i]=' ';
  istringstream timeISS(timeS);
  
  int auxYear;
  timeISS >> auxYear >> t.tm_mon >> t.tm_mday >> t.tm_hour >> t.tm_min >> t.tm_sec;
  t.tm_year = auxYear-1900;
  time_t tt = mktime (&t); //computes the day of the week
  t=*(localtime(&tt));
  
//   cout << asctime(&t) << endl;
  return;
}

int readTimeFromHeader(const char* fileName){
   fitsfile *fptr;         
  char card[FLEN_CARD]; 
  int status = 0,  nkeys, ii;  /* MUST initialize status */

  fits_open_file(&fptr, fileName, READONLY, &status);
  fits_get_hdrspace(fptr, &nkeys, NULL, &status);

  tm t;
  
  for (ii = 1; ii <= nkeys; ii++){ 
    fits_read_record(fptr, ii, card, &status); /* read keyword */
    
    if(strncmp(card,"UTSHUT",6) == 0){
      parseTime(card, t);
      cout << endl << asctime(&t) << "\t" << mktime (&t) << endl;
    }
  }
  fits_close_file(fptr, &status);

  if(status)          /* print any error messages */
    fits_report_error(stderr, status);
  return(status);
}

void printCopyHelp(const char *exeName, bool printFullHelp=false){
  
  if(printFullHelp){
    cout << bold;
    cout << endl;
    cout << "This program computes the median images from N input fit files.\n";
    cout << "It handles all the available HDUs. The HDUs in the output fit file\n";
    cout << "will be:\n";
    cout << " * float (32bits) for:  int8, int16, int32 and float input images\n";
    cout << " * double (64bits) for: double input images\n";
    cout << normal;
  }
  cout << "==========================================================================\n";
  cout << yellow;
  cout << "\nUsage:\n";
  cout << "  "   << exeName << " <input file 1> .. <input file N> -o <output filename>\n\n";
  cout << "or:\n";
  cout << "  "   << exeName << " -i <input list file > -o <output filename>\n\n";
  cout << "or:\n";
  cout << "  "   << exeName << " <input file 1> .. <input file N> -i <input list file > -o <output filename>\n";
  cout << "\nOptions:\n";
  cout << "  -q for quiet (no screen output)\n";
  cout << "  -s <HDU number> for processing a single HDU\n";
  cout << "  -m for MAD instead of Median\n";
  cout << "  -x for sigma from MAD instead of Median. sigma = 1.4826*MAD\n";
  cout << "  -a for stack (add all images) instead of Median\n\n";
  cout << normal;
  cout << blue;
  cout << "For any problems or bugs contact Javier Tiffenberg <javiert@fnal.gov>\n\n";
  cout << normal;
  cout << "==========================================================================\n\n";
}

string bitpix2TypeName(int bitpix){
  
  string typeName;
  switch(bitpix) {
      case BYTE_IMG:
          typeName = "BYTE(8 bits)";
          break;
      case SHORT_IMG:
          typeName = "SHORT(16 bits)";
          break;
      case LONG_IMG:
          typeName = "INT(32 bits)";
          break;
      case FLOAT_IMG:
          typeName = "FLOAT(32 bits)";
          break;
      case DOUBLE_IMG:
          typeName = "DOUBLE(64 bits)";
          break;
      default:
          typeName = "UNKNOWN";
  }
  return typeName;
}

void computeStack( vector< vector <double> > &vLinePix, vector<double> &vMedian ){

  int npix = vLinePix.size();
  vMedian.clear();
  vMedian.resize( npix );
  
  const int nElementsPerPix = vLinePix[0].size();
 
  #pragma omp parallel
  {
    #pragma omp for
    for(int c=0;c<npix;++c){
      double auxSum = 0;
      for(int i=0;i<nElementsPerPix;++i){
        auxSum += vLinePix[c][i];
      }
      vMedian[c] = auxSum; 
    }
  }

}

static bool isLargerThan(const int &x){
	return x>kOutlierThreshold;
}

void sortFilterAndComputeMedian( vector< vector <double> > &vLinePix, vector<double> &vMedian ){

  int npix = vLinePix.size();
  vMedian.clear();
  vMedian.resize( npix );
  
  
  #pragma omp parallel
  {
      #pragma omp for
      for(int c=0;c<npix;++c){
	  	vLinePix[c].erase(std::remove_if(vLinePix[c].begin(),vLinePix[c].end(), isLargerThan), vLinePix[c].end());
		
		const int nElementsPerPix = vLinePix[c].size();
		bool isOdd = !!(nElementsPerPix & 1);
		const int m = nElementsPerPix/2;
	  	
    	if( isOdd ){
	        std::nth_element(vLinePix[c].begin(), vLinePix[c].begin()+vLinePix[c].size()/2, vLinePix[c].end());
	        vMedian[c] = vLinePix[c][m];
	    }
	    else{
	        std::nth_element(vLinePix[c].begin(), vLinePix[c].begin()+vLinePix[c].size()/2, vLinePix[c].end());
	        const double highMedianVal = vLinePix[c][m];
	        std::nth_element(vLinePix[c].begin(), vLinePix[c].begin()+vLinePix[c].size()/2-1, vLinePix[c].end());
	        const double lowMedianVal  = vLinePix[c][m-1];
	        vMedian[c] = (highMedianVal+lowMedianVal)/2.;
	    }  
      }
	}
	
}

void sortFilterAndComputeMAD( vector< vector <double> > &vLinePix, vector<double> &vMAD, const string kMode = ""){
  
  const int npix = vLinePix.size();
  vector<double> vMedian( npix );
  sortFilterAndComputeMedian( vLinePix, vMedian );
  
  vMAD.clear();
  vMAD.resize( npix );
  
  
  #pragma omp parallel
  {
    #pragma omp for
    for(int c=0;c<npix;++c){
    	
		  const int nElementsPerPix = vLinePix[c].size();
		  bool isOdd = !!(nElementsPerPix & 1);
		  const int m = nElementsPerPix/2;

      for(int i=0;i<nElementsPerPix;++i){
        vLinePix[c][i] = fabs(vLinePix[c][i] - vMedian[c]);
      }

	    if( isOdd ){
        std::nth_element(vLinePix[c].begin(), vLinePix[c].begin()+vLinePix[c].size()/2, vLinePix[c].end());
        vMAD[c] = vLinePix[c][m];
	    }
	    else{
        std::nth_element(vLinePix[c].begin(), vLinePix[c].begin()+vLinePix[c].size()/2, vLinePix[c].end());
        const double highMADVal = vLinePix[c][m];
        std::nth_element(vLinePix[c].begin(), vLinePix[c].begin()+vLinePix[c].size()/2-1, vLinePix[c].end());
        const double lowMADVal  = vLinePix[c][m-1];
        vMAD[c] = (highMADVal+lowMADVal)/2.;
	    }
    }
    
  }

  if(kMode == "SIG") for(int c=0;c<npix;++c) vMAD[c] *= 1.4826;

}


bool checkInputFilesCompatibility(const vector<string> &inFileList){
  
  cout << "Start: checkInputFilesCompatibility\n";
  
  const int nFiles = inFileList.size();
  fitsfile *infptr_first;   /* FITS file pointers defined in fitsio.h */
  
  int status = 0;
  
  int nhdu_first = 0;
  
  const char* inF_first = inFileList.front().c_str();
  fits_open_file(&infptr_first, inF_first, READONLY, &status); /* Open the input file */
  if (status != 0) return(status);
  fits_get_num_hdus(infptr_first, &nhdu_first, &status);
  
  /* First check that all the files have the same number of hdus */
  for(int r=1; r<nFiles; ++r){
    
    /* Open the input file */
    fitsfile *infptr_other;
    int nhdu_other = 0;
    
    const char* inF_other = inFileList[r].c_str();
    fits_open_file(&infptr_other, inF_other, READONLY, &status);
    if (status != 0) return(status);
    fits_get_num_hdus(infptr_other, &nhdu_other, &status);
    fits_close_file(infptr_other, &status);
    
    if(nhdu_other != nhdu_first){
      fits_close_file(infptr_first,  &status);
      cerr << red << "\nInput files are not compatible. They don't have the same number of HDUs.\n\n" << normal;
      return false;
    }
  }
  
  
  for (int n=1; n<=nhdu_first; ++n)  /* Main loop through each extension */
  {
    /* get image dimensions and total number of pixels in image */
    int hdutype_first;
    int bitpix_first;
    int naxis_first = 0;
    long naxes_first[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    fits_movabs_hdu(infptr_first, n, &hdutype_first, &status);
    fits_get_img_param(infptr_first, 9, &bitpix_first, &naxis_first, naxes_first, &status);
    
    /* Loop on the input files */
    for(int r=1; r<nFiles; ++r){
      
      /* Open the input file */
      const char* inF_other = inFileList[r].c_str();
      fitsfile *infptr_other;
      int hdutype_other;
      int bitpix_other;
      int naxis_other = 0;
      long naxes_other[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
      
      fits_open_file(&infptr_other, inF_other, READONLY, &status);
      if (status != 0) return(status);
      fits_movabs_hdu(infptr_other, n, &hdutype_other, &status);
      fits_get_img_param(infptr_other, 9, &bitpix_other, &naxis_other, naxes_other, &status);
      
      if(bitpix_other!=bitpix_first){
        fits_close_file(infptr_first,  &status);
        fits_close_file(infptr_other, &status);
        cerr << red << "\nInput files are not compatible. At least one of the images has a different data type.\n\n" << normal;
        return false;
      }
      
      if(naxis_other!=naxis_first){
        fits_close_file(infptr_first,  &status);
        fits_close_file(infptr_other, &status);
        cerr << red << "\nInput files are not compatible. At least one of the images has a different number of HDUs.\n\n" << normal;
        return false;
      }
      
      for(int i=0; i<9; ++i){
        if(naxes_other[i]!=naxes_first[i]){
          fits_close_file(infptr_first,  &status);
          fits_close_file(infptr_other, &status);
          cerr << red << "\nInput files are not compatible. At least one of the images has a different size.\n\n" << normal;
          return false;
        }
      }
      
      fits_close_file(infptr_other, &status);
    }
    
  }
  
  fits_close_file(infptr_first,  &status);
  
  cout << "Start: checkInputFilesCompatibility\n";
  return true;
}

int copyStructure(const vector<string> inFileList, const char *outF, const int singleHdu){
    
  fitsfile  *outfptr; /* FITS file pointers defined in fitsio.h */
  fitsfile *infptr;   /* FITS file pointers defined in fitsio.h */
  
  int status = 0;
  int single = 0;
  
  int hdutype, bitpix, naxis = 0, nkeys;
  int nhdu = 0;
  long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
  long totpix = 0;
  //char card[81];
  
  ostringstream fileStructSummary;
  
  const char* inF = inFileList.front().c_str();
  fits_open_file(&infptr, inF, READONLY, &status); /* Open the input file */
  if (status != 0) return(status);
  
  fits_get_num_hdus(infptr, &nhdu, &status);
  
  if (singleHdu>0){
    if(singleHdu > nhdu){
      fits_close_file(infptr,  &status);
      cerr << red << "\nError: the file does not have the required HDU!\n\n" << normal;
      return -1000;
    }
    single = 1; /* Copy only a single HDU if a specific extension was given */
  }
  
  if(single)
    fileStructSummary << "The output file will contain HDU" << singleHdu << " of " << nhdu << " availables in the input files."<< endl;
  else
    fileStructSummary << "The output file will contain all the " << nhdu << " HDUs availables in the input files.\n";
    
  
  fits_create_file(&outfptr, outF, &status);/* Create the output file */
  if (status != 0) return(status);
  
  
  fileStructSummary << bold << "HDU   hdutype  #Axis  #Cols  #Rows   IN_datatype      OUT_datatype\n" << normal;
// HDU  hdutype #Axis #Cols #Rows datatype  
  for (int n=1; n<=nhdu; ++n)  /* Main loop through each extension */
  {

    if (single) n = singleHdu;
    
    /* get image dimensions and total number of pixels in image */
    fits_movabs_hdu(infptr, n, &hdutype, &status);
    for (int i = 0; i < 9; ++i) naxes[i] = 1;
    fits_get_img_param(infptr, 9, &bitpix, &naxis, naxes, &status);
    totpix = naxes[0] * naxes[1] * naxes[2] * naxes[3] * naxes[4] * naxes[5] * naxes[6] * naxes[7] * naxes[8];
    fileStructSummary  << setw (3) << n << "  "  << setw (8) << hdutype << "  " << setw (5) << naxis << "  " << setw (5) << naxes[0] << "  " << setw (5) << naxes[1] << "  " << setw (15) << bitpix2TypeName(bitpix) << "  " << setw (15) << bitpix2TypeName(abs(bitpix) == 64? -64 : -32);
    
    //continue;
    if (hdutype != IMAGE_HDU || naxis == 0 || totpix == 0){
      /* just copy tables and null images */
      fits_copy_hdu(infptr, outfptr, 0, &status);
      if (status != 0) return(status);
      fileStructSummary << magenta << "<- Not an image HDU" << normal;
    }
    else{
      /* create output image with the same size as the input image*/
      if(abs(bitpix) == 64) fits_create_img(outfptr, -64, naxis, naxes, &status);
      else{
        fits_create_img(outfptr, -32, naxis, naxes, &status);
      }
      if (status != 0) return(status);

      /* copy the relevant keywords (not the structural keywords) */
      fits_get_hdrspace(infptr, &nkeys, NULL, &status); 
//       for (int i = 1; i <= nkeys; ++i) {
//         fits_read_record(infptr, i, card, &status);
//         if (fits_get_keyclass(card) > TYP_CMPRS_KEY) fits_write_record(outfptr, card, &status);
//       }
      
    }
    fileStructSummary << endl;
    
    /* quit if only copying a single HDU */
    if (single) break;
  }
  fits_close_file(infptr, &status);
  fits_close_file(outfptr,  &status);

  if(gVerbosity){
    cout << bold << "Files structure summary:\n" << normal;
    cout << fileStructSummary.str();
    cout << green << "Structure copied.\n\n" << normal;
  }
  return status;
}


/* Compute the median image and write it to an existing output file 
 * that has the right structure (created by copyStructure function */
int computeMedianImages(const vector<string> inFileList, const char *outF, const int singleHdu, const string kMode){
  
  int status = 0;
  int single = 0;
  
  int nhdu = 0;
  const int nFiles = inFileList.size();
  
  /* Open the output file */
  fitsfile  *outfptr; /* FITS file pointers defined in fitsio.h */
  fits_open_file(&outfptr, outF, READWRITE, &status);
  if (status != 0) return(status);
  
  fits_get_num_hdus(outfptr, &nhdu, &status);
  
  if (singleHdu>0){
    single = 1; /* Copy only a single HDU if a specific extension was given */
  }
  
  /* Open the input files */
  vector<fitsfile *> infPtrs(nFiles);
//   for(int f=0; f<nFiles; ++f){
//     const char* inF = inFileList[f].c_str();
//     fits_open_file(&(infPtrs[f]), inF, READONLY, &status);
//     if (status != 0) return(status);
//   }
  
  const int kReadNLines = 100;
  const int nHDUsToProcess = (single>0)? 1 : nhdu;
  for (int n=1; n<=nhdu; ++n)  /* Main loop through each extension */
  {
    int nOut = n;
    if (single){
      n = singleHdu;
      nOut = 1;
    }
    
    /* get image dimensions and total number of pixels in image */
    int hdutype, bitpix, bytepix, naxis = 0, anynul;
    long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    double nulval = 0.;
    fits_movabs_hdu(outfptr, nOut, &hdutype, &status);
    for (int i = 0; i < 9; ++i) naxes[i] = 1;
    fits_get_img_param(outfptr, 9, &bitpix, &naxis, naxes, &status);
    long totpix = naxes[0] * naxes[1];
    
    /* Don't try to process data if the hdu is empty */    
    if (hdutype != IMAGE_HDU || naxis == 0 || totpix == 0){
      if(single) break;
      continue;
    }
    
    bytepix = abs(bitpix) / 8;
    if(bytepix!=4 && bytepix!=8) return -1000;
    
    char* outArray = new char[totpix*bytepix];
    
    for(int i=0;i<totpix*bytepix;++i) outArray[i] = 0;
    
    const long npixStart = naxes[0]*kReadNLines;
    long npix = npixStart;
    long pixelsLeft = totpix;
    
    if(gVerbosity){
      cout << bold << "\rProcessing HDU: " << n << normal << flush;
      showProgress(0,1);
    }
    
    const int nLoops = (const int)std::ceil(naxes[1]*1.0/kReadNLines);
    for(int l=0;l<nLoops;++l){
     
      if(l%10 && gVerbosity){
        if(single) showProgress(l,nLoops);
        else       showProgress((n-1)*nLoops+l,nLoops*nHDUsToProcess);
      }
      
      if(npix>pixelsLeft) npix=pixelsLeft;
      vector< vector <double> > vLinePix(npix, vector<double>(nFiles) );
      
      long firstpix[2];
      firstpix[0] = 1;
      firstpix[1] = kReadNLines*l+1;
      
      /* Loop on the input files */
      for(int r=0; r<nFiles; ++r){
        
        const char* inF = inFileList[r].c_str();
        fits_open_file(&(infPtrs[r]), inF, READONLY, &status);
        if (status != 0) return(status);
        
        double* lArray = new double[npix];
        
        /* Open the input file */
        fits_movabs_hdu(infPtrs[r], n, &hdutype, &status);
        
        /* Read the images as doubles, regardless of actual datatype. */
        fits_read_pix(infPtrs[r], TDOUBLE, firstpix, npix, &nulval, lArray, &anynul, &status);

        for(int c=0;c<npix;++c) vLinePix[c][r] = lArray[c];

        delete[] lArray;
        fits_close_file(infPtrs[r], &status);
      }
      
      vector<double> vMedian;
      
      if(kMode == "Median")
        sortFilterAndComputeMedian( vLinePix, vMedian );
      else if(kMode == "MAD")
        sortFilterAndComputeMAD( vLinePix, vMedian );
      else if(kMode == "SIG")
        sortFilterAndComputeMAD( vLinePix, vMedian, kMode );
      else if(kMode == "STACK")
        computeStack( vLinePix, vMedian );
      else{
	cerr << "Invalid mode!\n";
	return -1001;
      }
      
      if(bytepix==4){
        for(int c=0;c<npix;++c) reinterpret_cast<float *> (outArray)[l*npixStart + c] = vMedian[c];
      }
      else if(bytepix==8){
        for(int c=0;c<npix;++c) reinterpret_cast<double *>(outArray)[l*npixStart + c] = vMedian[c];
      }
      pixelsLeft-=npix;
    }
    
    fits_set_bscale(outfptr, 1., 0., &status);/* Set rescale to pixel = 1*pixel_value + 0 for output file*/
    if(bytepix==4)      fits_write_img(outfptr, TFLOAT, 1, totpix, reinterpret_cast<float *>(outArray), &status);
    else if(bytepix==8) fits_write_img(outfptr, TDOUBLE, 1, totpix, reinterpret_cast<double *>(outArray), &status);
    delete[] outArray;
    
    /* quit if only copying a single HDU */
    if (single) break;
  }
  
  /* Close the input files */
//   for(int f=0; f<nFiles; ++f){
//     fits_close_file(infPtrs[f], &status);
//   }
  
  /* Close the output file */
  fits_close_file(outfptr,  &status);
  
  if(gVerbosity){
    showProgress(1,1);
    cout << endl;
  }
  if(kMode == "Median")
    cout << green << "Median file computed.\n\n" << normal;
  else if(kMode == "MAD")
    cout << green << "MAD file computed.\n\n" << normal;
  
  return(status);
}

void checkArch(){
  if(sizeof(float)*CHAR_BIT!=32 || sizeof(double)*CHAR_BIT!=64){
    cout << red;
    cout << "\n ========================================================================================\n";
    cout << "   WARNING: the size of the float and double variables is non-standard in this computer.\n";
    cout << "   The program may malfunction or produce incorrect results\n";
    cout << " ========================================================================================\n";
    cout << normal;
  }
}

int processCommandLineArgs(const int argc, char *argv[], int &singleHdu, string &kMode, vector<string> &inFileList, string &outFile){
  
  if(argc == 1) return 1;
  
  bool outFileFlag = false;
  bool inListFileFlag = false;
  string inListFile = "";
  kMode="Median";
  int opt=0;
  while ( (opt = getopt(argc, argv, "amxi:o:s:qQhH?")) != -1) {
    switch (opt) {
    case 'o':
      if(!outFileFlag){
        outFile = optarg;
        outFileFlag = true;
      }
      else{
        cerr << red << "\nError, can not set more than one output file!\n\n" << normal;
        return 2;
      }
      break;
    case 'i':
      if(!inListFileFlag){
        inListFile = optarg;
        inListFileFlag = true;
      }
      else{
        cerr << red << "\nError, can not set more than one input list file!\n\n" << normal;
        return 2;
      }
      break;
    case 's':
      if(singleHdu<0){
        singleHdu = atoi(optarg);
      }
      else{
        cerr << red << "\nError, can not set more than one HDU!\n\n"  << normal;
        return 2;
      }
      break;
    case 'm':
      kMode = "MAD";
      break;
    case 'x':
      kMode = "SIG";
      break;
    case 'a':
      kMode = "STACK";
      break;
    case 'Q':
    case 'q':
      gVerbosity = 0;
      break;
    case 'h':
    case 'H':
    default: /* '?' */
      return 1;
    }
  }
  
  if(!outFileFlag){
    cerr << red << "\nError: output filename missing.\n" << normal;
    return 2;
  }

  inFileList.clear();
  
  if(inListFileFlag){
    cout << inListFile << endl;
    bool allOk = readListFile(inListFile, inFileList);
    if(!allOk){
      cerr << red << "\nError reading input list file!\n\n" << normal;
      return 1;
    }
  }
  
  for(int i=optind; i<argc; ++i){
    inFileList.push_back(argv[i]);
    if(!fileExist(argv[i])){
      cout << red << "\nError reading input file: " << argv[i] <<"\nThe file doesn't exist!\n\n" << normal;
      return 1;
    }
  }
  
  if(inFileList.size()==0){
    cerr << red << "Error: no input file(s) provided!\n\n" << normal;
    return 1;
  }
  
  return 0;
}

int main(int argc, char *argv[])
{
  
  checkArch(); //Check the size of the double and float variables.
  
  time_t start,end;
  double dif;
  time (&start);
  
  string outFile;
  vector<string> inFileList;
  int singleHdu=-1;
  string kMode;
  
  int returnCode = processCommandLineArgs( argc, argv, singleHdu, kMode, inFileList, outFile);
  if(returnCode!=0){
    if(returnCode == 1) printCopyHelp(argv[0],true);
    if(returnCode == 2) printCopyHelp(argv[0]);
    return returnCode;
  }
  
  if(gVerbosity){
    cout << endl << bold << "Mode: " << kMode << normal << endl;
    cout << bold << "\nWill read the following files:\n" << normal;
    for(unsigned int i=0; i<inFileList.size();++i) cout << "\t" << inFileList[i] << endl;
    cout << bold << "\nThe output will be saved in the file:\n\t" << normal << outFile << endl;
  }
  
  /* Check if input images are compatible in size, exposure time, etc.. */
  bool imagesAreCompatible = checkInputFilesCompatibility(inFileList);
  if(!imagesAreCompatible){
    cerr << red << "Will not continue\n\n" << normal;
    return -1;
  }
  
  /* Overwrite the output file if it already exist */
  if(fileExist(outFile.c_str())){
    cout << yellow << "\nThe output file exist. " << normal;
    deleteFile(outFile.c_str());
  }
  
  // readTimeFromHeader(argv[1]);
  
  /* Do the actual processing */
  int status = copyStructure( inFileList,  outFile.c_str(), singleHdu);
  if (status != 0){ 
    fits_report_error(stderr, status);
    return status;
  }
  status = computeMedianImages( inFileList,  outFile.c_str(), singleHdu, kMode);
  if (status != 0){ 
    fits_report_error(stderr, status);
    return status;
  }
  
  /* Report */
  time (&end);
  dif = difftime (end,start);
  if(gVerbosity) cout << green << "All done!\n" << bold << "-> It took me " << dif << " seconds to do it!\n\n" << normal;

  return status;
}
