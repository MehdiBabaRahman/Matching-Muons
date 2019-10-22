#include <iostream>
#include "Riostream.h"
#include <string>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <math.h>
#include "TTree.h"
#include "TBranch.h"
#include "TFrame.h"
#include "TCanvas.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TFile.h"
#include "TString.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "math.h"
#include "time.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include <algorithm>
#include "TLine.h"

using namespace std;

int Matching_muons_v1(){    

    Float_t selMu0_eta;
   // Float_t genMu0_eta;
    Float_t selMu0_phi;
    Float_t genMu0_phi;

    Float_t selMu1_eta;
    Float_t genMu1_eta;
    Float_t selMu1_phi;
    Float_t genMu1_phi;

    Float_t selMu2_eta;
    Float_t genMu2_eta;
    Float_t selMu2_phi;
    Float_t genMu2_phi;

    Float_t selMu3_eta;
    Float_t genMu3_eta;
    Float_t selMu3_phi;
    Float_t genMu3_phi;


    TFile *myFile = new TFile("ALL.root");
    TCanvas *cnv = new TCanvas();
    TH1F *sD_Lxy = new TH1F("dR", "dR", 100, -1E-4, 1E-2);
    sD_Lxy->SetXTitle("dR");
    sD_Lxy->SetYTitle("Number of Events");

    gDirectory->cd("cutFlowAnalyzerPXBL4PXFL3;1");
    gDirectory->pwd();
    TTree *myTree1 = nullptr;
    gDirectory->GetObject("Events;1",myTree1);

    //Number of entries
    int N = myTree1->GetEntries();
    cout << "Number of entries: " << N << endl;


    Float_t genMu0_eta;
    Int_t counter0, counter1 = 0 ;

    myTree1->SetBranchAddress("genMu0_eta", &genMu0_eta);
    myTree1->SetBranchAddress("selMu0_eta", &selMu0_eta);
    myTree1->SetBranchAddress("genMu1_eta", &genMu1_eta);
    myTree1->SetBranchAddress("selMu1_eta", &selMu1_eta);
    myTree1->SetBranchAddress("genMu2_eta", &genMu2_eta);
    myTree1->SetBranchAddress("selMu2_eta", &selMu2_eta);
    myTree1->SetBranchAddress("genMu3_eta", &genMu3_eta);
    myTree1->SetBranchAddress("selMu3_eta", &selMu3_eta);

    myTree1->SetBranchAddress("genMu0_phi", &genMu0_phi);
    myTree1->SetBranchAddress("selMu0_phi", &selMu0_phi);
    myTree1->SetBranchAddress("genMu1_phi", &genMu1_phi);
    myTree1->SetBranchAddress("selMu1_phi", &selMu1_phi);
    myTree1->SetBranchAddress("genMu2_phi", &genMu2_phi);
    myTree1->SetBranchAddress("selMu2_phi", &selMu2_phi);
    myTree1->SetBranchAddress("genMu3_phi", &genMu3_phi);
    myTree1->SetBranchAddress("selMu3_phi", &selMu3_phi);


    for(int ii; ii < N; ii++){


            myTree1->GetEntry(ii);
            
            Float_t a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, min;

            if (    
            selMu0_eta != -100 &&
            genMu0_eta != -100 && 
            selMu0_phi != -100 &&
            genMu0_phi != -100 &&

            selMu1_eta != -100 &&
            genMu1_eta != -100 &&
            selMu1_phi != -100 &&
            genMu1_phi != -100 &&

            selMu2_eta != -100 &&
            genMu2_eta != -100 &&
            selMu2_phi != -100 &&
            genMu2_phi != -100 &&

            selMu3_eta != -100 &&
            genMu3_eta != -100 &&
            selMu3_phi != -100 &&
            genMu3_phi != -100 ){

            

            a = ( (selMu0_eta - genMu0_eta) * (selMu0_eta - genMu0_eta) ) +  ( (selMu0_phi - genMu0_phi) * (selMu0_phi - genMu0_phi) ) ;
            b = ( (selMu0_eta - genMu1_eta) * (selMu0_eta - genMu1_eta) ) +  ( (selMu0_phi - genMu1_phi) * (selMu0_phi - genMu1_phi) ) ;
            c = ( (selMu0_eta - genMu2_eta) * (selMu0_eta - genMu2_eta) ) +  ( (selMu0_phi - genMu2_phi) * (selMu0_phi - genMu2_phi) ) ;
            d = ( (selMu0_eta - genMu3_eta) * (selMu0_eta - genMu3_eta) ) +  ( (selMu0_phi - genMu3_phi) * (selMu0_phi - genMu3_phi) ) ;


            e = ( (selMu1_eta - genMu0_eta) * (selMu1_eta - genMu0_eta) ) +  ( (selMu1_phi - genMu0_phi) * (selMu1_phi - genMu0_phi) ) ;
            f = ( (selMu1_eta - genMu1_eta) * (selMu1_eta - genMu1_eta) ) +  ( (selMu1_phi - genMu1_phi) * (selMu1_phi - genMu1_phi) ) ;
            g = ( (selMu1_eta - genMu2_eta) * (selMu1_eta - genMu2_eta) ) +  ( (selMu1_phi - genMu2_phi) * (selMu1_phi - genMu2_phi) ) ;
            h = ( (selMu1_eta - genMu3_eta) * (selMu1_eta - genMu3_eta) ) +  ( (selMu1_phi - genMu3_phi) * (selMu1_phi - genMu3_phi) ) ;

            i = ( (selMu2_eta - genMu0_eta) * (selMu2_eta - genMu0_eta) ) +  ( (selMu2_phi - genMu0_phi) * (selMu2_phi - genMu0_phi) ) ;
            j = ( (selMu2_eta - genMu1_eta) * (selMu2_eta - genMu1_eta) ) +  ( (selMu2_phi - genMu1_phi) * (selMu2_phi - genMu1_phi) ) ;
            k = ( (selMu2_eta - genMu2_eta) * (selMu2_eta - genMu2_eta) ) +  ( (selMu2_phi - genMu2_phi) * (selMu2_phi - genMu2_phi) ) ;
            l = ( (selMu2_eta - genMu2_eta) * (selMu2_eta - genMu2_eta) ) +  ( (selMu2_phi - genMu2_phi) * (selMu2_phi - genMu2_phi) ) ;

            m = ( (selMu3_eta - genMu0_eta) * (selMu3_eta - genMu0_eta) ) +  ( (selMu3_phi - genMu0_phi) * (selMu3_phi - genMu0_phi) ) ;
            n = ( (selMu3_eta - genMu1_eta) * (selMu3_eta - genMu1_eta) ) +  ( (selMu3_phi - genMu1_phi) * (selMu3_phi - genMu1_phi) ) ;
            o = ( (selMu3_eta - genMu2_eta) * (selMu3_eta - genMu2_eta) ) +  ( (selMu3_phi - genMu2_phi) * (selMu3_phi - genMu2_phi) ) ;
            p = ( (selMu3_eta - genMu3_eta) * (selMu3_eta - genMu3_eta) ) +  ( (selMu3_phi - genMu3_phi) * (selMu3_phi - genMu3_phi) ) ;
            min= std::min(std::min(std::min(std::min(std::min(std::min(std::min(std::min(std::min(std::min(std::min(std::min(std::min(a,b), c),d),e),f),g),h),k),l),m),n),o),p) ;
          //  cout << "minimum dR =  "<< sqrt(min) << endl;
            sD_Lxy->Fill(sqrt(min)); 

            }
            else {counter1 = counter1+1;}
             
            // else if ( 
            // selMu0_eta = -100 ||
            // genMu0_eta = -100 || 
            // selMu0_phi = -100 ||
            // genMu0_phi = -100 ||

            // selMu1_eta = -100 ||
            // genMu1_eta = -100 ||
            // selMu1_phi = -100 ||
            // genMu1_phi = -100 ||

            // selMu2_eta = -100 ||
            // genMu2_eta = -100 ||
            // selMu2_phi = -100 ||
            // genMu2_phi = -100 ||

            // selMu3_eta = -100 ||
            // genMu3_eta = -100 ||
            // selMu3_phi = -100 ||
            // genMu3_phi = -100 ){
            //        counter1 = counter1+1; 
            //          }



  


     }    

    cout << "not matched = " << counter1 << endl;   
    sD_Lxy->Draw();

    gStyle->SetOptFit(1);






    cnv->SaveAs("dR_min_Ntuples.pdf");

    return 0;
}

