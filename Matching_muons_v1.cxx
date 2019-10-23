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

   



    Int_t counter1;
    counter1 = 0;
    Int_t counter2;
    counter2 =0;
    Int_t counter3;
    counter3 =0;
    Int_t counter4;
    counter4 = 0;
    Int_t counter5;
    counter5 = 0;   

    Int_t counter6;
    counter6 = 0; 
    Int_t counter7;
    counter7 = 0; 
    Int_t counter8;
    counter8 = 0; 
    Int_t counter9;
    counter9 = 0; 
    Int_t counter10;
    counter10 = 0;
    Int_t counter11;
    counter11 = 0; 
    Int_t counter12;
    counter12 = 0; 
    Int_t counter13;
    counter13 = 0; 
    Int_t counter14;
    counter14 = 0; 
    Int_t counter15;
    counter15 = 0; 
    Int_t counter16;
    counter16 = 0; 
    Int_t counter17;
    counter17 = 0; 
    Int_t counter18;
    counter18 = 0; 
    Int_t counter19;
    counter19 = 0;
    Int_t counter20;
    counter20 = 0;



    TFile *myFile = new TFile("ALL.root");

    TCanvas *cnv1 = new TCanvas();
    TH1F *sD_Lxy1 = new TH1F("dR1", "dR1", 100, -1E-4, 5E-2);

    TCanvas *cnv2 = new TCanvas();
    TH1F *sD_Lxy2 = new TH1F("dR2", "dR2", 100, -1E-4, 5E-2);

    TCanvas *cnv3 = new TCanvas();
    TH1F *sD_Lxy3 = new TH1F("dR3", "dR3", 100, -1E-4, 5E-2);

    TCanvas *cnv4 = new TCanvas();
    TH1F *sD_Lxy4 = new TH1F("dR4", "dR4", 100, -1E-4, 5E-2);

    sD_Lxy1->SetXTitle("dR (best match for selMu0)");
    sD_Lxy1->SetYTitle("Number of Events");
    sD_Lxy2->SetXTitle("dR (best match for selMu1)");
    sD_Lxy2->SetYTitle("Number of Events ");
    sD_Lxy3->SetXTitle("dR (best match for selMu2)");
    sD_Lxy3->SetYTitle("Number of Events");
    sD_Lxy4->SetXTitle("dR(best match for selMu3)");
    sD_Lxy4->SetYTitle("Number of Events");

    gDirectory->cd("cutFlowAnalyzerPXBL4PXFL3;1");
    gDirectory->pwd();
    TTree *myTree1 = nullptr;
    gDirectory->GetObject("Events;1",myTree1);

    //Number of entries
    int N = myTree1->GetEntries();
    cout << "Number of entries: " << N << endl;

   
    Float_t genMu0_eta;



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



    for(int ii = 0; ii < N; ii++){


            myTree1->GetEntry(ii);
            
            Float_t a, b, c, d, min1;
            min1= min(min(min(a,b), c),d);
          
            if (    
            abs(selMu0_eta) <= 2.4 && 
            abs(genMu0_eta) <= 2.4 && 
            selMu0_phi != -100 &&
            genMu0_phi != -100 &&
            
            abs(genMu1_eta) <= 2.4 &&
            genMu1_phi != -100 &&

            
            abs(genMu2_eta) <= 2.4 &&
            genMu2_phi != -100 &&

            
            abs(genMu3_eta) <= 2.4 &&   
            genMu3_phi != -100 ){

            
            
            a = ( (selMu0_eta - genMu0_eta) * (selMu0_eta - genMu0_eta) ) +  ( (selMu0_phi - genMu0_phi) * (selMu0_phi - genMu0_phi) ) ;
            b = ( (selMu0_eta - genMu1_eta) * (selMu0_eta - genMu1_eta) ) +  ( (selMu0_phi - genMu1_phi) * (selMu0_phi - genMu1_phi) ) ;
            c = ( (selMu0_eta - genMu2_eta) * (selMu0_eta - genMu2_eta) ) +  ( (selMu0_phi - genMu2_phi) * (selMu0_phi - genMu2_phi) ) ;
            d = ( (selMu0_eta - genMu3_eta) * (selMu0_eta - genMu3_eta) ) +  ( (selMu0_phi - genMu3_phi) * (selMu0_phi - genMu3_phi) ) ;
            
            
            sD_Lxy1->Fill(sqrt(min1));

            }

            if ( min1 == a) {
            counter1 = counter1+1;
            }

            if (min1 == b) {
            counter2 = counter2+1;
            }

            if (min1 == c) {
            counter3 = counter3+1;
            
            }

            if (min1 == d) {
            counter4 = counter4+1;
            
            }

            else if ( min1 != a && min1 != b && min1 != c && min1 != d) { 
            counter5 = counter5+1;
            }



            Float_t e, f, g, h, min2;
            min2 = std::min(std::min(std::min(e,f), g),h);
    
            if (    
          
            abs(genMu0_eta) <= 2.4 && 
            genMu0_phi != -100 &&

            abs(selMu1_eta) <= 2.4 &&
            abs(genMu1_eta) <= 2.4 &&
            selMu1_phi != -100 &&
            genMu1_phi != -100 &&

 
            abs(genMu2_eta) <= 2.4 &&
            genMu2_phi != -100 &&

            abs(genMu3_eta) <= 2.4 &&
            genMu3_phi != -100 ){

            

            e = ( (selMu1_eta - genMu0_eta) * (selMu1_eta - genMu0_eta) ) +  ( (selMu1_phi - genMu0_phi) * (selMu1_phi - genMu0_phi) ) ;
            f = ( (selMu1_eta - genMu1_eta) * (selMu1_eta - genMu1_eta) ) +  ( (selMu1_phi - genMu1_phi) * (selMu1_phi - genMu1_phi) ) ;
            g = ( (selMu1_eta - genMu2_eta) * (selMu1_eta - genMu2_eta) ) +  ( (selMu1_phi - genMu2_phi) * (selMu1_phi - genMu2_phi) ) ;
            h = ( (selMu1_eta - genMu3_eta) * (selMu1_eta - genMu3_eta) ) +  ( (selMu1_phi - genMu3_phi) * (selMu1_phi - genMu3_phi) ) ;
            sD_Lxy2->Fill(sqrt(min2)); 
            
            }


            
            if ( min2 == e) {
            counter6 = counter6+1;
            }

            if (min2 == f) {
            counter7 = counter7+1;
            }

            if (min2 == g) {
            counter8 = counter8+1;
            
            }

            if (min2 == h) {
            counter9 = counter9+1;
            
            }

            else if ( min2 != e && min2 != f && min2 != g && min2 != h) { 
            counter10 = counter10+1;
            }            
            
            

            Float_t i, j, k, l, min3;
            min3 = std::min(std::min(std::min(i,j), k),l);
          

            if (    

            abs(genMu0_eta) <= 2.4 && 
            genMu0_phi != -100 &&

    
            abs(genMu1_eta) <= 2.4 &&
            genMu1_phi != -100 &&

            abs(selMu2_eta) <= 2.4 &&
            abs(genMu2_eta) <= 2.4 &&
            selMu2_phi != -100 &&
            genMu2_phi != -100 &&
    
            abs(genMu3_eta) <= 2.4 &&  
            genMu3_phi != -100 ){

            

            
           


            i = ( (selMu2_eta - genMu0_eta) * (selMu2_eta - genMu0_eta) ) +  ( (selMu2_phi - genMu0_phi) * (selMu2_phi - genMu0_phi) ) ;
            j = ( (selMu2_eta - genMu1_eta) * (selMu2_eta - genMu1_eta) ) +  ( (selMu2_phi - genMu1_phi) * (selMu2_phi - genMu1_phi) ) ;
            k = ( (selMu2_eta - genMu2_eta) * (selMu2_eta - genMu2_eta) ) +  ( (selMu2_phi - genMu2_phi) * (selMu2_phi - genMu2_phi) ) ;
            l = ( (selMu2_eta - genMu3_eta) * (selMu2_eta - genMu3_eta) ) +  ( (selMu2_phi - genMu3_phi) * (selMu2_phi - genMu3_phi) ) ;
            sD_Lxy3->Fill(sqrt(min3));
            }
            
            if ( min3 == i) {
            counter11 = counter11+1;
            }

            if (min3 == j) {
            counter12 = counter12+1;
            }

            if (min3 == k) {
            counter13 = counter13+1;
            
            }

            if (min3 == l) {
            counter14 = counter14+1;
            
            }

            else if ( min3 != i && min3 != j && min3 != k && min3 != l) { 
            counter15 = counter15+1;
            }            
                        
            


           

            Float_t m, n, o, p, min4;
            min4 = std::min(std::min(std::min(m,n), o),p);

            if (    

            abs(genMu0_eta) <= 2.4 && 
            genMu0_phi != -100 &&


            abs(genMu1_eta) <= 2.4 &&
            genMu1_phi != -100 &&

  
            abs(genMu2_eta) <= 2.4 &&
            genMu2_phi != -100 &&

            abs(selMu3_eta) <= 2.4 &&
            abs(genMu3_eta) <= 2.4 &&
            selMu3_phi != -100 &&
            genMu3_phi != -100 ){

            



            m = ( (selMu3_eta - genMu0_eta) * (selMu3_eta - genMu0_eta) ) +  ( (selMu3_phi - genMu0_phi) * (selMu3_phi - genMu0_phi) ) ;
            n = ( (selMu3_eta - genMu1_eta) * (selMu3_eta - genMu1_eta) ) +  ( (selMu3_phi - genMu1_phi) * (selMu3_phi - genMu1_phi) ) ;
            o = ( (selMu3_eta - genMu2_eta) * (selMu3_eta - genMu2_eta) ) +  ( (selMu3_phi - genMu2_phi) * (selMu3_phi - genMu2_phi) ) ;
            p = ( (selMu3_eta - genMu3_eta) * (selMu3_eta - genMu3_eta) ) +  ( (selMu3_phi - genMu3_phi) * (selMu3_phi - genMu3_phi) ) ;
            sD_Lxy4->Fill(sqrt(min4));
            }

            if ( min4 == m) {
            counter16 = counter16+1;
            }

            if (min4 == n) {
            counter17 = counter17+1;
            }

            if (min4 == o) {
            counter18 = counter18+1;
            
            }

            if (min4 == p) {
            counter19 = counter19+1;
            
            }

            else if ( min4 != m && min4 != n && min4 != o && min4 != p) { 
            counter20 = counter20+1;
            }            
            
            
            

    
            
             
    }

    cout << "leading sel muon matched to leading gen moun:  " << counter1 << endl;
    cout << "leading sel muon matched to sub-leading gen moun:  " << counter2 << endl;
    cout << "leading sel muon matched to 3rd gen moun:  " << counter3 << endl;
    cout << "leading sel muon matched to 4th gen moun:  " << counter4 << endl;
    cout << "leading sel moun not matched: " << counter5 << endl;

    cout << "sub-leading sel muon matched to leading gen moun:  " << counter6 << endl;
    cout << "sub-leading sel muon matched to sub-leading gen moun:  " << counter7 << endl;
    cout << "sub-leading sel muon matched to 3rd gen moun:  " << counter8 << endl;
    cout << "sub-leading sel muon matched to 4th gen moun:  " << counter9 << endl;
    cout << "sub-leading sel moun not matched: " << counter10 << endl;


    cout << "3rd sel muon matched to leading gen moun:  " << counter11 << endl;
    cout << "3rd sel muon matched to sub-leading gen moun:  " << counter12 << endl;
    cout << "3rd sel muon matched to 3rd gen moun:  " << counter13 << endl;
    cout << "3rd sel muon matched to 4th gen moun:  " << counter14 << endl;
    cout << "3rd sel moun not matched: " << counter15 << endl;


    cout << "4th sel muon matched to leading gen moun:  " << counter16 << endl;
    cout << "4th sel muon matched to sub-leading gen moun:  " << counter17 << endl;
    cout << "4th sel muon matched to 3rd gen moun:  " << counter18 << endl;
    cout << "4th sel muon matched to 4th gen moun:  " << counter19 << endl;
    cout << "4th sel moun not matched: " << counter20 << endl;

    cnv1->cd();   
    sD_Lxy1->Draw();
    cnv2->cd();
    sD_Lxy2->Draw();
    cnv3->cd();
    sD_Lxy3->Draw();
    cnv4->cd();
    sD_Lxy4->Draw();

    gStyle->SetOptFit(1);






    cnv1->SaveAs("dR1_min_Ntuples.pdf");
    cnv2->SaveAs("dR2_min_Ntuples.pdf");
    cnv3->SaveAs("dR3_min_Ntuples.pdf");
    cnv4->SaveAs("dR4_min_Ntuples.pdf");

    return 0;


}



            // if (min1 == a) {
            // counter1 = counter1+1;
            // }

            // if (min1 == b) {
            // counter2 = counter2+1;
            // }

            // if (min1 == c) {
            // counter3 = counter3+1;
            
            // }

            // if (min1 == d) {
            // counter4 = counter4+1;
            
            // }

            // else {
            // counter5 = counter5+1;
            // }



    // cout << "leading sel muon matched to leading gen moun:  " << counter1 << endl;
    // cout << "leading sel muon matched to sub-leading gen moun:  " << counter2 << endl;
    // cout << "leading sel muon matched to 3rd gen moun:  " << counter3 << endl;
    // cout << "leading sel muon matched to 4th gen moun:  " << counter4 << endl;
    // cout << "leading sel moun not matched: " << counter5 << endl;

    // cout << "total not matched:" << counter6 << endl;

