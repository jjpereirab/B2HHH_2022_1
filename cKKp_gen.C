#include <iostream>
#include <string>
#include <stdio.h>

//--------- DECLARING DATASET VARIABLES
TChain *B2HHH = new TChain("DecayTree");
Double_t B_FlightDistance;          Double_t B_VertexChi2;
Int_t H1_isMuon;         Int_t H2_isMuon;          Int_t H3_isMuon;
Int_t H1_Charge;         Int_t H2_Charge;          Int_t H3_Charge;
Double_t H1_IPChi2;      Double_t H2_IPChi2;       Double_t H3_IPChi2;
Double_t H1_PX;          Double_t H2_PX;           Double_t H3_PX;
Double_t H1_PY;          Double_t H2_PY;           Double_t H3_PY;
Double_t H1_PZ;          Double_t H2_PZ;           Double_t H3_PZ;
Double_t H1_ProbPi;      Double_t H2_ProbPi;       Double_t H3_ProbPi;
Double_t H1_ProbK;       Double_t H2_ProbK;        Double_t H3_ProbK;
int nentries, nbytes, i;  FILE *fp;

//--------- LITERATURE MASSES
double B_m  = 5279.38;   // B meson mass in MeV
double Pi_m = 139.57039; // MeV
double K_m  = 493.677;   // MeV
double D0m  = 1864.84;   // MeV

//--------- ADAPTATIVE BINNING OBJECTS
struct Bin{
    bool binOn = false;
    double xmin = 0.0, xmax = 0.0, ymin = 0.0, ymax = 0.0;
};
vector<Bin> vecBin(500);

void setBin(int i, bool Active,double x, double X, double y, double Y){
    vecBin[i].binOn = Active;
    vecBin[i].xmin = x;   vecBin[i].xmax = X;
    vecBin[i].ymin = y;   vecBin[i].ymax = Y;
}

//--------- CONTROL PANEL
double D0_thresh  = 30.;      double B_thresh   = 60.;   // D0 exclusion and B mass window thresholds
double above_this = 0.75;     double below_this = 0.45;   // PARTICLE ID EXCLUSION
double H1_pi = below_this;    double H1_ka = above_this;


//--------- HISTOGRAM DEFINITIONS
int dpxbins = 14;   int dpybins = 14;     //-------- DALITZ PLOT BINS
int infl    = 500;  int supl    = 6000;   //-------- range for normal dalitz plots

TH2F *whole         = new TH2F("whole","whole sym; Mka^{2}; Mpi^{2}",dpxbins,0,35, dpybins, 0, 35);
TH1F *h_Bmass_raw   = new TH1F("h_Bmass_raw","Bmass; m(MeV); Events",60,5000,5500);
TH1F *h_Bmass       = new TH1F("h_Bmass","Bmass; m(MeV); Events",60,5000,5500);
TH1F *h_PT          = new TH1F("h_PT","h_PT ; Pt(MeV); Events",60,1000,8000);

TH2F *h_Bplus       = new TH2F("h_Bplus"    ,"B^{+}; Mka^{2}; Mpi^{2}",dpxbins,0,35, dpybins, 0, 35);
TH2F *h_Bminus      = new TH2F("h_Bminus"   ,"B^{-}; Mka^{2}; Mpi^{2}",dpxbins,0,35, dpybins, 0, 35);
TH1F *h_Mpi_plus    = new TH1F("h_Mpi_plus" ,"Mpi ; m(MeV); Events",60,infl,supl);
TH1F *h_Mpi_minus   = new TH1F("h_Mpi_minus","Mpi minus; m(MeV); Events",60,infl,supl);
TH1F *h_MpiSQ_plus  = new TH1F("h_MpiSQ_plus" ,"Mpi^{2} ; m^{2}(GeV^{2}); Events",60,0,35);
TH1F *h_MpiSQ_minus = new TH1F("h_MpiSQ_minus","Mpi minus; m^{2}(GeV^{2}); Events",60,0,35);
TH1F *h_Mka_plus    = new TH1F("h_Mka_plus" ,"Mka ; m(MeV); Events",60,infl,supl);
TH1F *h_Mka_minus   = new TH1F("h_Mka_minus","Mka minus; m(MeV); Events",60,infl,supl);
TH1F *h_MkaSQ_plus  = new TH1F("h_MkaSQ_plus" ,"Mka^{2} ; m_{KK}^{2}(GeV^{2}); Events",20,0.8,2.8);
TH1F *h_MkaSQ_minus = new TH1F("h_MkaSQ_minus","Mka minus; m_{KK}^{2}(GeV^{2}); Events",20,0.8,2.8);



void func(TLorentzVector Pion, TLorentzVector Kaon1, TLorentzVector Kaon2, int arrayEnt)
{
    double Mpi = (Pion + Kaon1).M();
    double Mka = (Kaon1 + Kaon2).M();
    double PT = (Pion + Kaon1 + Kaon2).Pt();
    double Bmass = (Pion + Kaon1 + Kaon2).M(); h_Bmass_raw->Fill(Bmass,1);
    if(Bmass < B_m + arrayEnt && Bmass > B_m - arrayEnt) {  // B mass window
        if(PT > 1700) {
            if( (Mpi < D0m - D0_thresh || Mpi > D0m + D0_thresh) && (Mka < D0m - D0_thresh || Mka > D0m + D0_thresh) ){ // D0 exclusion
                if(H1_Charge == -1){ //------------ B+ CASE
                    h_Mpi_plus->Fill(Mpi,1);
                    h_Mka_plus->Fill(Mka,1);
                    h_Bplus->Fill(Mka*Mka/1e6,Mpi*Mpi/1e6);
                    h_MpiSQ_plus->Fill(Mpi*Mpi/1e6,1);
                    h_MkaSQ_plus->Fill(Mka*Mka/1e6,1);
                    h_Bmass->Fill(Bmass,1);
                    whole->Fill(Mka*Mka/1e6,Mpi*Mpi/1e6);
                    fprintf(fp, "%f,%f\n", Mka*Mka/1e6,Mpi*Mpi/1e6);
                }
                if(H1_Charge == 1){  //------------ B- CASE
                    h_Mpi_minus->Fill(Mpi,1);
                    h_Mka_minus->Fill(Mka,1);
                    h_Bminus->Fill(Mka*Mka/1e6,Mpi*Mpi/1e6);
                    h_MpiSQ_minus->Fill(Mpi*Mpi/1e6,1);
                    h_MkaSQ_minus->Fill(Mka*Mka/1e6,1);
                    h_Bmass->Fill(Bmass,1);
                    whole->Fill(Mka*Mka/1e6,Mpi*Mpi/1e6);
                    fprintf(fp, "%f,%f\n", Mka*Mka/1e6,Mpi*Mpi/1e6);
                }
            }
        }
    }
}


void cKKp_gen(){

  B2HHH->Add("data/B2HHH_MagnetUp.root");                            B2HHH->Add("data/B2HHH_MagnetDown.root");
  B2HHH->SetBranchAddress("B_FlightDistance", &B_FlightDistance);    B2HHH->SetBranchAddress("B_VertexChi2", &B_VertexChi2);
  B2HHH->SetBranchAddress("H1_isMuon", &H1_isMuon);     B2HHH->SetBranchAddress("H2_isMuon", &H2_isMuon);     B2HHH->SetBranchAddress("H3_isMuon", &H3_isMuon);
  B2HHH->SetBranchAddress("H1_Charge", &H1_Charge);     B2HHH->SetBranchAddress("H2_Charge", &H2_Charge);     B2HHH->SetBranchAddress("H3_Charge", &H3_Charge);
  B2HHH->SetBranchAddress("H1_IPChi2", &H1_IPChi2);     B2HHH->SetBranchAddress("H2_IPChi2", &H2_IPChi2);     B2HHH->SetBranchAddress("H3_IPChi2", &H3_IPChi2);
  B2HHH->SetBranchAddress("H1_PX", &H1_PX);             B2HHH->SetBranchAddress("H2_PX", &H2_PX);             B2HHH->SetBranchAddress("H3_PX", &H3_PX);
  B2HHH->SetBranchAddress("H1_PY", &H1_PY);             B2HHH->SetBranchAddress("H2_PY", &H2_PY);             B2HHH->SetBranchAddress("H3_PY", &H3_PY);
  B2HHH->SetBranchAddress("H1_PZ", &H1_PZ);             B2HHH->SetBranchAddress("H2_PZ", &H2_PZ);             B2HHH->SetBranchAddress("H3_PZ", &H3_PZ);
  B2HHH->SetBranchAddress("H1_ProbPi", &H1_ProbPi);     B2HHH->SetBranchAddress("H2_ProbPi", &H2_ProbPi);     B2HHH->SetBranchAddress("H3_ProbPi", &H3_ProbPi);
  B2HHH->SetBranchAddress("H1_ProbK", &H1_ProbK);       B2HHH->SetBranchAddress("H2_ProbK", &H2_ProbK);       B2HHH->SetBranchAddress("H3_ProbK", &H3_ProbK);
  nentries = (Int_t)B2HHH->GetEntries();
  fp = fopen("KKP.csv","w"); // Creating a CSV to write data

  int B_window[]={40,45,50,55,60};

  for(int ii = 0; ii< size(B_window); ++ii){

    for (i = 0; i < nentries; i++)
    {
        nbytes = B2HHH->GetEntry(i);

        if(H1_isMuon == 1 || H2_isMuon == 1 || H3_isMuon == 1) continue;      // NO MUONS
        if(H1_ProbPi > H1_pi || H1_ProbK < H1_ka) continue;                   // H1 MUST BE A KAON

        TLorentzVector Pion   = TLorentzVector();
        TLorentzVector Kaon1  = TLorentzVector();   Kaon1.SetXYZM(H1_PX, H1_PY, H1_PZ,K_m);
        TLorentzVector Kaon2  = TLorentzVector();

        if(H2_ProbPi > above_this && H2_ProbK < below_this){      //----- H2 IS THE PION
            if(H3_ProbK  > above_this && H3_ProbPi < below_this) {
                Pion.SetXYZM(H2_PX, H2_PY, H2_PZ,Pi_m);
                Kaon2.SetXYZM(H3_PX, H3_PY, H3_PZ,K_m);
                func(Pion,Kaon1,Kaon2,B_window[ii]);
            }
        }
        else if(H3_ProbPi > above_this && H3_ProbK < below_this){ //----- H3 IS THE PION
            if(H2_ProbK  > above_this && H2_ProbPi < below_this) {
                Pion.SetXYZM(H3_PX, H3_PY, H3_PZ,Pi_m);
                Kaon2.SetXYZM(H2_PX, H2_PY, H2_PZ,K_m);
                func(Pion,Kaon1,Kaon2,B_window[ii]);
            }
        }
    }

    // -------------------------------------------------------------------------- PLOT
    auto canvas = new TCanvas("canvas","canvas",1000,700);

    h_MkaSQ_plus->SetTitle(Form("Probs: %.2f-%.2f | B window: #pm%i MeV", above_this, below_this, B_window[ii]));
    h_MkaSQ_plus->SetLineWidth(2);  h_MkaSQ_plus->Draw("hist e");
    h_MkaSQ_minus->SetLineWidth(2); h_MkaSQ_minus->Draw("samehist e");
    h_MkaSQ_minus->SetLineColor(1);
    TLegend *leg = new TLegend(0.7,0.6,0.9,0.75);
    leg->AddEntry(h_MkaSQ_plus, "B+", "l"); leg->AddEntry(h_MkaSQ_minus, "B-", "l"); leg->Draw();
    // canvas->Print(Form("KKP_plots/KKP_focusAsym_Probs%.1f-%.1f_B%i.png", above_this, below_this, B_window[ii]),"Title: ");
    canvas->SaveAs(Form("KKP_plots/KKP_focusAsym_Probs%.2f-%.2f_B%i.png", above_this, below_this, B_window[ii]));


    // -------------------------------------------------------------------------- RESETTING HISTOGRAMS
    whole->Reset("ICESM");         h_PT->Reset("ICESM");
    h_Bmass_raw->Reset("ICESM");   h_Bmass->Reset("ICESM");
    h_Bplus->Reset("ICESM");       h_Bminus->Reset("ICESM");
    h_Mpi_plus->Reset("ICESM");    h_Mpi_minus->Reset("ICESM");
    h_MpiSQ_plus->Reset("ICESM");  h_MpiSQ_minus->Reset("ICESM");
    h_Mka_plus->Reset("ICESM");    h_Mka_minus->Reset("ICESM");
    h_MkaSQ_plus->Reset("ICESM");  h_MkaSQ_minus->Reset("ICESM");

  }

  fclose(fp);   fp = 0;   // Closing the CSV
  vecBin.clear();         // Clearing vector
}
