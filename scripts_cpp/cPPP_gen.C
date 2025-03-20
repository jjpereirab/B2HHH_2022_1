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
double min_prob_pi = 0.75;    double max_prob_ka = 0.25;   // PARTICLE ID EXCLUSION


//--------- HISTOGRAM DEFINITIONS
int dpxbins = 10;  int dpybins = 14;

TH2F *whole = new TH2F("whole","whole sym; Mlow^{2}; Mhigh^{2}",dpxbins,0,25, dpybins, 0, 35);
TH1F *h_Mlowplus_focused  = new TH1F("h_Mlowplus_focused", "m_{low}^{2}; m_{#pi#pi low}^{2} (GeV^{2}); Events",25,0,2.5);
TH1F *h_Mlowminus_focused = new TH1F("h_Mlowminus_focused","m_{low}^{2}; m_{#pi#pi low}^{2} (GeV^{2}); Events",25,0,2.5);

TH1F *h_Mlowplus   = new TH1F("h_Mlowplus" ,"m_{low}^{2}; m_{#pi#pi}^{2} (GeV^{2}); Events",80,0,15);
TH1F *h_Mlowminus  = new TH1F("h_Mlowminus","m_{low}^{2}; m_{#pi#pi}^{2} (GeV^{2}); Events",80,0,15);
TH1F *h_Mhighplus  = new TH1F("h_Mhighplus","m_{high}^{2}; m(GeV); Events",80,0,30);
TH1F *h_Mhighminus = new TH1F("h_Mhighminus","m_{high}^{2}; m(GeV); Events",80,0,30);
TH2F *h_Bplus      = new TH2F("h_Bplus","B^{+}; m12^{2}; m13^{2}",20,0,35, 20, 0, 35);
TH2F *h_Bminus     = new TH2F("h_Bminus","B^{-}; m12^{2}; m13^{2}",20,0,35, 20, 0, 35);
TH2F *h_Bplus_sym  = new TH2F("h_Bplus_sym","B^{+} sym; Mlow^{2}; Mhigh^{2}",dpxbins,0,25, dpybins, 0, 35);
TH2F *h_Bminus_sym = new TH2F("h_Bminus_sym","B^{-} sym; Mlow^{2}; Mhigh^{2}",dpxbins,0,25, dpybins, 0, 35);



void fun(double X, double Y){
    if(H1_Charge == -1) {
        h_Mlowplus->Fill(Y*Y/1e6);   h_Mhighplus->Fill(X*X/1e6);
        h_Bplus_sym->Fill(Y*Y/1e6,X*X/1e6);
        whole->Fill(Y*Y/1e6,X*X/1e6);
        fprintf(fp, "%f,%f\n", Y*Y/1e6,X*X/1e6);
        if(X*X > 15e6) {h_Mlowplus_focused->Fill(Y*Y/1e6);}
    }
    if(H1_Charge == 1) {
        h_Mlowminus->Fill(Y*Y/1e6);   h_Mhighminus->Fill(X*X/1e6);
        h_Bminus_sym->Fill(Y*Y/1e6,X*X/1e6);
        whole->Fill(Y*Y/1e6,X*X/1e6);
        fprintf(fp, "%f,%f\n", Y*Y/1e6,X*X/1e6);
        if(X*X > 15e6) {h_Mlowminus_focused->Fill(Y*Y/1e6);}
    }
}


void cPPP_gen(){

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
  fp = fopen("PPP.csv","w"); // Creating a CSV to write data

  int B_window[]={40,45,50,55,60};

  for(int ii = 0; ii< size(B_window); ++ii){

    for (i = 0; i < nentries; i++)
    {
        nbytes = B2HHH->GetEntry(i);

        if(H1_isMuon == 1 || H2_isMuon == 1 || H3_isMuon == 1) continue;
        if(H1_ProbPi < min_prob_pi || H2_ProbPi < min_prob_pi || H3_ProbPi < min_prob_pi) continue;
        if(H1_ProbK  > max_prob_ka || H2_ProbK  > max_prob_ka || H3_ProbK  > max_prob_ka) continue;

        // DEFINING AND FEEDING TLORENTZ VECTORS
        TLorentzVector Pion1  = TLorentzVector();   Pion1.SetXYZM(H1_PX, H1_PY, H1_PZ,Pi_m);
        TLorentzVector Pion2  = TLorentzVector();   Pion2.SetXYZM(H2_PX, H2_PY, H2_PZ,Pi_m);
        TLorentzVector Pion3  = TLorentzVector();   Pion3.SetXYZM(H3_PX, H3_PY, H3_PZ,Pi_m);

        // WORKING OUT THE RANDOMIZATION
        TRandom2 *rand = new TRandom2(i);          int r = rand->Binomial(1,0.5);
        TLorentzVector *vec_array = new TLorentzVector[2];
        vec_array[0] = Pion2;      vec_array[1] = Pion3;

        // 2-BODY INVARIANT MASSES (NOTE RANDOMNESS)
        double M12 = (Pion1 + vec_array[r]).M();     double M13 = (Pion1 + vec_array[1-r]).M();

        double PT = (Pion1 + Pion2 + Pion3).Pt();
        double Bmass = (Pion1 + Pion2 + Pion3).M();

        if(Bmass < B_m + B_window[ii] && Bmass > B_m - B_window[ii] && PT > 1700) { // B Pt and Mass cuts

            if( (M12 < D0m - D0_thresh || M12 > D0m + D0_thresh) && (M13 < D0m - D0_thresh || M13 > D0m + D0_thresh) ){ // D0 exclusion

                if(M12 > M13) fun(M12,M13);
                else          fun(M13,M12);

                if(H1_Charge == -1) h_Bplus->Fill(M12*M12/1e6,M13*M13/1e6);
                if(H1_Charge ==  1) h_Bminus->Fill(M12*M12/1e6,M13*M13/1e6);

            }
        }
    }

    // -------------------------------------------------------------------------- PLOT
    TCanvas *canvas = new TCanvas("canvas","canvas",1000,700);

    h_Mlowplus_focused->SetTitle(Form("Probs: %.2f-%.2f | B window: #pm%i MeV", min_prob_pi, max_prob_ka, B_window[ii]));
    h_Mlowplus_focused->SetLineWidth(2);   h_Mlowplus_focused->Draw("hist e");
    h_Mlowminus_focused->SetLineWidth(2);  h_Mlowminus_focused->Draw("samehist e");
    h_Mlowminus_focused->SetLineColor(1);

    int ymaxlim = h_Mlowminus_focused->GetBinContent(h_Mlowminus_focused->GetMaximumBin());
    h_Mlowplus_focused->SetMaximum(1.02*ymaxlim);

    TLegend *leg = new TLegend(0.7,0.6,0.9,0.75);
    leg->AddEntry(h_Mlowplus_focused, "B+", "l"); leg->AddEntry(h_Mlowminus_focused, "B-", "l"); leg->Draw();
    canvas->SaveAs(Form("PPP_plots/PPP_focusAsym_Probs%.2f-%.2f_B%i.png", min_prob_pi, max_prob_ka, B_window[ii]));


    // -------------------------------------------------------------------------- RESETTING HISTOGRAMS
    whole->Reset("ICESM");
    h_Mlowplus_focused->Reset("ICESM");  h_Mlowminus_focused->Reset("ICESM");
    h_Mlowplus->Reset("ICESM");          h_Mlowminus->Reset("ICESM");
    h_Mhighplus->Reset("ICESM");         h_Mhighminus->Reset("ICESM");
    h_Bplus->Reset("ICESM");             h_Bminus->Reset("ICESM");
    h_Bplus_sym->Reset("ICESM");         h_Bminus_sym->Reset("ICESM");

  }

  fclose(fp);   fp = 0;   // Closing the CSV
  vecBin.clear();         // Clearing vector
}
