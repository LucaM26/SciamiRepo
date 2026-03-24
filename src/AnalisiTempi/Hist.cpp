#include <cmath>
#include <vector>
#include <iostream>
#include <string>

#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TPaveText.h>

void histFit(TH1F *h, const std::string &fname) {
    int nBins = h->GetNbinsX();
    double Max_est = h->GetMaximum();
    double xmin = h->GetXaxis()->GetBinCenter(h->GetMinimumBin());
    double xmax = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
    double avg = h->GetMean();

    TF1 *f = new TF1("f", "[0] * exp(-[1]*x) + [2]", xmin, xmax);
    f->SetParameters(Max_est, 1/avg, 0);
    f->SetNpx(1000);
    h->Fit(f, "ILS Q"); 

    // Calcolo chi²
    double chi2 = f->GetChisquare();
    int ndf = f->GetNDF();

    // Disegno
    TCanvas *C = new TCanvas();
    h->SetStats(0);
    h->Draw();

    // Creazione della casella
    double x1 = 0.65, y1 = 0.55;
    double x2 = 0.9,  y2 = 0.85;
    TPaveText *pt = new TPaveText(x1, y1, x2, y2, "NDC");
    pt->SetFillColorAlpha(kWhite, 0.8); // bianco leggermente trasparente
    pt->SetLineColor(kBlack);
    pt->SetTextFont(42);
    pt->SetTextSize(0.03);
    pt->SetTextAlign(12); // left-top

    pt->AddText(Form("Parametri di fit:"));
    pt->AddText(Form("Parameter A : %.2f ± %.2f", f->GetParameter(0), f->GetParError(0)));
    pt->AddText(Form("Parameter R : %.2f ± %.2f Hz", f->GetParameter(1), f->GetParError(1)));
    pt->AddText(Form("Parameter B : %.2f ± %.2f", f->GetParameter(2), f->GetParError(2)));
    pt->AddText(Form("Chi2/NDF = %.2f", f->GetChisquare()/f->GetNDF()));

    pt->Draw(); 
    C->Update();

    // Salvataggio
    TFile file(fname.c_str(), "RECREATE");
    C->Write();
    file.Close();
    delete C;
}

void histMain(TTree* t, int n_ch) {

    // Variabili originali
    ULong64_t ch;
    ULong64_t time;

    t->SetBranchAddress("Channels", &ch);
    t->SetBranchAddress("Time", &time);

    ULong64_t N = t->GetEntries();

    std::vector<TH1F*> hist(n_ch);
    for (int i = 0; i < n_ch; i++)
        hist[i] = new TH1F(Form("h_ch%d", i), Form("Differenze temporali ch %d;time [s]", i), 100, 0, 4);

    // Riempimento istogrammi
    for (ULong64_t k = 0; k < N; k++){
        t->GetEntry(k);

        if (ch == 2147483648) continue;

        for (int bit = 0; bit < n_ch; bit++){
            if ((ch >> bit) & 1){
                ULong64_t start_time = time;
                ULong64_t stop_time = 0;

                for (ULong64_t h = k+1; h < N; h++){
                    t->GetEntry(h);
                    if ((ch >> bit) & 1){
                        stop_time = time;
                        break;
                    }
                }

                if (stop_time > start_time){
                    hist[bit]->Fill((stop_time - start_time) * 5e-9);
                }
            }
        }
    }

    for (int j = 0; j < n_ch; j++)
        histFit(hist[j], Form("h_ch%d.root", j));
}