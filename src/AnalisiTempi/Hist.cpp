#include <cmath>
#include <vector>
#include <iostream>
#include <string>

#include <TROOT.h>
#include <TTree.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TPaveText.h>
#include <TLine.h>
#include <TKey.h>
#include <TObject.h>
#include "TFitResultPtr.h"
#include "TFitResult.h"

#include "Hist.h"


void histFit(TH1F *h, const std::string &fname) {
    int nBins = h->GetNbinsX();
    double Max_est = h->GetMaximum();
    double xmin = h->GetXaxis()->GetXmin();
    double xmax = h->GetXaxis()->GetXmax();
    double avg = h->GetMean();

    TF1 *f = new TF1("f", "[0] * exp(-[1]*x) + [2]", xmin, xmax);
    f->SetParameters(Max_est, 1/avg, 0);
    f->SetNpx(1000);
    TFitResultPtr fit = h->Fit(f, "ILS Q"); 

    // Calcolo chi²
    double chi2 = f->GetChisquare();
    int ndf = f->GetNDF();

    //calcolo residui

    TGraph *gResidui = new TGraph();
    int point = 0;
    for (int i = 1; i <= nBins; i++) {
        double err = h->GetBinError(i);
        if (err > 0) {
            double x = h->GetXaxis()->GetBinCenter(i);
            double res = (h->GetBinContent(i) - f->Eval(x)) / err;
            gResidui->SetPoint(point++, x, res);
        }
    }

    // Disegno
    TCanvas *C = new TCanvas("c", "c", 800, 900);

    // Dividi in due pad verticali (asimmetrici)
    TPad *pad1 = new TPad("pad1", "fit",      0, 0.4, 1, 1.0);
    TPad *pad2 = new TPad("pad2", "residui",  0, 0.0, 1, 0.3);

    pad1->SetBottomMargin(0.02);
    pad2->SetTopMargin(0.01);
    pad2->SetBottomMargin(0.3);

    pad1->Draw();
    pad2->Draw();

    // Disegna il fit nel pad grande
    pad1->cd();
    h->SetStats(0);
    h->GetXaxis()->SetTitle("Differenze temporali [s]");
    h->GetYaxis()->SetTitle("Occorrenze [puro]");
    h->GetXaxis()->SetLabelSize(0);  
    h->GetYaxis()->SetLabelSize(0.075);
    h->GetYaxis()->SetTitleSize(0.075);
    h->GetXaxis()->SetTitle("");     
    h->GetXaxis()->SetTickLength(0);
    h->Draw("P E1");

    // Creazione della casella
    double x1 = 0.65, y1 = 0.55;
    double x2 = 0.9,  y2 = 0.85;
    TPaveText *pt = new TPaveText(x1, y1, x2, y2, "NDC");
    pt->SetFillColorAlpha(kWhite, 0.8); 
    pt->SetLineColor(kBlack);
    pt->SetTextFont(42);
    pt->SetTextSize(0.06);
    pt->SetTextAlign(12); 

    pt->AddText(Form("Parametri di fit:"));
    pt->AddText(Form("A : %.0f +/- %.0f", f->GetParameter(0), f->GetParError(0)));
    pt->AddText(Form("R : %.3f +/- %.3f Hz", f->GetParameter(1), f->GetParError(1)));
    pt->AddText(Form("B : %.1f +/- %.1f", f->GetParameter(2), f->GetParError(2)));
    pt->AddText(Form("D/NDF = %.2f", fit->MinFcnValue()/f->GetNDF()));

    pt->Draw(); 
    f->Draw("same");

    // Disegna i residui nel pad piccolo
    pad2->cd();
    gResidui->SetMarkerStyle(20);
    gResidui->SetMarkerSize(0.5);
    gResidui->GetXaxis()->SetTitle("Differenze temporali [s]");
    gResidui->GetYaxis()->SetTitle("Residui normalizzati [puro]");
    gResidui->GetXaxis()->SetLabelSize(0.1);
    gResidui->GetXaxis()->SetTitleSize(0.1);
    gResidui->GetYaxis()->SetLabelSize(0.1);
    gResidui->GetYaxis()->SetTitleSize(0.1);
    gResidui->Draw("AP");

    TLine *line = new TLine(xmin, 0, xmax, 0);
    line->SetLineColor(kRed);
    line->SetLineStyle(2);
    line->Draw();
    C->Update();

    // Salvataggio
    TFile file(fname.c_str(), "RECREATE");
    C->Write();
    file.Close();
    delete C;
}


void histMain(TTree* t, int n_ch) {

    gROOT->SetBatch(kTRUE);

    // Variabili originali
    ULong64_t ch;
    ULong64_t time;

    t->SetBranchAddress("Channels", &ch);
    t->SetBranchAddress("Time", &time);

    ULong64_t N = t->GetEntries();

    std::vector<double> t_max = {3.5, 2.5, 4, 3, 5, 4.5, 5, 5, 1.5, 1, 1, 1};
    std::vector<std::string> titles = {"Triple", "Doppie 1&2", "Doppie 2&3", "Doppie 1&3"};

    std::vector<TH1F*> hist(n_ch);
    for (int i = 0; i < n_ch; i++)
        hist[i] = new TH1F(Form("h_ch%d", i), titles[i % 4].c_str(), 100, 0, t_max[i]);

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


    std::vector<std::string> tel = {"T08", "T06", "T2004"};
    std::vector<double> R_est = {3.5, 6.3, 3.9, 10, 1.1, 2.4, 4.6, 1.2, 9.9, 17.6, 18.8, 10.3} ;
    for (int j = 0; j < n_ch; j += 4) {

        TFile *f = new TFile(Form("FileRoot/Exp%d.root", j), "RECREATE");
        TCanvas* C = new TCanvas(Form("C_%d", j), "", 800, 800);

        std::vector<TGraphErrors*> graphs;   // per salvare i grafici
        std::vector<TF1*> fits;              // per salvare i fit

        for (int k = 0; k <= 3; k++) {

            TH1F* h = hist[j + k];
            int nBins = h->GetNbinsX();

            // ---- Converti istogramma in TGraphErrors ----
            TGraphErrors* g = new TGraphErrors(nBins);
            for (int i = 1; i <= nBins; i++) {
                double x = h->GetBinCenter(i);
                double y = h->GetBinContent(i);
                double err = h->GetBinError(i);
                g->SetPoint(i - 1, x, y);
                g->SetPointError(i - 1, 0, err);
            }

            g->SetMarkerStyle(20); 
            g->SetMarkerSize(0.8);
            g->SetMarkerColor(k + 1); 
            g->SetLineColor(k + 1); 

            // ---- Fit manual senza redraw automatico ----
            double Max_est = h->GetMaximum();
            TF1* exp = new TF1(Form("f_%d", k), "[0] * exp(-[1]*x) + [2]", 0, t_max[j+k]);
            double b_est = 0;
            for (int i = 0; i < 20; i++){b_est += h->GetBinContent(nBins-i);}
            b_est /= 20;
            exp->SetParameters(Max_est, R_est[j+k], b_est);
            exp->SetLineColor(k+1);
            exp->SetLineWidth(2);
            TFitResultPtr fitResult = h->Fit(exp, "ILS"); 

            // devianza normalizzata
            double devianza_norm = 2 * fitResult->MinFcnValue()/exp->GetNDF();

            std::cout << "Devianza normalizzata: " << devianza_norm << std::endl;

            graphs.push_back(g);
            fits.push_back(exp);
        }

        
        for (size_t k = 0; k < graphs.size(); k++) {
            if (k == 0) {
                graphs[k]->GetXaxis()->SetTitle("Differenze temporali [s]");
                graphs[k]->GetYaxis()->SetTitle("Occorrenze [puro]");
                graphs[k]->SetTitle(Form("Differenze temporali %s", tel[j/4].c_str()));
                graphs[k]->Draw("AP E");  // primo grafico imposta assi
                fits[k]->Draw("SAME");
            } else {
                graphs[k]->Draw("P E SAME");  // sovrapposto
                fits[k]->Draw("SAME");
            }
        }

        // ---- Legenda ----
        TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        for (size_t k = 0; k < graphs.size(); k++)
            leg->AddEntry(graphs[k], titles[(j + k) % 4].c_str(), "lep");
        leg->Draw();

        C->Write();
        f->Close();

        // ---- Cleanup ----
        for (auto g : graphs) delete g;
        for (auto fitt : fits) delete fitt;
        delete C;
        delete f;
    }

    gROOT->SetBatch(kFALSE);
}



void histTriple(TTree* t, int n_t) {

    gROOT->SetBatch(kTRUE);

    // Variabili originali
    ULong64_t ch;
    ULong64_t time;

    t->SetBranchAddress("Channels", &ch);
    t->SetBranchAddress("Time", &time);

    ULong64_t N = t->GetEntries();

    std::vector<TH1F*> hist(n_t);
    std::vector<std::string> names = {"08-06", "06-2004", "08-2004"};
    for (int i = 0; i < n_t; i++){ hist[i] = new TH1F(Form("h_%d", i), Form("Differenze temporali %s;time [s]", names[i].c_str()), 100, 0, 5);}

    std::vector<std::pair<int, int>> bits = {{0, 4}, 
                                   {4, 8},
                                   {8, 0}};

    for(ULong64_t k = 0; k < N; k++){

        t->GetEntry(k);

        if (ch == 2147483648) continue;

        for (long unsigned int i = 0; i < bits.size(); i++){

            bool rel_first = (ch >> bits[i].first) & 1;

            bool rel_second = (ch >> bits[i].second) & 1;

            if (rel_first || rel_second){

                ULong64_t start_time = time;
                ULong64_t stop_time = 0;


                for (ULong64_t h = k+1; h < N; h++){

                    t->GetEntry(h);

                    bool rel_first_conf = (ch >> bits[i].first) & 1;

                    bool rel_second_conf = (ch >> bits[i].second) & 1;

                    if ((rel_first && rel_second_conf) || (rel_second && rel_first_conf)){
                        stop_time = time;
                        break;
                    }
                }

                if (stop_time > start_time){hist[i]->Fill((stop_time - start_time) * 5e-9);}
            }

        }

    }

    for (int j = 0; j < n_t; j++){

        histFit(hist[j], Form("FileRoot/h%d.root", j));}

        

    gROOT->SetBatch(kFALSE);

}
