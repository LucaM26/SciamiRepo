#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

#include <TGraph.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TAxis.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TLegend.h>

#include "Rate.h"


void Container::Fill(double y_val, double yerr_val){

        y.push_back(y_val);
        y_err.push_back(yerr_val);

}

TGraphErrors* rateGraph (std::vector<double> x, std::vector<double> y,std::vector<double> x_err, std::vector<double> y_err){

    TGraphErrors *g = new TGraphErrors(x.size(), x.data(), y.data(), x_err.data(), y_err.data());
    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.2);
    g->SetMinimum(0);
    g->SetTitle("Rate vs Tempo");
    g->GetXaxis()->SetTitle("Tempo [h]");
    g->GetYaxis()->SetTitle("Rate [Hz]");

    TF1 *f = new TF1("f", "[0]", 0, *std::max_element(x.begin(), x.end()));
    f->SetParameters(5);
    g->Fit("f", "S");
    return g;

};

void rateMain(TTree *t, int n_ch, double Dt){

    //Allestimento del Tree

    ULong64_t ch;
    ULong64_t time;

    t->SetBranchAddress("Channels", &ch);
    t->SetBranchAddress("Time", &time);

    //Numero di eventi

    ULong64_t N = t->GetEntries();

    //Variabili temporali e di conteggio

    double el = 0.;
    double tot = 0.;
    std::vector<double> timings;
    std::vector<double> timings_err;
    std::vector<double> counter(n_ch,0.);

    //Container

    std::vector<Container> data(n_ch);

    //Indice del while

    ULong64_t k = 0;

    //Passo del loop, 4 nella versione definitiva

    int step = n_ch/3;

    while (k < N) {

        t->GetEntry(k);

        if (ch == 2147483648) {

            el += 5;
            tot += 5./3600; 

        }

        for (int j = 0; j < n_ch; j++){

            if (((ch >> j) & 1) && (ch != 2147483648)){

                counter[j] += 1; 
            }
        }

        if (el >= Dt){

            el = 0.;

            timings.push_back(tot);
            timings_err.push_back(5./3600);

            for (int l = 0; l < n_ch; l++){

                data[l].Fill(counter[l]/Dt,std::sqrt(counter[l])/Dt);

                counter[l] = 0.;
            }
        }

        k++;

    }

    //Da inserire la correzione dell'efficienza

    int n_files = n_ch / step;

    std::vector<double> maxes = {15, 2, 10};

    for (int n = 0; n < n_files; n ++){

        TCanvas* c = new TCanvas(Form("c_setup_%d", n), "Rate vs Tempo", 800, 600);

        TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend->SetHeader("Tipo di evento","C");

        bool first = true;

        for (int ch_index = n*step; ch_index < (n+1)*step; ch_index++){

            TGraphErrors* g = rateGraph(timings, data[ch_index].y, timings_err, data[ch_index].y_err);

            g->SetMarkerColor(ch_index+1);   
            g->SetLineColor(ch_index+1);
            g->SetMaximum(maxes[n]);

            g->Draw(first ? "AP" : "P SAME");
            first = false;

            std::string label = (ch_index % step == 0 ? "Triple" : "Doppie");
            legend->AddEntry(g, label.c_str(), "p");
        }

        legend->Draw("SAME");
        c->Update();
        
        std::string name = "Rates_setup" + std::to_string(n+1) + ".root";
        TFile* f = new TFile(name.c_str(), "RECREATE");

        c->Write();
        f->Close();
        delete c;
    }

}




