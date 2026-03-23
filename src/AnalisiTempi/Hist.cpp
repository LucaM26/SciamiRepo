#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <array>
#include <cstring>

#include <TGraph.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TAxis.h>
#include <TF1.h>
#include <TTree.h>
#include <TFile.h>

void histFit(TH1F *h, std::string fname) {

    TF1 *f = new TF1("f", "[0] * exp(-[1]*x) + [2]", 0, 4);

    ULong64_t Max_est= h->GetMaximum();

    f->SetParameters(Max_est, 5, 0);

    h->Fit("f", "ILS");
    
    TCanvas *C = new TCanvas();
    h->Draw();
    C->Update();

    TFile* file = new TFile(fname.c_str(), "RECREATE");

    h->Write();
    file->Close();

    delete C;

}

void histMain(TTree* t, int n_ch){

    //Allestimento del Tree

    ULong64_t ch;
    ULong64_t time;

    t->SetBranchAddress("Channels", &ch);
    t->SetBranchAddress("Time", &time);

    //Numero di eventi

    ULong64_t N = t->GetEntries();

    //Creo gli istogrammi

    std::vector<TH1F*> hist(n_ch);

    for (int i = 0; i < n_ch; i++){

        hist[i] = new TH1F(Form("h_ch%d", i), Form("Differenze temporali ch %d;time", i), 100, 0, 4);
    }

    for (ULong64_t k = 0; k < N; k++){

        t->GetEntry(k);

        if (ch == 2147483648){continue;}

        for (int bit = 0; bit < n_ch; bit++){

            if ((ch >> bit) & 1){

                ULong64_t start = time;

                ULong64_t h = 1;

                while (k + h < N) {

                    t->GetEntry(k+h);

                    if ((ch >> bit) & 1) {break;}

                    h++;
                }

                ULong64_t stop = time;

                hist[bit]->Fill((stop-start)*5e-9);

            }

        }
    }

    for (int j = 0; j < n_ch; j++){

        histFit(hist[j], Form("h_ch%d.root", j));

    }

}