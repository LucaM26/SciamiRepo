#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <tuple>

#include <TROOT.h>
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

std::tuple<TGraphErrors*, float> rateGraph (const std::vector<double> &x,
                         const std::vector<double> &y,
                         const std::vector<double> &x_err,
                         const std::vector<double> &y_err){

    TGraphErrors *g = new TGraphErrors(x.size(), x.data(), y.data(), x_err.data(), y_err.data());

    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.1);
    g->SetTitle("Rate vs Tempo");

    g->GetXaxis()->SetTitle("Tempo [h]");
    g->GetYaxis()->SetTitle("Rate [Hz m^{-2} s^{-1}]");

    double xmax = *std::max_element(x.begin(), x.end());

    TF1 *f = new TF1("f", "[0]", 0, xmax);
    f->SetParameter(0, 5);
    g->Fit(f, "Q");

    float p = f->GetParameter(0);
    std::vector<float> diff;
    for (long unsigned int h = 0; h < y.size(); h++){
        diff.push_back(std::abs(100*(y[h]-p)/p));
    }

    float max = *std::max_element(diff.begin(), diff.end());

    std::tuple<TGraphErrors*, float> res;
    res = make_tuple(g, max);

    return res;
}

void rateMain(TTree *t, int n_ch, double Dt){

    gROOT->SetBatch(kTRUE);

    ULong64_t ch;
    ULong64_t time;

    t->SetBranchAddress("Channels", &ch);
    t->SetBranchAddress("Time", &time);

    ULong64_t N = t->GetEntries();

    double el = 0.;
    double tot = 0.;

    std::vector<double> timings;
    std::vector<double> timings_err;

    std::vector<double> counter(n_ch, 0.);
    std::vector<Container> data(n_ch);
    std::vector<Container> eff_data(n_ch);

    timings.reserve(N / 10);
    timings_err.reserve(N / 10);

    ULong64_t k = 0;
    int step = n_ch / 3;

    while (k < N) {

        t->GetEntry(k);

        if (ch == 2147483648) {
            el += 5.36870912;
            tot += 5.36870912/3600;
        }

        if (ch != 2147483648){

            ULong64_t mask = ch;

            while (mask){
                int j = __builtin_ctzll(mask);
                if (j < n_ch) counter[j] += 1;
                mask &= (mask - 1);
            }
        }

        if (el >= Dt){

            el = 0.;

            timings.push_back(tot);
            timings_err.push_back(5./3600);

            for (int l = 0; l < n_ch; l++){
                double rate = counter[l]/Dt;
                data[l].Fill(rate, std::sqrt(counter[l])/Dt);
                counter[l] = 0.;
            }
        }

        k++;
    }

    /*

    std::vector<std::vector<double>> C{
        {0.99, 0.98, 0.84},
        {1, 1, 1},
        {1, 1, 1}
    };

    std::vector<std::vector<double>> A = {
        {0.97, 1, 0.92},
        {0.36, 1, 0.42},
        {0.61, 1, 0.61}
    };

    std::vector<float> S = {0.1218, 0.08, 0.2};
    std::vector<float> O = {0.776, 0.171, 0.495};

    for (int m = 0; m < n_ch; m += 4){

        auto& T   = data[m].y;
        auto& D12 = data[m+1].y;
        auto& D23 = data[m+2].y;
        auto& D31 = data[m+3].y;

        const std::vector<double>& A_v = A[(m / 4) % 3];
        const std::vector<double>& C_v = C[(m / 4) % 3];

        for (size_t e = 0; e < T.size(); e++){

            if (T[e] > 0 && D12[e] > 0 && D23[e] > 0 && D31[e] > 0){

                double eff_T  = (T[e]/(D12[e] * C_v[2] * A_v[0])) * (T[e]/(D23[e] * C_v[0] * A_v[1])) * (T[e]/(D31[e] * C_v[1] * A_v[2]));
                double eff_12 = (T[e]/(D23[e] * C_v[0] * A_v[0])) * (T[e]/(D31[e] * C_v[1] * A_v[1]));
                double eff_23 = (T[e]/(D12[e] * C_v[2] * A_v[2])) * (T[e]/(D31[e] * C_v[1] * A_v[1]));
                double eff_31 = (T[e]/(D12[e] * C_v[2] * A_v[2])) * (T[e]/(D23[e] * C_v[0] * A_v[0]));

                eff_data[m].Fill(eff_T, 0.0);
                eff_data[m+1].Fill(eff_12, 0.0);
                eff_data[m+2].Fill(eff_23, 0.0);
                eff_data[m+3].Fill(eff_31, 0.0);

                data[m].y[e]   /= eff_T * O [(m / 4) % 3] * S[(m / 4) % 3];
                data[m+1].y[e] /= eff_12 * O [(m / 4) % 3] * S[(m / 4) % 3];
                data[m+2].y[e] /= eff_23 * O [(m / 4) % 3] * S[(m / 4) % 3];
                data[m+3].y[e] /= eff_31 * O [(m / 4) % 3] * S[(m / 4) % 3];
            }
        }
    }
    */

    int n_files = n_ch / step;
    std::vector<std::string> labels = {"1&2&3", "1&2", "2&3", "1&3"};

    for (int n = 0; n < n_files; n++) {

        TCanvas* c = new TCanvas(Form("c_setup_%d", n), "Rate vs Tempo", 900, 700);
        TLegend* legend = new TLegend(0.65, 0.7, 0.88, 0.88);
        legend->SetHeader("Tipo di evento", "C");

        double ymin = 1e9, ymax = -1e9;
        for (int ch_index = n*step; ch_index < (n+1)*step; ch_index++) {
            for (size_t i = 0; i < data[ch_index].y.size(); i++) {
                if (data[ch_index].y[i] > 0) {
                    ymin = std::min(ymin, data[ch_index].y[i]);
                    ymax = std::max(ymax, data[ch_index].y[i]);
                }
            }
        }

        ymin = std::max(0.0, ymin * 0.9);
        ymax = ymax * 1.1;

        bool first = true;

        for (int ch_index = n*step; ch_index < (n+1)*step; ch_index++) {

            std::tuple<TGraphErrors*, float> res = rateGraph(timings, data[ch_index].y, timings_err, data[ch_index].y_err);

            TGraphErrors* g = std::get<0>(res);

            g->SetMarkerColor((ch_index % step) + 1);
            g->SetLineColor((ch_index % step) + 1);

            g->SetMinimum(ymin);
            g->SetMaximum(ymax);

            if (first) {
                g->Draw("AP");
                first = false;
            } else {
                g->Draw("P SAME");
            }

            TF1* f = g->GetFunction("f");
            if (f) {
                double p = f->GetParameter(0);
                double ep = f->GetParError(0);
                legend->AddEntry(g, Form("%s : %.2f #pm %.2f Hz m^{-2} str^{-1} ", labels[ch_index % step].c_str(), p, ep), "p");
            } else {
                legend->AddEntry(g, labels[ch_index % step].c_str(), "p");
            }
        }

        legend->Draw();
        c->Update();

        std::string name = "FileRoot/Rates_setup" + std::to_string(n+1) + ".root";
        TFile f(name.c_str(), "RECREATE");
        c->Write();
        f.Close();

        delete c;
    }


}