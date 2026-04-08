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
    g->SetTitle("Flusso diff. vs Tempo");

    g->GetXaxis()->SetTitle("Tempo [HH:MM]");
    g->GetXaxis()->SetTimeDisplay(1);
    g->GetXaxis()->SetTimeFormat("%H:%M");
    g->GetXaxis()->SetTimeOffset(0,"gmt");
    g->GetXaxis()->SetTickLength(0.03);

    double x_min = *std::min_element(x.begin(), x.end());
    double x_max = *std::max_element(x.begin(), x.end());

    int n_major = std::ceil((x_max - x_min) / 3600. / 4.);
    g->GetXaxis()->SetNdivisions(n_major, kFALSE); 

    g->GetYaxis()->SetTitle("Flusso diff. [Hz m^{-2} sr^{-1}]");

    double xmax = *std::max_element(x.begin(), x.end());

    TF1 *f = new TF1("f", "[0]", 0, xmax);
    f->SetParameter(0, 5);
    g->Fit(f, "Q0");

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
    double start_tot =  TDatime(2026, 3, 25, 10, 38, 56).Convert();
    double tot = 0.;

    std::vector<double> timings;
    std::vector<double> timings_err;

    std::vector<double> counter(n_ch, 0.);
    std::vector<Container> data(n_ch);
    std::vector<Container> eff_data(n_ch-3);

    timings.reserve(N / 10);
    timings_err.reserve(N / 10);

    ULong64_t k = 0;
    int step = n_ch / 3;

    while (k < N) {

        t->GetEntry(k);

        if (ch == 2147483648) {

            el += 5.36870912;
            tot += 5.36870912;

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

            timings.push_back(tot + start_tot);
            timings_err.push_back(Dt/4);

            for (int l = 0; l < n_ch; l++){

                double rate = counter[l]/Dt;
                data[l].Fill(rate, std::sqrt(counter[l])/Dt);

            }

            for (int l = 0; l < 3; l++){ 

                double e1 = (counter[4 * l]/counter[4 * l + 2]);
                double e1_err = sqrt(e1 * (1 - e1)/counter[4 * l + 2]);
                double e2 = counter[4 * l]/counter[4 * l + 3];
                double e2_err = sqrt(e2 * (1 - e2)/counter[4 * l + 3]);
                double e3 = counter[4 * l]/counter[4 * l + 1];
                double e3_err = sqrt(e3 * (1 - e3)/counter[4 * l + 1]);

                eff_data[3*l].Fill(e1, e1_err);
                eff_data[3*l + 1].Fill(e2, e2_err);
                eff_data[3*l + 2].Fill(e3, e3_err);

            }

            for (int l = 0; l < n_ch; l++){

                counter[l] = 0.;

            }
        }

        k++;
    }


    std::vector<std::vector<double>> C{
        {0.99, 0.84, 0.98,},
        {1, 1, 1},
        {1, 1, 1}
    };

    std::vector<std::vector<double>> A = {
        {0.97, 1, 0.92},
        {0.36, 1, 0.42},
        {0.61, 1, 0.61}
    };

    std::vector<float> S = {0.1218, 0.08, 0.2};
    std::vector<float> dS = {0.0006, 0.0004, 0.0006};

    std::vector<float> O = {1.380, 0.173, 0.570};
    std::vector<float> dO = {0.002, 0.001, 0.002};

    std::map<int, int> eff_ind;
    eff_ind[0] = 0;
    eff_ind[4] = 3;
    eff_ind[8] = 6;

    for (int m = 0; m < n_ch; m += 4){

        const std::vector<double> A_v = A[(m / 4) % 3];
        const std::vector<double> C_v = C[(m / 4) % 3];

        for (size_t e = 0; e < (data[m].y).size(); e++){

            if (data[m].y[e] > 0 && data[m+1].y[e] > 0 && data[m+2].y[e] > 0 && data[m+3].y[e] > 0){

                eff_data[eff_ind[m]].y[e] /= (C_v[0] * A_v[0]);
                eff_data[eff_ind[m]].y_err[e] /= (C_v[0] * A_v[0]);
                eff_data[eff_ind[m] + 1].y[e] /= (C_v[1] * A_v[1]);
                eff_data[eff_ind[m] + 1].y_err[e] /= (C_v[1] * A_v[1]);
                eff_data[eff_ind[m] + 2].y[e] /= (C_v[2] * A_v[2]);
                eff_data[eff_ind[m] + 2].y_err[e] /= (C_v[2] * A_v[2]);

                
                double eff_T  = (eff_data[eff_ind[m]].y[e]) * (eff_data[eff_ind[m]+1].y[e]) * (eff_data[eff_ind[m]+2].y[e]);
                double effT_err = sqrt(pow(eff_data[eff_ind[m]].y_err[e]/eff_data[eff_ind[m]].y[e], 2) + pow(eff_data[eff_ind[m] + 1].y_err[e]/eff_data[eff_ind[m] + 1].y[e], 2) + pow(eff_data[eff_ind[m] + 2].y_err[e]/eff_data[eff_ind[m] + 2].y[e], 2)); 
                double eff_12 = (eff_data[eff_ind[m]].y[e]) * (eff_data[eff_ind[m]+1].y[e]);
                double eff12_err = sqrt(pow(eff_data[eff_ind[m]].y_err[e]/eff_data[eff_ind[m]].y[e],2) + pow(eff_data[eff_ind[m]+1].y_err[e]/eff_data[eff_ind[m]+1].y[e],2));
                double eff_23 = (eff_data[eff_ind[m]+2].y[e]) * (eff_data[eff_ind[m]+1].y[e]);
                double eff23_err =sqrt(pow(eff_data[eff_ind[m] + 2].y_err[e]/eff_data[eff_ind[m] + 2].y[e],2) + pow(eff_data[eff_ind[m]+1].y_err[e]/eff_data[eff_ind[m]+1].y[e],2)); 
                double eff_31 = (eff_data[eff_ind[m]+2].y[e]) * (eff_data[eff_ind[m]].y[e]);
                double eff31_err = sqrt(pow(eff_data[eff_ind[m] + 2].y_err[e]/eff_data[eff_ind[m] + 2].y[e],2) + pow(eff_data[eff_ind[m]].y_err[e]/eff_data[eff_ind[m]].y[e],2));

                double relT = sqrt(pow(data[m].y_err[e]/data[m].y[e],2) + pow(dO [(m / 4) % 3]/O [(m / 4) % 3],2) + pow(dS[(m / 4) % 3]/S[(m / 4) % 3],2) + pow(effT_err/eff_T,2));
                double rel12 = sqrt(pow(data[m + 1].y_err[e]/data[m + 1].y[e],2) + pow(dO [(m / 4) % 3]/O [(m / 4) % 3],2) + pow(dS[(m / 4) % 3]/S[(m / 4) % 3],2) + pow(eff12_err/eff_12,2));
                double rel23 = sqrt(pow(data[m + 2].y_err[e]/data[m + 2].y[e],2) + pow(dO [(m / 4) % 3]/O [(m / 4) % 3],2) + pow(dS[(m / 4) % 3]/S[(m / 4) % 3],2) + pow(eff23_err/eff_23,2));
                double rel31 = sqrt(pow(data[m + 3].y_err[e]/data[m + 3].y[e],2) + pow(dO [(m / 4) % 3]/O [(m / 4) % 3],2) + pow(dS[(m / 4) % 3]/S[(m / 4) % 3],2) + pow(eff31_err/eff_31,2));

                data[m].y[e]   /= eff_T * O [(m / 4) % 3] * S[(m / 4) % 3];
                data[m + 1].y[e] /= eff_12 * O [(m / 4) % 3] * S[(m / 4) % 3];
                data[m + 2].y[e] /= eff_23 * O [(m / 4) % 3] * S[(m / 4) % 3];
                data[m + 3].y[e] /= eff_31 * O [(m / 4) % 3] * S[(m / 4) % 3];

                data[m].y_err[e] = data[m].y[e] * relT;
                data[m + 1].y_err[e] = data[m + 1].y[e] * rel12;
                data[m + 2].y_err[e] = data[m + 2].y[e] * rel23;
                data[m + 3].y_err[e] = data[m + 3].y[e] * rel31;

            }
        }
    }


    
    int n_files = n_ch / step;
    std::vector<std::string> labels = {"1&2&3", "1&2", "2&3", "1&3"};

    for (int n = 0; n < n_files; n++) {

        TCanvas* c = new TCanvas(Form("c_setup_%d", n), "Flux vs Tempo", 900, 700);
        TLegend* legend = new TLegend(0.65, 0.7, 0.88, 0.88);
        legend->SetHeader("Evento", "C");

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
            std::cout << std::get<1>(res) << std::endl;

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
                legend->AddEntry(g, Form("%s : %.3f #pm %.3f Hz m^{-2} sr^{-1} ", labels[ch_index % step].c_str(), p, ep), "p");
            } else {
                legend->AddEntry(g, labels[ch_index % step].c_str(), "p");
            }
        }

        legend->Draw();
        c->Update();

        std::string name = "FileRoot/Flux_setup" + std::to_string(n+1) + ".root";
        TFile f(name.c_str(), "RECREATE");
        c->Write();
        f.Close();

        delete c;
    }


    /*
    std::vector<std::string> labels_eff = {"01", "02", "03"};

    for (int n = 0; n < n_files; n++){

        TCanvas* c = new TCanvas(Form("c_setup_%d", n), "Efficienza vs Tempo", 900, 700);
        TLegend* legend = new TLegend(0.65, 0.7, 0.88, 0.88);
        legend->SetHeader("PMT", "C");

        bool first = true;

        for (int ch_index = n*(step-1); ch_index < (n+1)*(step-1); ch_index++) {

            std::tuple<TGraphErrors*, float> res = rateGraph(timings, eff_data[ch_index].y, timings_err, eff_data[ch_index].y_err);

            TGraphErrors* g = std::get<0>(res);

            g->SetMarkerColor((ch_index % (step-1)) + 1);
            g->SetLineColor((ch_index % (step-1)) + 1);

            g->SetMinimum(0);
            g->SetMaximum(1);

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
                legend->AddEntry(g, Form("%s : %.4f #pm %.4f ", labels_eff[ch_index % (step-1)].c_str(), p, ep), "p");
            } else {
                legend->AddEntry(g, labels_eff[ch_index % (step-1)].c_str(), "p");
            }
        }

        legend->Draw();
        c->Update();

        std::string name = "FileRoot/PMT_setup" + std::to_string(n+1) + ".root";
        TFile f(name.c_str(), "RECREATE");
        c->Write();
        f.Close();

        delete c;

        
    }
    */

    gROOT->SetBatch(kFALSE);

}