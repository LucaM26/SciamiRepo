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

struct RateData{

    std::vector<double> T_08, T_06, T_04, 
                        D12_08, D12_06, D13_04,
                        D23_08, D23_06, D23_04
                        D13_08, D23_06, D23_04;

    void RateData::Fill(double &Val1,double &Val2, double &Val3, 
                        double &Val4, double &Val5, double &Val6,
                        double &Val7, double &Val8, double &Val9,
                        double &Val10, double &Val11, double &Val12){

                            T_08.push_back(Val1);
                            D12_08.push_pack(Val2)
                            D23_08.push_back(Val3);
                            D13_08.push_back(Val4);
                            T_06.push_back(Val5);
                            D12_06.push_back(Val6);
                            D23_06.push_back(Val7);
                            D13_06.push_back(Val8);
                            T_04.push_back(Val9);
                            D12_04.push_back(Val10);
                            D23_04.push_back(Val11);
                            D13_04.push_back(Val12);

                        }
    struct Iterator {
        RateData* rd;
        size_t index;    // indice di Fill
        size_t sub;      // quale vettore (0..11)

        Iterator(RateData* r, size_t i, size_t s) : rd(r), index(i), sub(s) {}

        double operator*() const {
            switch(sub){
                case 0: return rd->T_08[index];
                case 1: return rd->D12_08[index];
                case 2: return rd->D23_08[index];
                case 3: return rd->D13_08[index];
                case 4: return rd->T_06[index];
                case 5: return rd->D12_06[index];
                case 6: return rd->D23_06[index];
                case 7: return rd->D13_06[index];
                case 8: return rd->T_04[index];
                case 9: return rd->D12_04[index];
                case 10: return rd->D23_04[index];
                case 11: return rd->D13_04[index];
            }
            return 0.0;
        }

        Iterator& operator++() {
            if(sub < 11) sub++;
            else { sub=0; index++; }
            return *this;
        }

        bool operator!=(const Iterator& other) const {
            return index != other.index || sub != other.sub;
        }
    };

    Iterator begin() { return Iterator(this,0,0); }
    Iterator end() {
        if(T_08.empty()) return Iterator(this,0,0);
        return Iterator(this,T_08.size(),0);
    }
};



std::tuple<TGraphErrors*, float> rateGraph (const std::vector<double> &x, std::vector<double> &y, const std::vector<double> &x_err, const std::vector<double> &y_err){

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

std::tuple<std::vector<double>, std::vector<double>, RateData, RateData> rateFill(TTree *t, int n_ch, double Dt){

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

    RateData RD;
    RateData RD_err;

    timings.reserve(N / 10);
    timings_err.reserve(N / 10);

    ULong64_t k = 0;
    int n_pmt = n_ch / 3;

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
            timings_err.push_back(5.36870912/3600);

            RD.Fill(counter[0]/Dt, counter[1]/Dt, counter[2]/Dt, counter[3]/Dt,
                    counter[4]/Dt, counter[5]/Dt, counter[6]/Dt, counter[7]/Dt,
                    counter[8]/Dt, counter[9]/Dt, counter[10]/Dt, counter[11]/Dt);

            RD_err.Fill(sqrt(counter[0])/Dt, sqrt(counter[1])/Dt, sqrt(counter[2])/Dt, sqrt(counter[3])/Dt,
                        sqrt(counter[4])/Dt, sqrt(counter[5])/Dt, sqrt(counter[6])/Dt, sqrt(counter[7])/Dt,
                        sqrt(counter[8])/Dt, sqrt(counter[9])/Dt, sqrt(counter[10])/Dt, sqrt(counter[11])/Dt);

            for (int l = 0; l < n_ch; l++){ counter[l] = 0.; }
        }

        k++;
    }

    std::tuple<std::vector<double>, std::vector<double> ,RateData, RateData> results = std::make_tuple(timings, timings_err, RD, RD_err);

    return results;

}

void RateMain(TTree *t, int n_ch, double Dt){

    std::tuple<std::vector<double>, std::vector<double>, RateData, RateData> res = rateFill(t, n_ch, Dt);

    for (auto v : RD){

        

    }

}





