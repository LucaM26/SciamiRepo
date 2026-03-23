#ifndef RATE_H
#define RATE_H

#include <TTree.h>
#include <TGraphErrors.h>

struct Container {

    std::vector<double> y, y_err;

    void Fill(double y_val, double yerr_val);

};

TGraphErrors* rateGraph (std::vector<double> x, std::vector<double> y,std::vector<double> x_err, std::vector<double> y_err);

void rateMain(TTree *t, int n_ch, double Dt);

#endif