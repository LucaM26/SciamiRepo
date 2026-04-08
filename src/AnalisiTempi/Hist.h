#ifndef HIST_H
#define HIST_H

#include <string>
#include <TH1F.h>
#include <TTree.h>

void hist(TH1F * h, const std::string &fname);

void histMain(TTree *t, int n_ch);

void histTriple(TTree* t, int n_t);


#endif