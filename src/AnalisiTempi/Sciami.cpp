#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <array>
#include <cstring>

#include <TGraph.h>
#include <TCanvas.h>
#include <TChain.h>

#include "Rate.h"
#include "Hist.h"

//Crea il file .root con i dati del FIFO

void FileToTree (const char * fname, const std::string &outname ,const char * option = "RECREATE") {

    //Apertura del file con test

    std::ifstream data(fname);

    if (!data) {cout << "Errore di apertura \n" << endl;}

    //Inizializzazione del File di output e del Tree

    TFile OutFile((outname + ".root").c_str(), option);

    TTree Tree("T", "Contiene i dati del DE10-NANO");

    //Variabili, canali scattati e tempi in digit

    ULong64_t ch = 0;

    ULong64_t time = 0;

    //Creazione dei branch

    Tree.Branch("Channels", &ch, "Channel/l");

    Tree.Branch("Time", &time, "Time/l");

    //Riempimento del TTree

    while (data >> ch >> time){

        Tree.Fill() ;

        //Per controllo puoi printare le entries del file, ma è sconsigliato per file lunghi

        /*

        std::cout << ch << "," << time <<endl;

        */

    }

    Tree.Write();
    OutFile.Close();
}

void createChain(const char * f_1, const char * f_2){

    TChain chain("T");
    chain.Add(f_1);
    chain.Add(f_2);

    std::cout << "Numero di eventi nella chain: " << chain.GetEntries() << std::endl;

    TFile fOut("Events.root", "RECREATE");

    TTree *treeOut = chain.CloneTree(-1);

    treeOut->Write();  
    fOut.Close();

}

//Fa partire un'analisi a scelta; le keyword sono "rate" e "hist"

void Sciami(const char * path, const char* kw){

    //Apertura file e ottenimento del TTree

    TFile *f = TFile::Open(path);
    TTree *t = (TTree*)f->Get("T");

    //Stringa di scelta
    if (strcmp(kw, "rate") == 0) {rateMain(t, 12, 1200);}

    if (strcmp(kw, "hist") == 0) {histMain(t, 12);}

    if (strcmp(kw, "hist2") == 0) {histTriple(t, 3);}

    f->Close();
    delete f;

}