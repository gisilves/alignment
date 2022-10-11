#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include "TFile.h"
#include "TTree.h"
#include "TH3F.h"

#include "Functions.h"

using namespace std;

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        cout << "Usage: ./Viewer input_file #stations" << endl;
        exit(-1);
    }

    int msd_stations = atoi(argv[2]);
    TH3F *viewer = new TH3F("3d_view", "3D View", 200, -5, 5, 200, -5, 5, 200, -20, 20);
    // variable declaration

    TFile *in_file = new TFile(argv[1], "R");
    if (in_file->IsZombie())
    {
        cout << "\n\nError reading input root filen\n";
        exit(-1);
    }
    gDirectory->cd();

    TTree *myTree = (TTree *)gDirectory->Get("tree");

    Double_t tmp_val[msd_stations][3], l_distances[msd_stations];
    // myTree->SetBranchAddress("Xmeas", tmp_val);
    if (myTree->SetBranchAddress("MSDmeas", tmp_val) < 0)
    {
        cout << "ERROR: MSDmeas branch is not present, checking for Xmeas branch" << endl;
        if (myTree->SetBranchAddress("Xmeas", tmp_val) < 0)
        {
            cout << "ERROR: neither Xmeas nor MSDmeas branch is present" << endl;
            return -1;
        }
        else
        {
            cout << "Found Xmeas branch" << endl;
            myTree->SetBranchAddress("Xmeas", tmp_val);
        }
    }
    else
    {
        cout << "Found MSDmeas branch" << endl;
        myTree->SetBranchAddress("MSDmeas", tmp_val);
    }

    for (int idx = 0; idx < myTree->GetEntries(); idx++)
    {
        myTree->GetEntry(idx);
        for (int i = 0; i < msd_stations; i++)
        {
            //  cout << tmp_val[i][0] << " " << tmp_val[i][1] << " " << tmp_val[i][2] << endl;
            viewer->Fill(tmp_val[i][0], tmp_val[i][1], tmp_val[i][2]);
        }
    }

    viewer->Draw();
    viewer->SaveAs("viewer.root");
    return 0;
}
