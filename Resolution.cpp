#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>

#include "TFile.h"
#include "TTree.h"

#include "Functions.h"

using namespace std;

int main(int argc, char **argv)
{
  if (argc < 4)
  {
    cout << "Usage: ./Resolution input_file #stations #verbose" << endl;
    exit(-1);
  }

  int msd_stations = atoi(argv[2]);
  bool verbose = atoi(argv[3]);

  Double_t correct_position[msd_stations - 1][2];
  Double_t rotation_angle[msd_stations - 1][2];

  TFile *in_file = new TFile(argv[1], "R");
  if (in_file->IsZombie())
  {
    cout << "\n\nError reading input root filen\n";
    exit(-1);
  }
  gDirectory->cd();

  Int_t pos = TString(in_file->GetName()).First(".");

  cout << "Computing resolution with run " << TString(in_file->GetName())(pos - 8, 8) << endl;

  TTree *myTree = (TTree *)gDirectory->Get("tree");

  Double_t tmp_val[msd_stations][3], l_distances[msd_stations];

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
      myTree->SetBranchAddress("Xmeas", tmp_val);
    }
  }
  else
  {
    myTree->SetBranchAddress("MSDmeas", tmp_val);
  }

  myTree->GetEntry(0);

  for (int i = 0; i < msd_stations; i++)
  {
    l_distances[i] = tmp_val[i][2];
  }

  Resolution(myTree, correct_position, rotation_angle, l_distances, msd_stations, verbose);

  cout << "\n\nProgram completed correctly \n\n";
  return 0;
}
