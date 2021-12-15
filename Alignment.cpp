#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include "TFile.h"
#include "TTree.h"

#include "Functions.h"

using namespace std;

int main(int argc, char **argv)
{
  if (argc < 5)
  {
    cout << "Usage: ./Alignment input_file #stations #iterations #verbose" << endl;
    exit(-1);
  }

  int msd_stations = atoi(argv[2]);
  bool verbose = atoi(argv[4]);

  cout << "Trying to aling " << msd_stations << " msd_stations" << endl;
  //variable declaration
  Double_t correct_position[msd_stations - 1][2];

  TFile *in_file = new TFile(argv[1], "R");
  if (in_file->IsZombie())
  {
    cout << "\n\nError reading input root filen\n";
    exit(-1);
  }
  gDirectory->cd();

  TTree *myTree = (TTree *)gDirectory->Get("tree");

  Double_t tmp_val[msd_stations][3], l_distances[msd_stations];
  myTree->SetBranchAddress("Xmeas", tmp_val);
  myTree->GetEntry(0);
  
  for (int i = 0; i < msd_stations; i++)
  {
    l_distances[i] = tmp_val[i][2];
  }

  /////////////////// ALIGNMENT METHOD 1 ///////////////////

  cout << "\n\n/////////////////// ALIGNMENT METHOD 1 ///////////////////\n\n";
  Alignment_1(myTree, correct_position, msd_stations, verbose);

  cout << "\n\n'Correct Position' matrix after the first alignment process: " << endl;
  for (int idx_l = 0; idx_l < msd_stations - 1; idx_l++)
  {
    for (int idx_p = 0; idx_p < 2; idx_p++)
      cout << setw(5) << correct_position[idx_l][idx_p] << "\t";
    cout << endl;
  }

  /////////////////// ALIGNMENT METHOD 2 ///////////////////

  cout << "\n\n/////////////////// ALIGNMENT METHOD 2 ///////////////////\n\n";
  Alignment_2(myTree, correct_position, msd_stations, verbose);

  cout << "\n\n'Correct Position' matrix after the second alignment process: " << endl;
  for (int idx_l = 0; idx_l < msd_stations - 1; idx_l++)
  {
    for (int idx_p = 0; idx_p < 2; idx_p++)
      cout << setw(5) << correct_position[idx_l][idx_p] << "\t";
    cout << endl;
  }

  /////////////////// ALIGNMENT METHOD 3 ///////////////////

  cout << "\n\n/////////////////// ALIGNMENT METHOD 3 ///////////////////\n\n";
  Alignment_3(myTree, correct_position, l_distances, msd_stations, atoi(argv[3]), verbose);

  cout << "\n\n'Correct Position' matrix after the third alignment process: " << endl;
  for (int idx_l = 0; idx_l < msd_stations - 1; idx_l++)
  {
    for (int idx_p = 0; idx_p < 2; idx_p++)
      cout << setw(5) << correct_position[idx_l][idx_p] << "\t";
    cout << endl;
  }

  cout << "\n\nProgram completed correctly \n\n";
  return 0;
}
