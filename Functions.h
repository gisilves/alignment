#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>

using namespace std;

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TProfile.h"
#include "TMath.h"

void Alignment_1(TTree *data_tree, Double_t correct_position[][2], int msd_stations, bool verbose);
void Alignment_2(TTree *data_tree, Double_t correct_position[][2], int msd_stations, bool verbose);
void Alignment_3(TTree *data_tree, Double_t correct_position[][2], Double_t rotation_angle[][2], Double_t l_distances[], int msd_stations, int iterations, bool verbose);
void Resolution(TTree *data_tree, Double_t correct_position[][2], Double_t rotation_angle[][2], Double_t l_distances[], int msd_stations, bool verbose);
void compute_line(Double_t l_distances[], Double_t center_of_gravity[], Double_t &ang_coeff, Double_t &intercept, int msd_stations, int to_exclude);

void Alignment_1(TTree *data_tree, Double_t correct_position[][2], int msd_stations, bool verbose)
{
  Int_t n_entries = data_tree->GetEntries();
  Double_t tmp_val[msd_stations][3];
  Bool_t good_layer[msd_stations][2], proceed_track;

  TH1F *h_distribution[msd_stations][2];
  for (int idx_l = 0; idx_l < msd_stations; idx_l++)
  {
    TString h_distribution_x_name = "h_distribution_X_layer_";
    TString h_distribution_y_name = "h_distribution_Y_layer_";
    h_distribution_x_name += idx_l;
    h_distribution_y_name += idx_l;
    h_distribution[idx_l][0] = new TH1F(h_distribution_x_name, h_distribution_x_name, 250, 0, -1);
    h_distribution[idx_l][1] = new TH1F(h_distribution_y_name, h_distribution_y_name, 250, 0, -1);
  }

  if (data_tree->SetBranchAddress("MSDmeas", tmp_val) < 0)
  {
    cout << "ERROR: MSDmeas branch is not present, checking for Xmeas branch" << endl;
    if (data_tree->SetBranchAddress("Xmeas", tmp_val) < 0)
    {
      cout << "ERROR: neither Xmeas nor MSDmeas branch is present" << endl;
      return;
    }
    else
    {
      cout << "Found Xmeas branch" << endl;
      data_tree->SetBranchAddress("Xmeas", tmp_val);
    }
  }
  else
  {
    cout << "Found MSDmeas branch" << endl;
    data_tree->SetBranchAddress("MSDmeas", tmp_val);
  }

  TFile *results = new TFile("out_files/Alignment_1_results.root", "RECREATE");

  for (int idx = 0; idx < n_entries; idx++)
  {
    data_tree->GetEntry(idx);

    for (int idx_l = 0; idx_l < msd_stations; idx_l++)
      for (int idx_p = 0; idx_p < 2; idx_p++)
        good_layer[idx_l][idx_p] = false;

    proceed_track = true;

    for (int idx_l = 0; idx_l < msd_stations; idx_l++)
      for (int idx_p = 0; idx_p < 2; idx_p++)
        if (tmp_val[idx_l][idx_p] != (-999))
          good_layer[idx_l][idx_p] = true;

    for (int idx_l = 0; idx_l < msd_stations; idx_l++)
      for (int idx_p = 0; idx_p < 2; idx_p++)
        if (good_layer[idx_l][idx_p] == false)
        {
          proceed_track = false;
          break;
        }

    if (proceed_track == true)
    {

      for (int idx_l = 0; idx_l < msd_stations; idx_l++)
        for (int idx_p = 0; idx_p < 2; idx_p++)
          h_distribution[idx_l][idx_p]->Fill(tmp_val[idx_l][idx_p]);

    } //end if "proceed_track==true"
  }   //end for on entries

  TF1 *gaus = new TF1("gaus", "gaus", -30, 30);
  for (int idx_l = 0; idx_l < msd_stations; idx_l++)
    for (int idx_p = 0; idx_p < 2; idx_p++)
    {
      h_distribution[idx_l][idx_p]->Fit(gaus, "QR", "", h_distribution[idx_l][idx_p]->GetMean() - h_distribution[idx_l][idx_p]->GetRMS(), h_distribution[idx_l][idx_p]->GetMean() - h_distribution[idx_l][idx_p]->GetRMS());
      h_distribution[idx_l][idx_p]->Write();
    }

  for (int idx_l = 1; idx_l < msd_stations; idx_l++)
    for (int idx_p = 0; idx_p < 2; idx_p++)
      if (h_distribution[idx_l][idx_p]->GetMean() > h_distribution[0][idx_p]->GetMean())
        correct_position[idx_l - 1][idx_p] = -TMath::Abs(h_distribution[idx_l][idx_p]->GetMean() - h_distribution[0][idx_p]->GetMean());
      else
        correct_position[idx_l - 1][idx_p] = TMath::Abs(h_distribution[0][idx_p]->GetMean() - h_distribution[idx_l][idx_p]->GetMean());

  for (int idx_l = 0; idx_l < msd_stations; idx_l++)
    for (int idx_p = 0; idx_p < 2; idx_p++)
      h_distribution[idx_l][idx_p]->Delete();

  gaus->Delete();

  results->Close();

} //end Alignment_1 function

void Alignment_2(TTree *data_tree, Double_t correct_position[][2], int msd_stations, bool verbose)
{

  Int_t n_entries = data_tree->GetEntries();
  Double_t tmp_val[msd_stations][3];
  Bool_t good_layer[msd_stations][2], proceed_track;

  TH1F *h_difference[msd_stations - 1][2];

  for (int idx_l = 0; idx_l < msd_stations; idx_l++)
  {
    TString h_difference_x_name = "h_difference_X_layer_";
    TString h_difference_y_name = "h_difference_Y_layer_";
    h_difference_x_name += idx_l;
    h_difference_y_name += idx_l;
    h_difference[idx_l][0] = new TH1F(h_difference_x_name, h_difference_x_name, 250, 0, -1);
    h_difference[idx_l][1] = new TH1F(h_difference_y_name, h_difference_y_name, 250, 0, -1);
  }

  TFile *results = new TFile("out_files/Alignment_2_results.root", "RECREATE");

  if (data_tree->SetBranchAddress("MSDmeas", tmp_val) < 0)
  {
    cout << "ERROR: MSDmeas branch is not present, checking for Xmeas branch" << endl;
    if (data_tree->SetBranchAddress("Xmeas", tmp_val) < 0)
    {
      cout << "ERROR: neither Xmeas nor MSDmeas branch is present" << endl;
      return;
    }
    else
    {
      cout << "Found Xmeas branch" << endl;
      data_tree->SetBranchAddress("Xmeas", tmp_val);
    }
  }
  else
  {
    cout << "Found MSDmeas branch" << endl;
    data_tree->SetBranchAddress("MSDmeas", tmp_val);
  }
  for (int idx = 0; idx < n_entries; idx++)
  {
    data_tree->GetEntry(idx);

    for (int idx_l = 0; idx_l < msd_stations; idx_l++)
      for (int idx_p = 0; idx_p < 2; idx_p++)
        good_layer[idx_l][idx_p] = false;

    proceed_track = true;

    for (int idx_l = 0; idx_l < msd_stations; idx_l++)
      for (int idx_p = 0; idx_p < 2; idx_p++)
        if (tmp_val[idx_l][idx_p] != (-999))
          good_layer[idx_l][idx_p] = true;

    for (int idx_l = 0; idx_l < msd_stations; idx_l++)
      for (int idx_p = 0; idx_p < 2; idx_p++)
        if (good_layer[idx_l][idx_p] == false)
        {
          proceed_track = false;
          break;
        }

    if (proceed_track == true)
    {

      //Correct the position using the informations obtained at Step-1

      for (int idx_l = 1; idx_l < msd_stations; idx_l++)
        for (int idx_p = 0; idx_p < 2; idx_p++)
          tmp_val[idx_l][idx_p] += correct_position[idx_l - 1][idx_p];

      for (int idx_l = 1; idx_l < msd_stations; idx_l++)
        for (int idx_p = 0; idx_p < 2; idx_p++)
          h_difference[idx_l - 1][idx_p]->Fill(tmp_val[idx_l][idx_p] - tmp_val[idx_l - 1][idx_p]);

    } //end if "proceed_track==true"
  }   //end for on entries

  TF1 *gaus = new TF1("gaus", "gaus", -30, 30);
  for (int idx_l = 1; idx_l < msd_stations; idx_l++)
    for (int idx_p = 0; idx_p < 2; idx_p++)
    {
      h_difference[idx_l - 1][idx_p]->Fit(gaus, "QR", "", h_difference[idx_l - 1][idx_p]->GetMean() - h_difference[idx_l - 1][idx_p]->GetRMS(), h_difference[idx_l - 1][idx_p]->GetMean() - h_difference[idx_l - 1][idx_p]->GetRMS());
      h_difference[idx_l - 1][idx_p]->Write();
    }

  if (verbose)
    cout << "\n\n";

  for (int idx_l = 1; idx_l < msd_stations; idx_l++)
    for (int idx_p = 0; idx_p < 2; idx_p++)
      if (h_difference[idx_l - 1][idx_p]->GetMean() > 0)
      {
        correct_position[idx_l - 1][idx_p] += -h_difference[idx_l - 1][idx_p]->GetMean();
        if (verbose)
          cout << "layer " << idx_l << " shifted in direction " << idx_p << " of " << h_difference[idx_l - 1][idx_p]->GetMean() << " mm" << endl;
      }
      else
      {
        correct_position[idx_l - 1][idx_p] += -h_difference[idx_l - 1][idx_p]->GetMean();
        if (verbose)
          cout << "layer " << idx_l << " shifted in direction " << idx_p << " of " << h_difference[idx_l - 1][idx_p]->GetMean() << " mm" << endl;
      }

  for (int idx_l = 1; idx_l < msd_stations; idx_l++)
    for (int idx_p = 0; idx_p < 2; idx_p++)
      h_difference[idx_l - 1][idx_p]->Delete();

  gaus->Delete();

  results->Close();

} //end Alignment_2 function

void Alignment_3(TTree *data_tree, Double_t correct_position[][2], Double_t rotation_angle[][2], Double_t l_distances[], int msd_stations, int iterations, bool verbose)
{
  if (verbose)
    cout << "Trying to align " << msd_stations << " msd_stations" << endl;

  Int_t n_entries = data_tree->GetEntries(), iteration_number = iterations;
  Double_t tmp_val[msd_stations][3], tmp_rotation[msd_stations][2], center_of_gravity[msd_stations], intercept, ang_coeff, residual, delta_angle[msd_stations - 1][2];
  Bool_t good_layer[msd_stations][2], proceed_track;

  if (data_tree->SetBranchAddress("MSDmeas", tmp_val) < 0)
  {
    cout << "ERROR: MSDmeas branch is not present, checking for Xmeas branch" << endl;
    if (data_tree->SetBranchAddress("Xmeas", tmp_val) < 0)
    {
      cout << "ERROR: neither Xmeas nor MSDmeas branch is present" << endl;
      return;
    }
    else
    {
      cout << "Found Xmeas branch" << endl;
      data_tree->SetBranchAddress("Xmeas", tmp_val);
    }
  }
  else
  {
    cout << "Found MSDmeas branch" << endl;
    data_tree->SetBranchAddress("MSDmeas", tmp_val);
  }
  TFile *results = new TFile("out_files/Alignment_3_results.root", "RECREATE");
  TF1 *line = new TF1("line", "pol1", -30, 30);

  for (int idx_l = 0; idx_l < msd_stations - 1; idx_l++)
    for (int idx_p = 0; idx_p < 2; idx_p++)
      rotation_angle[idx_l][idx_p] = 0;

  for (int idx_it = 0; idx_it < iteration_number; idx_it++)
  {

    TProfile h_rotation[msd_stations][2];
    TH1F h_residual[msd_stations][2];
    for (int idx_l = 0; idx_l < msd_stations; idx_l++)
    {
      TString h_residual_x_name = "h_residual_X_layer_";
      TString h_residual_y_name = "h_residual_Y_layer_";
      TString h_rotation_x_name = "h_rotation_X_layer_";
      TString h_rotation_y_name = "h_rotation_Y_layer_";
      h_residual_x_name += idx_l;
      h_residual_y_name += idx_l;
      h_rotation_x_name += idx_l;
      h_rotation_y_name += idx_l;

      h_residual_x_name += "_iteration_";
      h_residual_y_name += "_iteration_";
      h_rotation_x_name += "_iteration_";
      h_rotation_y_name += "_iteration_";

      h_residual_x_name += idx_it;
      h_residual_y_name += idx_it;
      h_rotation_x_name += idx_it;
      h_rotation_y_name += idx_it;

      h_residual[idx_l][0] = TH1F(h_residual_x_name, h_residual_x_name, 500, 0, -1);
      h_residual[idx_l][1] = TH1F(h_residual_y_name, h_residual_y_name, 500, 0, -1);
      h_rotation[idx_l][0] = TProfile(h_rotation_x_name, h_rotation_x_name, 1000, -2.5, 2.5, -5, 5);
      h_rotation[idx_l][1] = TProfile(h_rotation_y_name, h_rotation_y_name, 1000, -2.5, 2.5, -5, 5);
    }

    for (int idx = 0; idx < n_entries; idx++)
    {
      data_tree->GetEntry(idx);

      for (int idx_l = 0; idx_l < msd_stations; idx_l++)
        for (int idx_p = 0; idx_p < 2; idx_p++)
          good_layer[idx_l][idx_p] = false;

      proceed_track = true;

      for (int idx_l = 0; idx_l < msd_stations; idx_l++)
        for (int idx_p = 0; idx_p < 2; idx_p++)
          if (tmp_val[idx_l][idx_p] != (-999))
            good_layer[idx_l][idx_p] = true;

      for (int idx_l = 0; idx_l < msd_stations; idx_l++)
        for (int idx_p = 0; idx_p < 2; idx_p++)
          if (good_layer[idx_l][idx_p] == false)
          {
            proceed_track = false;
            break;
          }

      if (proceed_track == true)
      {
        //Correct measures
        for (int idx_l = 1; idx_l < msd_stations; idx_l++)
        {
          tmp_val[idx_l][0] += correct_position[idx_l - 1][0];
          tmp_val[idx_l][1] += correct_position[idx_l - 1][1];
          if (idx_it != 0)
          {
            tmp_val[idx_l][0] = TMath::Cos(rotation_angle[idx_l - 1][0]) * tmp_val[idx_l][0] - TMath::Sin(rotation_angle[idx_l - 1][0]) * tmp_val[idx_l][1];
            tmp_val[idx_l][1] = TMath::Sin(rotation_angle[idx_l - 1][1]) * tmp_val[idx_l][0] + TMath::Cos(rotation_angle[idx_l - 1][1]) * tmp_val[idx_l][1];
          }
        }

        for (int idx_p = 0; idx_p < 2; idx_p++)
        {
          for (int idx_l = 0; idx_l < msd_stations; idx_l++)
            center_of_gravity[idx_l] = tmp_val[idx_l][idx_p];

          compute_line(l_distances, center_of_gravity, ang_coeff, intercept, msd_stations, -1);

          for (int idx_l = 0; idx_l < msd_stations; idx_l++)
          {
            residual = center_of_gravity[idx_l] - (ang_coeff * l_distances[idx_l] + intercept);
            h_residual[idx_l][idx_p].Fill(residual);
            if (idx_p == 0)
              h_rotation[idx_l][idx_p].Fill(tmp_val[idx_l][1], residual);
            else
              h_rotation[idx_l][idx_p].Fill(tmp_val[idx_l][0], residual);
          }
        } //end for on positions
      }   //end if "proceed_track==true"
    }     //end for on entries

    if (verbose)
      cout << "\n\nHere the shifting informations: \n";
    for (int idx_l = 1; idx_l < msd_stations; idx_l++)
      for (int idx_p = 0; idx_p < 2; idx_p++)
        if (h_residual[idx_l][idx_p].GetMean() > h_residual[0][idx_p].GetMean())
        {
          correct_position[idx_l - 1][idx_p] += -(h_residual[idx_l][idx_p].GetMean() - h_residual[0][idx_p].GetMean());
          if (verbose)
            cout << "layer " << idx_l << " shifted in direction " << idx_p << " by " << -(h_residual[idx_l][idx_p].GetMean() - h_residual[0][idx_p].GetMean()) << " (mm)" << endl;
        }
        else
        {
          correct_position[idx_l - 1][idx_p] += h_residual[0][idx_p].GetMean() - h_residual[idx_l][idx_p].GetMean();
          if (verbose)
            cout << "layer " << idx_l << " shifted in direction " << idx_p << " by " << h_residual[0][idx_p].GetMean() - h_residual[idx_l][idx_p].GetMean() << " (mm)" << endl;
        }

    //Filling tmp_rotation matrix...
    for (int idx_l = 0; idx_l < msd_stations; idx_l++)
      for (int idx_p = 0; idx_p < 2; idx_p++)
      {
        h_rotation[idx_l][idx_p].Fit("line", "QR", "", h_rotation[idx_l][idx_p].GetMean() - 2 * h_rotation[idx_l][idx_p].GetRMS(), h_rotation[idx_l][idx_p].GetMean() + 2 * h_rotation[idx_l][idx_p].GetRMS());
        if (idx_p == 0)
          tmp_rotation[idx_l][idx_p] = TMath::ASin(-line->GetParameter(1));
        else
          tmp_rotation[idx_l][idx_p] = TMath::ASin(line->GetParameter(1));
      }

    //Filling rotation_angle matrix, after have computed the angle to rotate each layer respect layer "0"

    for (int idx_l = 1; idx_l < msd_stations; idx_l++)
      for (int idx_p = 0; idx_p < 2; idx_p++)
      {
        delta_angle[idx_l - 1][idx_p] = rotation_angle[idx_l - 1][idx_p];
        rotation_angle[idx_l - 1][idx_p] += (tmp_rotation[0][idx_p] - tmp_rotation[idx_l][idx_p]);
        delta_angle[idx_l - 1][idx_p] -= rotation_angle[idx_l - 1][idx_p];
      }

    if (verbose)
      cout << "\n\nHere the 'tmp_rotation' matrix\n";
    for (int idx_l = 0; idx_l < msd_stations; idx_l++)
    {
      for (int idx_p = 0; idx_p < 2; idx_p++)
        if (verbose)
          cout << setw(5) << tmp_rotation[idx_l][idx_p] << "\t";
      if (verbose)
        cout << endl;
    }

    cout << "Iteration step " << idx_it << endl;

    if (idx_it == (iteration_number - 1) || idx_it == 0)
    {
      cout << "\n\nHere the 'rotation_angle' matrix: \n";

      for (int idx_l = 0; idx_l < msd_stations; idx_l++)
      {
        cout << setw(5) << rotation_angle[idx_l][0] << "\t" << rotation_angle[idx_l][1] << endl;

        for (int idx_p = 0; idx_p < 2; idx_p++)
        {
          h_residual[idx_l][idx_p].Fit("gaus", "QR", "", h_residual[idx_l][idx_p].GetMean() - h_residual[idx_l][idx_p].GetRMS(), h_residual[idx_l][idx_p].GetMean() - h_residual[idx_l][idx_p].GetRMS());
          h_residual[idx_l][idx_p].Write();
          h_rotation[idx_l][idx_p].Write();
        }
      }
    }

  } //end for on iterations

  results->Close();
  line->Delete();

} //end Alignment_3 function

void Resolution(TTree *data_tree, Double_t correct_position[][2], Double_t rotation_angle[][2], Double_t l_distances[], int msd_stations, bool verbose)
{

  Int_t n_entries = data_tree->GetEntries();
  Double_t tmp_val[msd_stations][3], tmp_rotation[msd_stations][2], center_of_gravity[msd_stations], intercept, ang_coeff, residual;
  Bool_t good_layer[msd_stations][2], proceed_track;

  if (data_tree->SetBranchAddress("MSDmeas", tmp_val) < 0)
  {
    cout << "ERROR: MSDmeas branch is not present, checking for Xmeas branch" << endl;
    if (data_tree->SetBranchAddress("Xmeas", tmp_val) < 0)
    {
      cout << "ERROR: neither Xmeas nor MSDmeas branch is present" << endl;
      return;
    }
    else
    {
      cout << "Found Xmeas branch" << endl;
      data_tree->SetBranchAddress("Xmeas", tmp_val);
    }
  }
  else
  {
    cout << "Found MSDmeas branch" << endl;
    data_tree->SetBranchAddress("MSDmeas", tmp_val);
  }
  TFile *results = new TFile("out_files/Resolution_biased_results.root", "RECREATE");
  TF1 *line = new TF1("line", "pol1", -30, 30);
  TF1 *fittedgaus;

  TH1F h_residual[msd_stations][2];
  for (int idx_l = 0; idx_l < msd_stations; idx_l++)
  {
    TString h_residual_x_name = "h_residual_X_layer_";
    TString h_residual_y_name = "h_residual_Y_layer_";
    h_residual_x_name += idx_l;
    h_residual_y_name += idx_l;
    h_residual[idx_l][0] = TH1F(h_residual_x_name, h_residual_x_name, 500, 0, -1);
    h_residual[idx_l][1] = TH1F(h_residual_y_name, h_residual_y_name, 500, 0, -1);
  }

  for (int idx = 0; idx < n_entries; idx++)
  {
    data_tree->GetEntry(idx);

    for (int idx_l = 0; idx_l < msd_stations; idx_l++)
      for (int idx_p = 0; idx_p < 2; idx_p++)
        good_layer[idx_l][idx_p] = false;

    proceed_track = true;

    for (int idx_l = 0; idx_l < msd_stations; idx_l++)
      for (int idx_p = 0; idx_p < 2; idx_p++)
        if (tmp_val[idx_l][idx_p] != (-999))
          good_layer[idx_l][idx_p] = true;

    for (int idx_l = 0; idx_l < msd_stations; idx_l++)
      for (int idx_p = 0; idx_p < 2; idx_p++)
        if (good_layer[idx_l][idx_p] == false)
        {
          proceed_track = false;
          break;
        }

    if (proceed_track == true)
    {

      //Correct measures
      for (int idx_l = 1; idx_l < msd_stations; idx_l++)
      {
        tmp_val[idx_l][0] += correct_position[idx_l - 1][0];
        tmp_val[idx_l][1] += correct_position[idx_l - 1][1];

        tmp_val[idx_l][0] = TMath::Cos(rotation_angle[idx_l - 1][0]) * tmp_val[idx_l][0] - TMath::Sin(rotation_angle[idx_l - 1][0]) * tmp_val[idx_l][1];
        tmp_val[idx_l][1] = TMath::Sin(rotation_angle[idx_l - 1][1]) * tmp_val[idx_l][0] + TMath::Cos(rotation_angle[idx_l - 1][1]) * tmp_val[idx_l][1];
      }

      for (int idx_p = 0; idx_p < 2; idx_p++)
      {
        for (int idx_l = 0; idx_l < msd_stations; idx_l++)
          center_of_gravity[idx_l] = tmp_val[idx_l][idx_p];

        for (int idx_l = 0; idx_l < msd_stations; idx_l++)
        {
          compute_line(l_distances, center_of_gravity, ang_coeff, intercept, msd_stations, idx_l);

          residual = center_of_gravity[idx_l] - (ang_coeff * l_distances[idx_l] + intercept);
          h_residual[idx_l][idx_p].Fill(residual);
        }
      } //end for on positions
    }   //end if "proceed_track==true"
  }     //end for on entries

  Double_t alpha = 1;
  Double_t mean_z = 0.;
  Double_t mean_z2 = 0.;

  for (int idx_p = 0; idx_p < 2; idx_p++)
  {
    mean_z = 0;
    mean_z2 = 0;

    for (int idx_l = 0; idx_l < msd_stations; idx_l++)
    {
      mean_z += l_distances[idx_l];
      mean_z2 += l_distances[idx_l] * l_distances[idx_l];
    }

    mean_z = mean_z / msd_stations;
    mean_z2 = mean_z2 / msd_stations;

    for (int idx_l = 0; idx_l < msd_stations; idx_l++)
    {
      h_residual[idx_l][idx_p].Fit("gaus", "QR", "", h_residual[idx_l][idx_p].GetMean() - h_residual[idx_l][idx_p].GetRMS(), h_residual[idx_l][idx_p].GetMean() - h_residual[idx_l][idx_p].GetRMS());
      fittedgaus = (TF1 *)h_residual[idx_l][idx_p].GetListOfFunctions()->FindObject("gaus");
      h_residual[idx_l][idx_p].Write();

      alpha = 1 + (mean_z * mean_z + l_distances[idx_l] - 2 * l_distances[idx_l] * mean_z) / (msd_stations * (mean_z2 - mean_z * mean_z));

      cout << "Resolution for detector " << 2 * idx_l + idx_p + 1 << " is " << fittedgaus->GetParameter(2) << endl;
      cout << "Resolution for detector (after correction) " << 2 * idx_l + idx_p + 1 << " is " << fittedgaus->GetParameter(2) / sqrt(alpha) << endl;

    }

    cout << endl;
  }

  results->Close();
  line->Delete();

} //end Resolution_biased function

void compute_line(Double_t l_distances[], Double_t center_of_gravity[], Double_t &ang_coeff, Double_t &intercept, int msd_stations, int to_exclude)
{

  Double_t A = 0, B = 0, C = 0, D = 0, E = 0, F = 0;

  for (int idx_l = 0; idx_l < msd_stations; idx_l++)
  {
    if (idx_l == to_exclude)
      continue;

    A += l_distances[idx_l];
    B += 1;
    C += center_of_gravity[idx_l];
    D += pow(l_distances[idx_l], 2);
    E += l_distances[idx_l] * center_of_gravity[idx_l];
    F += pow(center_of_gravity[idx_l], 2);
  }

  intercept = (D * C - E * A) / (D * B - pow(A, 2));
  ang_coeff = (E * B - C * A) / (D * B - pow(A, 2));
}
