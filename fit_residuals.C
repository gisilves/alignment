#include <Riostream.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include "TH2.h"
#include "TGraph.h"

double gaussian(double x, double mean, double sigma)
{
    return exp(-0.5 * pow((x - mean) / sigma, 2));
}

// fitfunction is a sum of two gaussians
double fitfunction(double *x, double *par)
{
    return par[0] * gaussian(x[0], par[1], par[2]) + par[3] * gaussian(x[0], par[4], par[5]);
}

double fitgauss(double *x, double *par)
{
    return par[0] * gaussian(x[0], par[1], par[2]);
}

void fit_residuals()
{
    TFile *f = new TFile("out_files/Resolution_unbiased_results.root");
    // get histo
    TH1F *h = (TH1F *)f->Get("h_residual_X_layer_1");
    TF1 *f1 = new TF1("f1", fitfunction, -0.1, 0.1, 6);
    f1->SetParameters(1000, 0, 0.001, 100, 0, 0.01);

    h->Fit("f1", "R");
    // retrieve fit function
    TF1 *fit = h->GetFunction("f1");
    // get fit parameters
    double par[6];
    fit->GetParameters(par);

    // compute integral of fit function
    double integral = fit->Integral(-0.1, 0.1);

    //set log Y scale
    gPad->SetLogy();
    h->Draw();

    // Draw gaussian 1
    TF1 *gaus1 = new TF1("gaus1", fitgauss, -0.1, 0.1, 3);
    gaus1->SetParameters(par[0], par[1], par[2]);
    gaus1->SetLineColor(kBlue);
    gaus1->Draw("same");

    // Draw gaussian 2
    TF1 *gaus2 = new TF1("gaus2", fitgauss, -0.1, 0.1, 3);
    gaus2->SetParameters(par[3], par[4], par[5]);
    gaus2->SetLineColor(kGreen);
    gaus2->Draw("same");

    cout << "integral = " << integral << endl;
    cout << "integral_g1 = " << gaus1->Integral(-0.1, 0.1) << " corresponding to " << gaus1->Integral(-0.1, 0.1) / integral * 100 << " %" << endl;
    cout << "integral_g2 = " << gaus2->Integral(-0.1, 0.1) << " corresponding to " << gaus2->Integral(-0.1, 0.1) / integral * 100 << " %" << endl;

    // Resolution is given by weighed average of the gaussian sigmas
    cout << "Resolution from first gaussian fit = " << par[2] << endl;
    cout << "Resolution from second gaussian fit = " << par[5] << endl;
    double res = (par[2] * gaus1->Integral(-0.1, 0.1) + par[5] * gaus2->Integral(-0.1, 0.1)) / integral;
    cout << "Final resolution = " << res << endl;

}