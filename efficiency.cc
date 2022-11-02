#include <iostream>
#include <fstream>
#include <vector>
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TEfficiency.h"

void efficiency(){
    TFile *file = new TFile("img1006/hist1006_2.root");
    
    TH1D *probe_muon_endcap_momentum_thr1_cut_hist = (TH1D*)file->Get("probe_muon_endcap_momentum_thr1_cut_hist");
    TH1D *probe_muon_endcap_momentum_thr2_cut_hist = (TH1D*)file->Get("probe_muon_endcap_momentum_thr2_cut_hist");
    TH1D *probe_muon_endcap_momentum_thr3_cut_hist = (TH1D*)file->Get("probe_muon_endcap_momentum_thr3_cut_hist");
    TH1D *probe_muon_endcap_momentum_thr4_cut_hist = (TH1D*)file->Get("probe_muon_endcap_momentum_thr4_cut_hist");
    TH1D *probe_muon_endcap_momentum_thr5_cut_hist = (TH1D*)file->Get("probe_muon_endcap_momentum_thr5_cut_hist");
    TH1D *probe_muon_endcap_momentum_thr6_cut_hist = (TH1D*)file->Get("probe_muon_endcap_momentum_thr6_cut_hist");
    TH1D *probe_muon_barrel_momentum_thr1_cut_hist = (TH1D*)file->Get("probe_muon_barrel_momentum_thr1_cut_hist");
    TH1D *probe_muon_barrel_momentum_thr2_cut_hist = (TH1D*)file->Get("probe_muon_barrel_momentum_thr2_cut_hist");
    TH1D *probe_muon_barrel_momentum_thr3_cut_hist = (TH1D*)file->Get("probe_muon_barrel_momentum_thr3_cut_hist");
    TH1D *probe_muon_barrel_momentum_thr4_cut_hist = (TH1D*)file->Get("probe_muon_barrel_momentum_thr4_cut_hist");
    TH1D *probe_muon_barrel_momentum_thr5_cut_hist = (TH1D*)file->Get("probe_muon_barrel_momentum_thr5_cut_hist");
    TH1D *probe_muon_barrel_momentum_thr6_cut_hist = (TH1D*)file->Get("probe_muon_barrel_momentum_thr6_cut_hist");
    TH1D *probe_muon_endcap_momentum_cut_hist = (TH1D*)file->Get("probe_muon_endcap_momentum_cut_hist");
    TH1D *probe_muon_endcap_momentum_hist = (TH1D*)file->Get("probe_muon_endcap_momentum_hist");
    TH1D *probe_muon_endcap_deltaEta_cut_hist = (TH1D*)file->Get("probe_muon_endcap_deltaEta_cut_hist");
    TH1D *probe_muon_endcap_deltaEta_hist = (TH1D*)file->Get("probe_muon_endcap_deltaEta_hist");
    TH1D *probe_muon_endcap_deltaPhi_cut_hist = (TH1D*)file->Get("probe_muon_endcap_deltaPhi_cut_hist");
    TH1D *probe_muon_endcap_deltaPhi_hist = (TH1D*)file->Get("probe_muon_endcap_deltaPhi_hist");
    TH1D *probe_muon_endcap_deltaR_cut_hist = (TH1D*)file->Get("probe_muon_endcap_deltaR_cut_hist");
    TH1D *probe_muon_endcap_deltaR_hist = (TH1D*)file->Get("probe_muon_endcap_deltaR_hist");
    TH1D *probe_muon_barrel_momentum_cut_hist = (TH1D*)file->Get("probe_muon_barrel_momentum_cut_hist");
    TH1D *probe_muon_barrel_momentum_hist = (TH1D*)file->Get("probe_muon_barrel_momentum_hist");
    TH1D *probe_muon_barrel_deltaEta_cut_hist = (TH1D*)file->Get("probe_muon_barrel_deltaEta_cut_hist");
    TH1D *probe_muon_barrel_deltaEta_hist = (TH1D*)file->Get("probe_muon_barrel_deltaEta_hist");
    TH1D *probe_muon_barrel_deltaPhi_cut_hist = (TH1D*)file->Get("probe_muon_barrel_deltaPhi_cut_hist");
    TH1D *probe_muon_barrel_deltaPhi_hist = (TH1D*)file->Get("probe_muon_barrel_deltaPhi_hist");
    TH1D *probe_muon_barrel_deltaR_cut_hist = (TH1D*)file->Get("probe_muon_barrel_deltaR_cut_hist");
    TH1D *probe_muon_barrel_deltaR_hist = (TH1D*)file->Get("probe_muon_barrel_deltaR_hist");
    TH1D *L1_probe_muon_momentum_cut_hist = (TH1D*)file->Get("L1_probe_muon_momentum_cut_hist");
    TH1D *probe_muon_momentum_hist = (TH1D*)file->Get("probe_muon_momentum_hist");
    TH1D *probe_thrNum1_pt = (TH1D*)file->Get("probe_thrNum1_pt");
    TH1D *probe_thrNum2_pt = (TH1D*)file->Get("probe_thrNum2_pt");
    TH1D *probe_thrNum3_pt = (TH1D*)file->Get("probe_thrNum3_pt");
    TH1D *probe_thrNum4_pt = (TH1D*)file->Get("probe_thrNum4_pt");
    TH1D *probe_thrNum5_pt = (TH1D*)file->Get("probe_thrNum5_pt");
    TH1D *probe_thrNum6_pt = (TH1D*)file->Get("probe_thrNum6_pt");
    TH1D *probe_muon_deltaR_thr1_hist = (TH1D*)file->Get("probe_muon_deltaR_thr1_hist");
    TH1D *probe_muon_deltaR_thr2_hist = (TH1D*)file->Get("probe_muon_deltaR_thr2_hist");
    TH1D *probe_muon_deltaR_thr3_hist = (TH1D*)file->Get("probe_muon_deltaR_thr3_hist");
    TH1D *probe_muon_deltaR_thr4_hist = (TH1D*)file->Get("probe_muon_deltaR_thr4_hist");
    TH1D *probe_muon_deltaR_thr5_hist = (TH1D*)file->Get("probe_muon_deltaR_thr5_hist");
    TH1D *probe_muon_deltaR_thr6_hist = (TH1D*)file->Get("probe_muon_deltaR_thr6_hist");

    TCanvas *canvas1 = new TCanvas();
    TCanvas *canvas2 = new TCanvas();
    TCanvas *canvas3 = new TCanvas();
    TCanvas *canvas4 = new TCanvas();
    TCanvas *canvas5 = new TCanvas();
    TCanvas *canvas6 = new TCanvas();
    TCanvas *canvas7 = new TCanvas();
    TCanvas *canvas8 = new TCanvas();
    TCanvas *canvas9 = new TCanvas();
    TCanvas *canvas10 = new TCanvas();
    TCanvas *canvas11 = new TCanvas();
    TCanvas *canvas12 = new TCanvas();
    TCanvas *canvas13 = new TCanvas();
    
    // efficiency

    TEfficiency *pEff_endcap = new TEfficiency(*probe_muon_endcap_momentum_cut_hist, *probe_muon_endcap_momentum_hist);
    pEff_endcap->SetTitle("efficiency in endcap;pt[GeV];efficiency");

    //separate by thr Number
    TEfficiency *pEff_endcap_thr1 = new TEfficiency(*probe_muon_endcap_momentum_thr1_cut_hist, *probe_muon_endcap_momentum_hist);
    pEff_endcap->SetTitle("efficiency in endcap;pt[GeV];efficiency");

    TEfficiency *pEff_endcap_thr2 = new TEfficiency(*probe_muon_endcap_momentum_thr2_cut_hist, *probe_muon_endcap_momentum_hist);
    pEff_endcap->SetTitle("efficiency in endcap;pt[GeV];efficiency");

    TEfficiency *pEff_endcap_thr3 = new TEfficiency(*probe_muon_endcap_momentum_thr3_cut_hist, *probe_muon_endcap_momentum_hist);
    pEff_endcap->SetTitle("efficiency in endcap;pt[GeV];efficiency");

    TEfficiency *pEff_endcap_thr4 = new TEfficiency(*probe_muon_endcap_momentum_thr4_cut_hist, *probe_muon_endcap_momentum_hist);
    pEff_endcap->SetTitle("efficiency in endcap;pt[GeV];efficiency");

    TEfficiency *pEff_endcap_thr5 = new TEfficiency(*probe_muon_endcap_momentum_thr5_cut_hist, *probe_muon_endcap_momentum_hist);
    pEff_endcap->SetTitle("efficiency in endcap;pt[GeV];efficiency");

    TEfficiency *pEff_endcap_thr6 = new TEfficiency(*probe_muon_endcap_momentum_thr6_cut_hist, *probe_muon_endcap_momentum_hist);
    pEff_endcap->SetTitle("efficiency in endcap;pt[GeV];efficiency");

    TEfficiency *eta_Eff_endcap = new TEfficiency(*probe_muon_endcap_deltaEta_cut_hist, *probe_muon_endcap_deltaEta_hist);
    eta_Eff_endcap->SetTitle("efficiency in endcap;eta;efficiency");

    TEfficiency *phi_Eff_endcap = new TEfficiency(*probe_muon_endcap_deltaPhi_cut_hist, *probe_muon_endcap_deltaPhi_hist);
    phi_Eff_endcap->SetTitle("efficiency in endcap;phi;efficiency");

    TEfficiency *R_Eff_endcap = new TEfficiency(*probe_muon_endcap_deltaR_cut_hist, *probe_muon_endcap_deltaR_hist);
    R_Eff_endcap->SetTitle("efficiency in endcap;pt[GeV];efficiency");

    TEfficiency *pEff_barrel = new TEfficiency(*probe_muon_barrel_momentum_cut_hist, *probe_muon_barrel_momentum_hist);
    pEff_barrel->SetTitle("efficiency in barrel;pt[GeV];efficiency");

    TEfficiency *eta_Eff_barrel = new TEfficiency(*probe_muon_barrel_deltaEta_cut_hist, *probe_muon_barrel_deltaEta_hist);
    eta_Eff_barrel->SetTitle("efficiency in barrel;eta;efficiency");

    TEfficiency *phi_Eff_barrel = new TEfficiency(*probe_muon_barrel_deltaPhi_cut_hist, *probe_muon_barrel_deltaPhi_hist);
    phi_Eff_barrel->SetTitle("efficiency in barrel;phi;efficiency");

    TEfficiency *R_Eff_barrel = new TEfficiency(*probe_muon_barrel_deltaR_cut_hist, *probe_muon_barrel_deltaR_hist);
    R_Eff_barrel->SetTitle("efficiency in barrel;pt[GeV];efficiency");

    TEfficiency *pEff_L1 = new TEfficiency(*L1_probe_muon_momentum_cut_hist, *probe_muon_momentum_hist);
    pEff_L1->SetTitle("L1 efficiency;pt[GeV];efficiency");

    //Draw
    canvas1->cd();
    probe_muon_endcap_momentum_thr1_cut_hist->Draw();
    probe_muon_endcap_momentum_thr2_cut_hist->SetLineColor(3);
    probe_muon_endcap_momentum_thr2_cut_hist->Draw("same");
    probe_muon_endcap_momentum_thr3_cut_hist->SetLineColor(4);
    probe_muon_endcap_momentum_thr3_cut_hist->Draw("same");
    probe_muon_endcap_momentum_thr4_cut_hist->SetLineColor(5);
    probe_muon_endcap_momentum_thr4_cut_hist->Draw("same");
    probe_muon_endcap_momentum_thr5_cut_hist->SetLineColor(6);
    probe_muon_endcap_momentum_thr5_cut_hist->Draw("same");
    probe_muon_endcap_momentum_thr6_cut_hist->SetLineColor(7);
    probe_muon_endcap_momentum_thr6_cut_hist->Draw("same");
    //canvas1->SaveAs("img1004/probe_endcap_momentum.png");

    canvas2->cd();
    probe_thrNum1_pt->Draw();
    probe_thrNum2_pt->SetLineColor(3);
    probe_thrNum2_pt->Draw("same");
    probe_thrNum3_pt->SetLineColor(4);
    probe_thrNum3_pt->Draw("same");
    probe_thrNum4_pt->SetLineColor(5);
    probe_thrNum4_pt->Draw("same");
    probe_thrNum5_pt->SetLineColor(6);
    probe_thrNum5_pt->Draw("same");
    probe_thrNum6_pt->SetLineColor(7);
    probe_thrNum6_pt->Draw("same");
    //canvas2->SaveAs("img1004/probe_momentum.png");

    canvas3->cd();
    pEff_L1->Draw();
    //canvas3->SaveAs("img1004/efficiency_L1.png");
    canvas4->cd();
    pEff_barrel->Draw();
    //canvas4->SaveAs("img1004/efficiency_barrel.png");
    canvas5->cd();
    pEff_endcap->Draw("AP");
    //canvas5->SaveAs("img1004/efficiency_endcap.png");
    canvas6->cd();
    pEff_endcap_thr1->Draw("AP");
    canvas7->cd();
    pEff_endcap_thr2->Draw("AP");
    canvas8->cd();
    pEff_endcap_thr3->Draw("AP");
    canvas9->cd();
    pEff_endcap_thr4->Draw("AP");
    canvas10->cd();
    pEff_endcap_thr5->Draw("AP");
    canvas11->cd();
    pEff_endcap_thr6->Draw("AP");
    canvas12->cd();
    probe_muon_deltaR_thr1_hist->Draw();
    probe_muon_deltaR_thr2_hist->SetLineColor(3);
    probe_muon_deltaR_thr2_hist->Draw("same");
    probe_muon_deltaR_thr3_hist->SetLineColor(4);
    probe_muon_deltaR_thr3_hist->Draw("same");
    probe_muon_deltaR_thr4_hist->SetLineColor(5);
    probe_muon_deltaR_thr4_hist->Draw("same");
    probe_muon_deltaR_thr5_hist->SetLineColor(6);
    probe_muon_deltaR_thr5_hist->Draw("same");
    probe_muon_deltaR_thr6_hist->SetLineColor(7);
    probe_muon_deltaR_thr6_hist->Draw("same"); 
    canvas13->cd();
    probe_muon_barrel_momentum_thr1_cut_hist->Draw();
    probe_muon_barrel_momentum_thr2_cut_hist->SetLineColor(3);
    probe_muon_barrel_momentum_thr2_cut_hist->Draw("same");
    probe_muon_barrel_momentum_thr3_cut_hist->SetLineColor(4);
    probe_muon_barrel_momentum_thr3_cut_hist->Draw("same");
    probe_muon_barrel_momentum_thr4_cut_hist->SetLineColor(5);
    probe_muon_barrel_momentum_thr4_cut_hist->Draw("same");
    probe_muon_barrel_momentum_thr5_cut_hist->SetLineColor(6);
    probe_muon_barrel_momentum_thr5_cut_hist->Draw("same");
    probe_muon_barrel_momentum_thr6_cut_hist->SetLineColor(7);
    probe_muon_barrel_momentum_thr6_cut_hist->Draw("same");
    
}