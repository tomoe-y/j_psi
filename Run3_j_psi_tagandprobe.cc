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
#include "Dict.cc"

void Run3_j_psi_tagandprobe(){
    TChain* chain = new TChain("physics");
    //chain->Add("/mnt/susy11/data04/atlas/data16_13TeV/periodL/mu_sample/user.junpei.00311481.physics_Main.merge.NTUP_MCP.1.f758_m1714.00-00-32_L1TGCNtuple_derivated.02-04-00.root");
    //chain->Add("/mnt/susy11/data04/atlas/data16_13TeV/periodL/mu_sample/*.root");
    chain->Add("/home/toyamash/create_Ntuple/run/L1TGCNtuple.root");

    chain->SetBranchStatus("*", 0);
    chain->SetBranchStatus("mu_m", 1);
    chain->SetBranchStatus("mu_pt", 1);
    chain->SetBranchStatus("mu_eta", 1);
    chain->SetBranchStatus("mu_phi", 1);
    chain->SetBranchStatus("mu_charge", 1);
    chain->SetBranchStatus("mu_author", 1);
    chain->SetBranchStatus("mu_muonType", 1);
    chain->SetBranchStatus("trigger_info_chain", 1);
    //chain->SetBranchStatus("trigger_info_isPassed", 1);
    chain->SetBranchStatus("trigger_info_etaVec", 1);
    chain->SetBranchStatus("trigger_info_phiVec", 1);
    chain->SetBranchStatus("trigger_info_ptVec", 1);
    chain->SetBranchStatus("trig_L1_mu_eta", 1);
    chain->SetBranchStatus("trig_L1_mu_phi", 1);
    chain->SetBranchStatus("trig_L1_mu_RoINumber", 1);
    chain->SetBranchStatus("trig_L1_mu_thrNumber", 1);
    chain->SetBranchStatus("trig_L1_mu_source", 1);
    chain->SetBranchStatus("mu_ext_b_targetEtaVec", 1);
    chain->SetBranchStatus("mu_ext_b_targetPhiVec", 1);
    chain->SetBranchStatus("vxp_type", 1);
    chain->SetBranchStatus("vxp_chi2", 1);

    std::vector<float> *mu_m = 0;
    std::vector<float> *mu_pt = 0;
    std::vector<float> *mu_eta = 0;
    std::vector<float> *mu_phi = 0;
    std::vector<int> *mu_charge = 0;
    std::vector<int> *mu_author = 0;
    std::vector<int> *mu_muonType = 0;
    std::vector<std::string> *trigger_info_chain = 0;
    //std::vector<int> *trigger_info_isPassed = 0;
    std::vector<std::vector<float>> *trigger_info_etaVec = 0;
    std::vector<std::vector<float>> *trigger_info_phiVec = 0;
    std::vector<std::vector<float>> *trigger_info_ptVec = 0;
    std::vector<float> *trig_L1_mu_eta = 0;
    std::vector<float> *trig_L1_mu_phi = 0;
    std::vector<short> *trig_L1_mu_RoINumber = 0;
    std::vector<short> *trig_L1_mu_thrNumber = 0;
    std::vector<short> *trig_L1_mu_source = 0;
    std::vector<std::vector<float>> *mu_ext_b_targetEtaVec = 0;
    std::vector<std::vector<float>> *mu_ext_b_targetPhiVec = 0;
    std::vector<int> *vxp_type = 0;
    std::vector<float> *vxp_chi2 = 0;


    chain->SetBranchAddress("mu_m", &mu_m);
    chain->SetBranchAddress("mu_pt", &mu_pt);
    chain->SetBranchAddress("mu_eta", &mu_eta);
    chain->SetBranchAddress("mu_phi", &mu_phi);
    chain->SetBranchAddress("mu_charge", &mu_charge);
    chain->SetBranchAddress("mu_author", &mu_author);
    chain->SetBranchAddress("mu_muonType", &mu_muonType);
    chain->SetBranchAddress("trigger_info_chain", &trigger_info_chain);
    //chain->SetBranchAddress("trigger_info_isPassed", &trigger_info_isPassed);
    chain->SetBranchAddress("trigger_info_phiVec", &trigger_info_phiVec);
    chain->SetBranchAddress("trigger_info_etaVec", &trigger_info_etaVec);
    chain->SetBranchAddress("trigger_info_ptVec", &trigger_info_ptVec);
    chain->SetBranchAddress("trig_L1_mu_eta", &trig_L1_mu_eta);
    chain->SetBranchAddress("trig_L1_mu_phi", &trig_L1_mu_phi);
    chain->SetBranchAddress("trig_L1_mu_RoINumber", &trig_L1_mu_RoINumber);
    chain->SetBranchAddress("trig_L1_mu_thrNumber", &trig_L1_mu_thrNumber);
    chain->SetBranchAddress("trig_L1_mu_source", &trig_L1_mu_source);
    chain->SetBranchAddress("mu_ext_b_targetEtaVec", &mu_ext_b_targetEtaVec);
    chain->SetBranchAddress("mu_ext_b_targetPhiVec", &mu_ext_b_targetPhiVec);
    chain->SetBranchAddress("vxp_type", vxp_type);
    chain->SetBranchAddress("vxp_chi2", vxp_chi2);

    TH1D *mass_hist = new TH1D("mass_hist", "mass_hist", 1000, 0, 100000);
    TH1D *cut_pair_mass_hist = new TH1D("cut_pair_mass_hist", "cut_pair_mass_hist", 1000, 0, 10000);
    TH1D *pt_hist = new TH1D("pt_hist", "pt_hist", 1000, 0, 100000);
    TH1D *cut_pt_hist = new TH1D("cut_pt_hist", "cut_pt_hist", 1000, 0, 100000);
    TH1D *deltaR_hist = new TH1D("deltaR_hist", "deltaR_hist", 100, 0, 5);
    TH1D *cut_deltaR_hist = new TH1D("cut_deltaR_hist", "cut_deltaR_hist", 100, 0, 5);
    TH1D *deltaPhi_hist = new TH1D("deltaPhi_hist", "deltaPhi_hist", 100, 0, 5);
    TH1D *cut_deltaPhi_hist = new TH1D("cut_deltaPhi_hist", "cut_deltaPhi_hist", 100, 0, 5);
    TH1D *tag_probe_deltaR_hist = new TH1D("tag_probe_deltaR_hist", "tag_probe_deltaR_hist", 100, 0, 1);
    TH1D *hlt_deltaR_hist = new TH1D("hlt_deltaR_hist", "hlt_deltaR_hist", 100, 0, 3);
    TH2D *hlt_deltaR_pt_hist = new TH2D("hlt_deltaR_pt_hist", "hlt_deltaR_pt_hist", 100, 0, 100000, 100, 0, 0.01);
    TH1D *tag_muon_momentum_hist = new TH1D("tag_muon_momentum_hist", "tag_muon_momentum_hist", 1000, 0, 100);
    TH1D *tag_muon_momentum_cut_hist = new TH1D("tag_muon_momentum_cut_hist", "tag_muon_momentum_cut_hist", 1000, 0, 100);
    TH1D *probe_muon_momentum_hist = new TH1D("probe_muon_momentum_hist", "probe_muon_momentum_hist", 1000, 0, 100);
    TH1D *probe_muon_eta_hist = new TH1D("probe_muon_eta_hist", "probe_muon_eta_hist", 100, -4, 4);
    TH1D *probe_muon_phi_hist = new TH1D("probe_muon_phi_hist", "probe_muon_phi_hist", 100, 0, 5);
    TH1D *HLT_probe_muon_momentum_cut_hist = new TH1D("HLT_probe_muon_momentum_cut_hist", "HLT_probe_muon_momentum_cut_hist", 1000, 0, 100);
    TH1D *L1_probe_muon_momentum_cut_hist = new TH1D("L1_probe_muon_momentum_cut_hist", "L1_probe_muon_momentum_cut_hist", 1000, 0, 100);
    TH1D *probe_muon_endcap_momentum_hist = new TH1D("probe_muon_endcap_momentum_hist", "probe_muon_endcap_momentum_hist", 1000, 0, 100);
    TH1D *probe_muon_endcap_deltaEta_hist = new TH1D("probe_muon_endcap_deltaEta_hist", "probe_muon_endcap_deltaEta_hist", 100, -4, 4);
    TH1D *probe_muon_endcap_deltaPhi_hist = new TH1D("probe_muon_endcap_deltaPhi_hist", "probe_muon_endcap_deltaPhi_hist", 100, 0, 4);
    TH1D *probe_muon_endcap_deltaR_hist = new TH1D("probe_muon_endcap_deltaR_hist", "probe_muon_endcap_deltaR_hist", 100, 0, 5);
    TH1D *probe_muon_endcap_momentum_cut_hist = new TH1D("probe_muon_endcap_momentum_cut_hist", "probe_muon_endcap_momentum_cut_hist", 1000, 0, 100);
    TH1D *probe_muon_endcap_momentum_thr1_cut_hist = new TH1D("probe_muon_endcap_momentum_thr1_cut_hist", "probe_muon_endcap_momentum_thr1_cut_hist", 1000, 0, 100);
    TH1D *probe_muon_endcap_momentum_thr2_cut_hist = new TH1D("probe_muon_endcap_momentum_thr2_cut_hist", "probe_muon_endcap_momentum_thr2_cut_hist", 1000, 0, 100);
    TH1D *probe_muon_endcap_momentum_thr3_cut_hist = new TH1D("probe_muon_endcap_momentum_thr3_cut_hist", "probe_muon_endcap_momentum_thr3_cut_hist", 1000, 0, 100);
    TH1D *probe_muon_endcap_momentum_thr4_cut_hist = new TH1D("probe_muon_endcap_momentum_thr4_cut_hist", "probe_muon_endcap_momentum_thr4_cut_hist", 1000, 0, 100);
    TH1D *probe_muon_endcap_momentum_thr5_cut_hist = new TH1D("probe_muon_endcap_momentum_thr5_cut_hist", "probe_muon_endcap_momentum_thr5_cut_hist", 1000, 0, 100);
    TH1D *probe_muon_endcap_momentum_thr6_cut_hist = new TH1D("probe_muon_endcap_momentum_thr6_cut_hist", "probe_muon_endcap_momentum_thr6_cut_hist", 1000, 0, 100);
    TH1D *probe_muon_barrel_momentum_thr1_cut_hist = new TH1D("probe_muon_barrel_momentum_thr1_cut_hist", "probe_muon_barrel_momentum_thr1_cut_hist", 1000, 0, 100);
    TH1D *probe_muon_barrel_momentum_thr2_cut_hist = new TH1D("probe_muon_barrel_momentum_thr2_cut_hist", "probe_muon_barrel_momentum_thr2_cut_hist", 1000, 0, 100);
    TH1D *probe_muon_barrel_momentum_thr3_cut_hist = new TH1D("probe_muon_barrel_momentum_thr3_cut_hist", "probe_muon_barrel_momentum_thr3_cut_hist", 1000, 0, 100);
    TH1D *probe_muon_barrel_momentum_thr4_cut_hist = new TH1D("probe_muon_barrel_momentum_thr4_cut_hist", "probe_muon_barrel_momentum_thr4_cut_hist", 1000, 0, 100);
    TH1D *probe_muon_barrel_momentum_thr5_cut_hist = new TH1D("probe_muon_barrel_momentum_thr5_cut_hist", "probe_muon_barrel_momentum_thr5_cut_hist", 1000, 0, 100);
    TH1D *probe_muon_barrel_momentum_thr6_cut_hist = new TH1D("probe_muon_barrel_momentum_thr6_cut_hist", "probe_muon_barrel_momentum_thr6_cut_hist", 1000, 0, 100);
    TH1D *probe_muon_endcap_deltaEta_cut_hist = new TH1D("probe_muon_endcap_deltaEta_cut_hist", "probe_muon_endcap_deltaEta_cut_hist", 100, -4, 4);
    TH1D *probe_muon_endcap_deltaPhi_cut_hist = new TH1D("probe_muon_endcap_deltaPhi_cut_hist", "probe_muon_endcap_deltaPhi_cut_hist", 100, 0, 4);
    TH1D *probe_muon_endcap_deltaR_cut_hist = new TH1D("probe_muon_endcap_deltaR_cut_hist", "probe_muon_endcap_deltaR_cut_hist", 100, 0, 5);
    TH1D *probe_muon_deltaR_thr1_hist = new TH1D("probe_muon_deltaR_thr1_hist", "probe_muon_deltaR_thr1_hist", 100, 0, 0.5);
    TH1D *probe_muon_deltaR_thr2_hist = new TH1D("probe_muon_deltaR_thr2_hist", "probe_muon_deltaR_thr2_hist", 100, 0, 0.5);
    TH1D *probe_muon_deltaR_thr3_hist = new TH1D("probe_muon_deltaR_thr3_hist", "probe_muon_deltaR_thr3_hist", 100, 0, 0.5);
    TH1D *probe_muon_deltaR_thr4_hist = new TH1D("probe_muon_deltaR_thr4_hist", "probe_muon_deltaR_thr4_hist", 100, 0, 0.5);
    TH1D *probe_muon_deltaR_thr5_hist = new TH1D("probe_muon_deltaR_thr5_hist", "probe_muon_deltaR_thr5_hist", 100, 0, 0.5);
    TH1D *probe_muon_deltaR_thr6_hist = new TH1D("probe_muon_deltaR_thr6_hist", "probe_muon_deltaR_thr6_hist", 100, 0, 0.5);
    TH1D *probe_muon_barrel_momentum_hist = new TH1D("probe_muon_barrel_momentum_hist", "probe_muon_barrel_momentum_hist", 1000, 0, 100);
    TH1D *probe_muon_barrel_deltaEta_hist = new TH1D("probe_muon_barrel_deltaEta_hist", "probe_muon_barrel_deltaEta_hist", 100, -4, 4);
    TH1D *probe_muon_barrel_deltaPhi_hist = new TH1D("probe_muon_barrel_deltaPhi_hist", "probe_muon_barrel_deltaPhi_hist", 100, 0, 4);
    TH1D *probe_muon_barrel_deltaR_hist = new TH1D("probe_muon_barrel_deltaR_hist", "probe_muon_barrel_deltaR_hist", 100, 0, 100000);
    TH1D *probe_muon_barrel_momentum_cut_hist = new TH1D("probe_muon_barrel_momentum_cut_hist", "probe_muon_barrel_momentum_cut_hist", 1000, 0, 100);
    TH1D *probe_muon_barrel_deltaEta_cut_hist = new TH1D("probe_muon_barrel_deltaEta_cut_hist", "probe_muon_barrel_deltaEta_cut_hist", 100, -4, 4);
    TH1D *probe_muon_barrel_deltaPhi_cut_hist = new TH1D("probe_muon_barrel_deltaPhi_cut_hist", "probe_muon_barrel_deltaPhi_cut_hist", 100, 0, 4);
    TH1D *probe_muon_barrel_deltaR_cut_hist = new TH1D("probe_muon_barrel_deltaR_cut_hist", "probe_muon_barrel_deltaR_cut_hist", 100, 0, 5);
    TH1D *L1_probe_DeltaR_hist = new TH1D("L1_probe_DeltaR_hist", "L1_probe_DeltaR_hist", 100, 0, 4);
    TH1D *probe_thrNum1_pt = new TH1D("probe_thrNum1_pt", "probe_thrNum1_pt", 1000, 0, 100);
    TH1D *probe_thrNum2_pt = new TH1D("probe_thrNum2_pt", "probe_thrNum2_pt", 1000, 0, 100);
    TH1D *probe_thrNum3_pt = new TH1D("probe_thrNum3_pt", "probe_thrNum3_pt", 1000, 0, 100);
    TH1D *probe_thrNum4_pt = new TH1D("probe_thrNum4_pt", "probe_thrNum4_pt", 1000, 0, 100);
    TH1D *probe_thrNum5_pt = new TH1D("probe_thrNum5_pt", "probe_thrNum5_pt", 1000, 0, 100);
    TH1D *probe_thrNum6_pt = new TH1D("probe_thrNum6_pt", "probe_thrNum6_pt", 1000, 0, 100);
    TH2D *L1req_pt_hist = new TH2D("L1req_pt_hist", "L1req_pt_hist", 100, 0, 20, 100, 0, 1);
    TH2D *L1probe_deltaR_pt_hist = new TH2D("L1probe_deltaR_pt_hist", "L1probe_deltaR_pt_hist", 100, 0, 20, 100, 0, 1);

    //int entry = 10;
    int entry = chain->GetEntries();
    std::cout << entry << std::endl;

    std::cout << "start selecting muon..." << std::endl;

    //std::vector<std::vector<std::pair<int, int>>> tag_and_probe_mu_pair;

    int counts = 0;

    std::ofstream deltaR_text;
    std::string filename ="img1115/deltaR.txt";
    deltaR_text.open(filename, std::ios::out);

    TFile hist_file("img1115/hist1115.root", "RECREATE");

    for(int i = 0; i < entry; i++){
        //counts++;
        //std::cout << "counts = " << counts << std::endl;
        
        std::vector<std::pair<int, int>> mu_pair_number;
        
        chain->GetEntry(i);
        
        int trig_chain = 0;
        bool flag_eventselection = 0;

        for(int j = 0; j < trigger_info_chain->size(); j++){
            // for j_psi
            //if(trigger_info_chain->at(j) == "HLT_mu20_2mu0noL1_JpsimumuFS"){
            // for j/psi in Run3
            if(trigger_info_chain->at(j) == "HLT_mu24_ivarmedium_mu4_probe_L1MU14FCH"){
            // z mumu
            //if(trigger_info_chain->at(j) == "HLT_mu26_ivarmedium"){
                //if(trigger_info_isPassed->at(j) == 1){
                    trig_chain = j;
                    flag_eventselection = 1;
                    // std::cout << "there are j psi" << std::endl;
                //}
            }
        }
				
        //if( vxp_chi2->at(0) > 20. ) continue;

        //if(flag_eventselection){
        for(int j = 0; j < mu_m->size(); j++){
            //mu select
            bool flag_mu1_author = mu_author->at(j) == 1;
            bool flag_mu1_Type = mu_muonType->at(j) == 0;
            bool flag_mu1_eta = fabsf(mu_eta->at(j)) < 2.5;

            if(flag_mu1_author && flag_mu1_Type && flag_mu1_eta){
                TLorentzVector mu1;
                mu1.SetPtEtaPhiM(mu_pt->at(j), mu_eta->at(j), mu_phi->at(j), mu_m->at(j));
                for(int k = 0; k < mu_m->size(); k++){
                    if (j == k) continue;
                    //mu select & charge cut
                    int charge = mu_charge->at(j) * mu_charge->at(k);
                    bool flag_mu2_author = mu_author->at(k) == 1;
                    bool flag_mu2_Type = mu_muonType->at(k) == 0;
                    bool flag_charge = charge == -1;
                    bool flag_mu2_eta = fabsf(mu_eta->at(k)) < 2.5;
                    
                    if (flag_mu2_author && flag_mu2_Type && flag_charge && flag_mu2_eta){
                        TLorentzVector mu2;
                        mu2.SetPtEtaPhiM(mu_pt->at(k), mu_eta->at(k), mu_phi->at(k), mu_m->at(k));

                        TLorentzVector mu_pair = mu1 + mu2;
                        float pair_mass = mu_pair.M();
                        float pair_deltaR = mu2.DeltaR(mu1);
                        float pair_deltaPhi = mu2.DeltaPhi(mu1);

                        mass_hist->Fill(pair_mass);
                        pt_hist->Fill(mu_pt->at(k));
                        deltaR_hist->Fill(pair_deltaR);
                        deltaPhi_hist->Fill(pair_deltaPhi);

                        // deltaR cut & mass cut & deltaPhi cut
                        bool flag_pair_DeltaR = pair_deltaR > 0.25;
                        //jpsi mass
                        bool flag_pair_mass = pair_mass > 2700 && 3500 > pair_mass;
                        //bool flag_pair_mass = pair_mass > 80000 && 100000 > pair_mass;
                        //bool flag_pair_deltaPhi = fabsf(pair_deltaPhi) > 0.14 && 3.0 >fabsf(pair_deltaPhi);
                        bool flag_pair_deltaPhi = 3.0 >fabsf(pair_deltaPhi);

                        //deltaR cut
                        // req_L1_tag_deltaR
                        float deltaR_req_tag = 0;
                        if (mu_pt->at(j) > 10000){
                            deltaR_req_tag = 0.08;
                        }
                        else{
                            deltaR_req_tag = ( mu_pt->at(j) / 1000 ) * ( -0.01 ) + 0.18;
                        }
                        // req_L1_probe_deltaR
                        float deltaR_req_probe = 0;
                        if (mu_pt->at(k) > 10000){
                            deltaR_req_probe = 0.08;
                        }
                        else{
                            deltaR_req_probe = ( mu_pt->at(k) / 1000 ) * ( -0.01 ) + 0.18;
                        }

                        bool flag_deltaR_mu_L1 = ( deltaR_req_tag + deltaR_req_probe ) < pair_deltaR;


                        if (flag_pair_DeltaR && flag_pair_mass && flag_pair_deltaPhi && flag_deltaR_mu_L1){
                            mu_pair_number.push_back(std::make_pair(j, k));
                            cut_pair_mass_hist->Fill(pair_mass);
                            cut_pt_hist->Fill(mu_pt->at(k));
                            cut_deltaR_hist->Fill(pair_deltaR);
                            cut_deltaPhi_hist->Fill(pair_deltaPhi);
                        }
                    }
                }
            }
        }
        //    }
        
        
        //tag
        //std::cout << "start tagging" << std::endl;

				
        for (int j = 0; j < mu_pair_number.size(); j++){
            std::cout << "mu_pair_size= " << mu_pair_number.size() << std::endl;
            std::cout << "mu pair number = " << j << std::endl;

            int tag_muon_number = mu_pair_number.at(j).first;
            int probe_muon_number = mu_pair_number.at(j).second;

            if (mu_ext_b_targetEtaVec->at(probe_muon_number).size() == 0) continue;
            if (mu_ext_b_targetPhiVec->at(probe_muon_number).size() == 0) continue;
            if (trigger_info_ptVec->size() == 0) continue;

            float tag_muon_pt = mu_pt->at(tag_muon_number) / 1000;
            float probe_muon_pt = mu_pt->at(probe_muon_number) / 1000;

            std::cout << "trigger info size = " << trigger_info_ptVec->at(trig_chain).size() << std::endl;
            
            TVector3 tag_muon;
            tag_muon.SetPtEtaPhi(mu_pt->at(tag_muon_number), mu_eta->at(tag_muon_number), mu_phi->at(tag_muon_number));

            TVector3 probe_muon;
            probe_muon.SetPtEtaPhi(mu_pt->at(probe_muon_number), mu_eta->at(probe_muon_number), mu_phi->at(probe_muon_number));

            float tag_probe_deltaR = tag_muon.DeltaR(probe_muon);
            
            //tag matching with HLT
            for (int k = 0; k < trigger_info_ptVec->at(trig_chain).size(); k++){
                std::cout << "trig info number = " << k << std::endl;
                TVector3 hlt_muon;
                hlt_muon.SetPtEtaPhi(trigger_info_ptVec->at(trig_chain).at(k), trigger_info_etaVec->at(trig_chain).at(k), trigger_info_phiVec->at(trig_chain).at(k));

                float tag_DeltaR = tag_muon.DeltaR(hlt_muon);
                hlt_deltaR_hist->Fill(tag_DeltaR);
                hlt_deltaR_pt_hist->Fill(tag_muon_pt, tag_DeltaR);
                tag_muon_momentum_hist->Fill(tag_muon_pt);

                //pt cut
                bool flag_mu_tag_pt = tag_muon_pt > 24;

                bool flag_trigger_menu = trigger_info_chain->at(trig_chain) == "HLT_mu24_ivarmedium_mu4_probe_L1MU14FCH";

                //probe matching with HLT
                for (int l = 0; l < trigger_info_ptVec->at(trig_chain).size(); l++){
                    std::cout << "probe trig info number = " << l << std::endl;
                    TVector3 hlt_muon_for_probe;
                    hlt_muon_for_probe.SetPtEtaPhi(trigger_info_ptVec->at(trig_chain).at(l), trigger_info_etaVec->at(trig_chain).at(l), trigger_info_phiVec->at(trig_chain).at(l));

                    float probe_DeltaR = probe_muon.DeltaR(hlt_muon_for_probe);

                    bool flag_trigger_menu_for_probe = trigger_info_chain->at(trig_chain) == "HLT_mu24_ivarmedium_mu4_probe_L1MU14FCH";

                    if(tag_DeltaR < 0.01 && probe_DeltaR < 0.01 && flag_trigger_menu && flag_mu_tag_pt && flag_trigger_menu_for_probe){
                        tag_muon_momentum_cut_hist->Fill(tag_muon_pt);
                        probe_muon_momentum_hist->Fill(probe_muon_pt);
                        probe_muon_eta_hist->Fill(mu_eta->at(probe_muon_number));
                        probe_muon_phi_hist->Fill(mu_phi->at(probe_muon_number));
                        tag_probe_deltaR_hist->Fill(tag_probe_deltaR);
/*
                        //endcap
                        if( fabsf(mu_eta->at(probe_muon_number)) > 1.05 && fabsf(mu_eta->at(probe_muon_number)) < 2.4 ){
                            probe_muon_endcap_momentum_hist->Fill(probe_muon_pt);
                            probe_muon_endcap_deltaR_hist->Fill(tag_probe_deltaR);
                            probe_muon_endcap_deltaEta_hist->Fill(mu_eta->at(probe_muon_number));
                            probe_muon_endcap_deltaPhi_hist->Fill(mu_phi->at(probe_muon_number));
                        }
    */
                        // barrel
                        if( fabsf(mu_eta->at(probe_muon_number)) < 1.05 ){
                            probe_muon_barrel_momentum_hist->Fill(probe_muon_pt);
                            probe_muon_barrel_deltaR_hist->Fill(tag_probe_deltaR);
                            probe_muon_barrel_deltaEta_hist->Fill(mu_eta->at(probe_muon_number));
                            probe_muon_barrel_deltaPhi_hist->Fill(mu_phi->at(probe_muon_number));
                        }
                        //L1 probe
                        for (int m = 0; m < trig_L1_mu_eta->size(); m++){
                            std::cout << "proving..." << std::endl;
                            //RoI Number
                            if (trig_L1_mu_RoINumber->at(m) == -1) continue;

                            //tag and probe to only RPC
                            if (fabsf(mu_eta->at(probe_muon_number)) > 1.05 && fabsf(mu_eta->at(probe_muon_number)) < 2.4) continue;
                            
                            //threshold Number
                            int thr_Num = trig_L1_mu_thrNumber->at(m);

                            std::cout << mu_ext_b_targetEtaVec->at(probe_muon_number).at(2) << std::endl;

                            // L1 efficiency
                            float L1_probe_deltaEta = mu_ext_b_targetEtaVec->at(probe_muon_number).at(2) - trig_L1_mu_eta->at(m);
                            float L1_probe_deltaPhi = TVector2::Phi_mpi_pi(mu_ext_b_targetPhiVec->at(probe_muon_number).at(2) - trig_L1_mu_phi->at(m));

                            float L1_probe_DeltaR = TMath::Sqrt(L1_probe_deltaEta*L1_probe_deltaEta + L1_probe_deltaPhi*L1_probe_deltaPhi);

                            // req_L1_deltaR
                            float deltaR_L1_req_mu = 0;
                            
                            if (probe_muon_pt > 10){
                                deltaR_L1_req_mu = 0.08;
                            }
                            else{
                                deltaR_L1_req_mu = ( probe_muon_pt ) * ( -0.01 ) + 0.18;
                            }
                            L1_probe_DeltaR_hist->Fill(L1_probe_DeltaR);
                            deltaR_text << deltaR_L1_req_mu << " " << m << " " << tag_muon_number << " " <<  probe_muon_number << " " << L1_probe_DeltaR << std::endl;

                            L1req_pt_hist->Fill(probe_muon_pt, deltaR_L1_req_mu);
                            L1probe_deltaR_pt_hist->Fill(probe_muon_pt, L1_probe_DeltaR);

                            deltaR_L1_req_mu = 0.1;

                            // L1 matching
                            if(L1_probe_DeltaR < deltaR_L1_req_mu){
                                L1_probe_muon_momentum_cut_hist->Fill(probe_muon_pt);

                                if (thr_Num >= 1) {
                                    probe_thrNum1_pt->Fill(probe_muon_pt);
                                    probe_muon_deltaR_thr1_hist->Fill(L1_probe_DeltaR);
                                }
                                if (thr_Num >= 2) {
                                    probe_thrNum2_pt->Fill(probe_muon_pt);
                                    probe_muon_deltaR_thr2_hist->Fill(L1_probe_DeltaR);
                                }
                                if (thr_Num >= 3) {
                                    probe_thrNum3_pt->Fill(probe_muon_pt);
                                    probe_muon_deltaR_thr3_hist->Fill(L1_probe_DeltaR);
                                }
                                if (thr_Num >= 4) {
                                    probe_thrNum4_pt->Fill(probe_muon_pt);
                                    probe_muon_deltaR_thr4_hist->Fill(L1_probe_DeltaR);
                                }
                                if (thr_Num >= 5) {
                                    probe_thrNum5_pt->Fill(probe_muon_pt);
                                    probe_muon_deltaR_thr5_hist->Fill(L1_probe_DeltaR);
                                }
                                if (thr_Num >= 6) {
                                    probe_thrNum6_pt->Fill(probe_muon_pt);
                                    probe_muon_deltaR_thr6_hist->Fill(L1_probe_DeltaR);
                                }
/*
                                //endcap
                                if( fabsf(mu_eta->at(probe_muon_number)) > 1.05 && fabsf(mu_eta->at(probe_muon_number)) < 2.4){
                                    probe_muon_endcap_momentum_cut_hist->Fill(probe_muon_pt);
                                    probe_muon_endcap_deltaR_cut_hist->Fill(tag_probe_deltaR);
                                    probe_muon_endcap_deltaEta_cut_hist->Fill(mu_eta->at(probe_muon_number));
                                    probe_muon_endcap_deltaPhi_cut_hist->Fill(mu_phi->at(probe_muon_number));
                                    if (thr_Num >= 1) {
                                        probe_muon_endcap_momentum_thr1_cut_hist->Fill(probe_muon_pt);
                                    }
                                    if (thr_Num >= 2) {
                                        probe_muon_endcap_momentum_thr2_cut_hist->Fill(probe_muon_pt);
                                    }
                                    if (thr_Num >= 3) {
                                        probe_muon_endcap_momentum_thr3_cut_hist->Fill(probe_muon_pt);
                                    }
                                    if (thr_Num >= 4) {
                                        probe_muon_endcap_momentum_thr4_cut_hist->Fill(probe_muon_pt);
                                    }
                                    if (thr_Num >= 5) {
                                        probe_muon_endcap_momentum_thr5_cut_hist->Fill(probe_muon_pt);
                                    }
                                    if (thr_Num >= 6) {
                                        probe_muon_endcap_momentum_thr6_cut_hist->Fill(probe_muon_pt);
                                    }
                                }
*/
                                // barrel
                                if( fabsf(mu_eta->at(probe_muon_number)) < 1.05){
                                    probe_muon_barrel_momentum_cut_hist->Fill(probe_muon_pt);
                                    probe_muon_barrel_deltaR_cut_hist->Fill(tag_probe_deltaR);
                                    probe_muon_barrel_deltaEta_cut_hist->Fill(tag_probe_deltaR);
                                    probe_muon_barrel_deltaPhi_cut_hist->Fill(mu_phi->at(probe_muon_number));
                                    if (thr_Num >= 1) {
                                        probe_muon_barrel_momentum_thr1_cut_hist->Fill(probe_muon_pt);
                                    }
                                    if (thr_Num >= 2) {
                                        probe_muon_barrel_momentum_thr2_cut_hist->Fill(probe_muon_pt);
                                    }
                                    if (thr_Num >= 3) {
                                        probe_muon_barrel_momentum_thr3_cut_hist->Fill(probe_muon_pt);
                                    }
                                    if (thr_Num >= 4) {
                                        probe_muon_barrel_momentum_thr4_cut_hist->Fill(probe_muon_pt);
                                    }
                                    if (thr_Num >= 5) {
                                        probe_muon_barrel_momentum_thr5_cut_hist->Fill(probe_muon_pt);
                                    }
                                    if (thr_Num >= 6) {
                                        probe_muon_barrel_momentum_thr6_cut_hist->Fill(probe_muon_pt);
                                    }
                                }
                            }
                            break;
                        }
                    }
                }
            }
        }
    }
    deltaR_text.close();

    hist_file.cd();

    mass_hist->Write();
    deltaR_hist->Write();
    deltaPhi_hist->Write();
    cut_pair_mass_hist->Write();
    pt_hist->Write();
    cut_pt_hist->Write();
    cut_deltaR_hist->Write();
    cut_deltaPhi_hist->Write();
    tag_probe_deltaR_hist->Write();
    tag_muon_momentum_cut_hist->Write();
    probe_muon_momentum_hist->Write();
    L1_probe_muon_momentum_cut_hist->Write();
    probe_muon_barrel_momentum_hist->Write();
    probe_muon_barrel_momentum_cut_hist->Write();
    //probe_muon_endcap_momentum_hist->Write();
    //probe_muon_endcap_momentum_cut_hist->Write();
    hlt_deltaR_hist->Write();
    hlt_deltaR_pt_hist->Write();
    //probe_muon_endcap_deltaEta_hist->Write();
    //probe_muon_endcap_deltaPhi_hist->Write();
    probe_muon_barrel_deltaEta_hist->Write();
    //probe_muon_endcap_deltaPhi_hist->Write();
    //probe_muon_endcap_deltaEta_hist->Write();
    //probe_muon_endcap_deltaEta_cut_hist->Write();
    //probe_muon_endcap_deltaPhi_hist->Write();
    //probe_muon_endcap_deltaPhi_cut_hist->Write();
    probe_muon_barrel_deltaEta_hist->Write();
    probe_muon_barrel_deltaEta_cut_hist->Write();
    probe_muon_barrel_deltaPhi_hist->Write();
    probe_muon_barrel_deltaPhi_cut_hist->Write();
    L1_probe_DeltaR_hist->Write();
    probe_muon_eta_hist->Write();
    probe_muon_phi_hist->Write();
    //probe_muon_endcap_deltaR_hist->Write();
    //probe_muon_endcap_deltaR_cut_hist->Write();
    probe_muon_barrel_deltaR_hist->Write();
    probe_muon_barrel_deltaR_cut_hist->Write();
    probe_thrNum1_pt->Write();
    probe_thrNum2_pt->Write();
    probe_thrNum3_pt->Write();
    probe_thrNum4_pt->Write();
    probe_thrNum5_pt->Write();
    probe_thrNum6_pt->Write();
    //probe_muon_endcap_momentum_thr1_cut_hist->Write();
    //probe_muon_endcap_momentum_thr2_cut_hist->Write();
    //probe_muon_endcap_momentum_thr3_cut_hist->Write();
    //probe_muon_endcap_momentum_thr4_cut_hist->Write();
    //probe_muon_endcap_momentum_thr5_cut_hist->Write();
    //probe_muon_endcap_momentum_thr6_cut_hist->Write();
    probe_muon_deltaR_thr1_hist->Write();
    probe_muon_deltaR_thr2_hist->Write();
    probe_muon_deltaR_thr3_hist->Write();
    probe_muon_deltaR_thr4_hist->Write();
    probe_muon_deltaR_thr5_hist->Write();
    probe_muon_deltaR_thr6_hist->Write();
    probe_muon_barrel_momentum_thr1_cut_hist->Write();
    probe_muon_barrel_momentum_thr2_cut_hist->Write();
    probe_muon_barrel_momentum_thr3_cut_hist->Write();
    probe_muon_barrel_momentum_thr4_cut_hist->Write();
    probe_muon_barrel_momentum_thr5_cut_hist->Write();
    probe_muon_barrel_momentum_thr6_cut_hist->Write();
    L1req_pt_hist->Write();
    L1probe_deltaR_pt_hist->Write();
    hist_file.Close();
    std::cout << "end succeeded" << std::endl;
}
