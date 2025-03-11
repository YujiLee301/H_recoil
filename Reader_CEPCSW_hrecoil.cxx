#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "HggTwoSidedCBPdf.h"
#include "HggTwoSidedCBPdf.cxx"

#include "TLorentzVector.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "EvtProcessor.h"
#include "CEPCStyle.h"

#include <TCanvas.h>
#include <TH1D.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooPlot.h>
#include <TLatex.h>

float upmass=130;
float downmass=120;

void FitAndPlot(TH1D* hist, TString filename) {
    // 创建画布
    TCanvas* c = new TCanvas("c", "Gaussian Fit", 800, 600);
    
    // 设置直方图样式
    hist->SetLineColor(kBlack);
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(0.8);
    
    // 执行高斯拟合
    hist->Fit("gaus", "Q"); // "Q" 表示静默模式（不打印拟合结果）
    
    // 获取拟合函数和参数
    TF1* fitFunc = hist->GetFunction("gaus");
    double mean = fitFunc->GetParameter(1);
    double sigma = fitFunc->GetParameter(2);
    
    // 绘制直方图和拟合曲线
    hist->Draw("E"); // "E" 绘制误差条
    fitFunc->Draw("same");
    
    // 添加统计信息标签
    TLatex latex;
    latex.SetNDC(true); // 使用归一化坐标系（与画布尺寸无关）
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.65, 0.8, Form("#mu = %.2f", mean));
    latex.DrawLatex(0.65, 0.75, Form("#sigma = %.2f", sigma));
    
    // 更新画布
    c->Update();
    c->Draw();
    c->SaveAs(filename);
}

void RooFitGaussian(RooDataSet data, TString filename) {
    SetCEPCStyle();
    TCanvas* c = new TCanvas("c_roofit", "RooFit Gaussian Fit", 800, 600);

    RooRealVar x("x", "invMass", 100, 160);


    
    RooRealVar meanCB("meanCB","",125.,120.,130.);
    RooRealVar sigmaCB("sigmaCB","",0.5,0.,5.);
    RooRealVar alphaCBLo("alphaCBLo","",1.,0.,5.);
    RooRealVar alphaCBHi("alphaCBHi","",1.,0.,10.);
    RooRealVar nCBLo("nCBLo","",2.,0.,5.);
    RooRealVar nCBHi("nCBHi","",2.,0.,10.);

    RooAbsPdf *pdf =new HggTwoSidedCBPdf("pdf", "pdf", x, meanCB, sigmaCB, alphaCBLo, nCBLo, alphaCBHi, nCBHi);
    // RooAbsPdf *pdf =new RooGaussian("pdf", "pdf", invMass, meanCB, sigmaCB);

    RooFitResult *result=pdf->fitTo(data);

    RooPlot *frame=x.frame(50);
    
    
    data.plotOn(frame);
    pdf->plotOn(frame);
    frame->GetYaxis()->SetTitle("Events");
    frame->GetXaxis()->SetTitle("m_{#mu#mu}");
    frame->Draw();
    TLatex *tex = new TLatex();
    FormatLatex(tex);
    tex->DrawLatexNDC(0.75, 0.86, "#it{#bf{CEPC}} Ref-TDR");
    tex->DrawLatexNDC(0.75, 0.8, "ZH#rightarrow #nu#nu#mu#mu");
    // tex->DrawLatexNDC(0.75, 0.66, (TString)"Output: "+ (TString)to_string(evt_h1)+(TString)" events");
    tex->DrawLatexNDC(0.17,0.58,Form("mean = %.3f #pm %.3f",meanCB.getVal(),meanCB.getError()));
    tex->DrawLatexNDC(0.17,0.51,Form("#sigma = %.3f #pm %.3f",sigmaCB.getVal(),sigmaCB.getError()));
    tex->DrawLatexNDC(0.17,0.44,Form("Resolution = %.2f%%",100*sigmaCB.getVal()/meanCB.getVal()));
    tex->DrawLatexNDC(0.17,0.37,Form("#chi^{2}/ndf = %.2f",frame->chiSquare(6)));
    c->SaveAs(filename);
}

void processing(TString outputTag);         // For Float variable
// void fun2(TString a, int b); // For int variable
void Reader_CEPCSW_hrecoil_zz_bkg()
{   
    gROOT->ProcessLine(".L EvtProcessor.cxx");
    SetCEPCStyle();
    processing("Hrecoil");

}
void processing(TString outputTag){
    TChain *Flowzz_l0mumu = new TChain("events", "events");
    TChain *Flow4mu = new TChain("events", "events");
    TChain *Flowsl0mu_down= new TChain("events", "events");
    TChain *Flowsl0mu_up = new TChain("events", "events");
    TChain *Flowzzorww_l0 = new TChain("events", "events");
    TChain *Flowll = new TChain("events", "events");
    TChain *Flowsig = new TChain("events", "events");
    TChain *Flowsznu_l0mumu = new TChain("events", "events");
    std::ifstream fileListzz_l0mumu("sample_input_zz_l0mumu.txt");
    std::ifstream fileList4mu("sample_input_zz_l04mu.txt");
    std::ifstream fileListsl0mu_down("sample_input_zz_sl0mu_down.txt");
    std::ifstream fileListsl0mu_up("sample_input_zz_sl0mu_up.txt");
    std::ifstream fileListzzorww_l0("sample_input_zzorww_l0mumu.txt");
    std::ifstream fileListll("sample_input_ll.txt");
    std::ifstream fileListsig("sample_input_sig.txt");
    std::ifstream fileListsznu_l0mumu("sample_input_sznu_l0mumu.txt");
    float Lumi = 20000; //fb-1

    float xszz_l0mumu=19.38;
    float xssznu_l0mumu=43.42;
    float xs4mu=15.56;
    float xssl0mu_down=136.14;
    float xssl0mu_up=87.39;
    float xszzorww_l0=221.10;
    float xsll=5332.71;
    float xssig=6.77;
    float NUMzz_l0mumu=0; float DENzz_l0mumu=0;
    float NUM4mu=0; float DEN4mu=0;
    float NUMsl0mu_down=0; float DENsl0mu_down=0;
    float NUMsl0mu_up=0; float DENsl0mu_up=0;
    float NUMzzorww_l0=0; float DENzzorww_l0=0;
    float NUMll=0; float DENll=0;
    float NUMsig=0; float DENsig=0;
    float NUMsznu_l0mumu=0; float DENsznu_l0mumu=0;

    std::string linezz_l0mumu;
    std::string line4mu;
    std::string linesl0mu_down;
    std::string linesl0mu_up;
    std::string linezzorww_l0;
    std::string linell;
    std::string linesig;
    std::string linesznu_l0mumu;
    while (std::getline(fileListzz_l0mumu, linezz_l0mumu)) {
        size_t startzz_l0mumu = linezz_l0mumu.find_first_not_of(" \t");
        size_t endzz_l0mumu= linezz_l0mumu.find_last_not_of(" \t");
        if (startzz_l0mumu == std::string::npos || linezz_l0mumu.empty()) {
            continue;
        }
        std::string fileNamezz_l0mumu = linezz_l0mumu.substr(startzz_l0mumu, endzz_l0mumu - startzz_l0mumu + 1);

        Flowzz_l0mumu->Add(fileNamezz_l0mumu.c_str());
    }

    fileListzz_l0mumu.close();
    TTreeReader R_Flowzz_l0mumu(Flowzz_l0mumu);
    TTreeReaderArray<Float_t> PFO_Ezz_l0mumu(R_Flowzz_l0mumu,  "CyberPFO.energy");
    TTreeReaderArray<Float_t> PFO_Czz_l0mumu(R_Flowzz_l0mumu,  "CyberPFO.charge");
    TTreeReaderArray<Float_t> PFO_xzz_l0mumu(R_Flowzz_l0mumu,  "CyberPFO.momentum.x");
    TTreeReaderArray<Float_t> PFO_yzz_l0mumu(R_Flowzz_l0mumu,  "CyberPFO.momentum.y");
    TTreeReaderArray<Float_t> PFO_zzz_l0mumu(R_Flowzz_l0mumu,  "CyberPFO.momentum.z");

    TTreeReaderArray<double> Truth_mzz_l0mumu(R_Flowzz_l0mumu,  "MCParticle.mass");
    TTreeReaderArray<int> Truth_PDGzz_l0mumu(R_Flowzz_l0mumu,   "MCParticle.PDG");
    TTreeReaderArray<int> Truth_statuszz_l0mumu(R_Flowzz_l0mumu, "MCParticle.generatorStatus");
    TTreeReaderArray<Float_t> Truth_xzz_l0mumu(R_Flowzz_l0mumu,  "MCParticle.momentum.x");
    TTreeReaderArray<Float_t> Truth_yzz_l0mumu(R_Flowzz_l0mumu,  "MCParticle.momentum.y");
    TTreeReaderArray<Float_t> Truth_zzz_l0mumu(R_Flowzz_l0mumu,  "MCParticle.momentum.z");
    RooRealVar invMass("invMass", "invMass", 120, 130);
    RooRealVar weight("weight", "Event Weight", 0, 1000000);
    RooDataSet* data = new RooDataSet("data","data",  RooArgSet(invMass,weight));
    //RooRealVar invMassZ ("invMassZ", "invMassZ", 80, 100);
    //RooDataSet dataZ("dataZ","dataZ", invMassZ);
    //TH1D* MH; TH1D* MZ;
    //MH = new TH1D("H_recoil_Mass","H_recoil_Mass",80,100,160);
    //MZ = new TH1D("Z_Mass","Z_Mass",80,0,400);
    std::vector<float> MH_zz_l0mumu;
    while (R_Flowzz_l0mumu.Next()){
        DENzz_l0mumu++;
        int PFOsizezz_l0mumu=PFO_Ezz_l0mumu.GetSize();
        int MCsizezz_l0mumu =Truth_mzz_l0mumu.GetSize();
        if (PFOsizezz_l0mumu==0) continue;
        std::vector<Float_t> indices_PFO(PFOsizezz_l0mumu);
        std::vector<Float_t> indices_Truth(MCsizezz_l0mumu);

        std::vector<Float_t> PFO_E_zz_l0mumu(PFOsizezz_l0mumu);
        std::vector<Float_t> PFO_C_zz_l0mumu(PFOsizezz_l0mumu);
        std::vector<Float_t> PFO_x_zz_l0mumu(PFOsizezz_l0mumu);
        std::vector<Float_t> PFO_y_zz_l0mumu(PFOsizezz_l0mumu);
        std::vector<Float_t> PFO_z_zz_l0mumu(PFOsizezz_l0mumu);

        std::vector<double> Truth_m_zz_l0mumu(MCsizezz_l0mumu);
        std::vector<int> Truth_PDG_zz_l0mumu(MCsizezz_l0mumu);
        std::vector<int> Truth_status_zz_l0mumu(MCsizezz_l0mumu);
        std::vector<Float_t> Truth_x_zz_l0mumu(MCsizezz_l0mumu);
        std::vector<Float_t> Truth_y_zz_l0mumu(MCsizezz_l0mumu);
        std::vector<Float_t> Truth_z_zz_l0mumu(MCsizezz_l0mumu);

        for (int kzz_l0mumu=0;kzz_l0mumu<PFOsizezz_l0mumu;kzz_l0mumu++){
            PFO_E_zz_l0mumu[kzz_l0mumu]=PFO_Ezz_l0mumu[kzz_l0mumu];
            PFO_C_zz_l0mumu[kzz_l0mumu]=PFO_Czz_l0mumu[kzz_l0mumu];
            PFO_x_zz_l0mumu[kzz_l0mumu]=PFO_xzz_l0mumu[kzz_l0mumu];
            PFO_y_zz_l0mumu[kzz_l0mumu]=PFO_yzz_l0mumu[kzz_l0mumu];
            PFO_z_zz_l0mumu[kzz_l0mumu]=PFO_zzz_l0mumu[kzz_l0mumu];
        }

        for (int jzz_l0mumu=0;jzz_l0mumu<MCsizezz_l0mumu;jzz_l0mumu++){
            Truth_m_zz_l0mumu[jzz_l0mumu]=Truth_mzz_l0mumu[jzz_l0mumu];
            Truth_PDG_zz_l0mumu[jzz_l0mumu]=Truth_PDGzz_l0mumu[jzz_l0mumu];
            Truth_status_zz_l0mumu[jzz_l0mumu]=Truth_statuszz_l0mumu[jzz_l0mumu];
            Truth_x_zz_l0mumu[jzz_l0mumu]=Truth_xzz_l0mumu[jzz_l0mumu];
            Truth_y_zz_l0mumu[jzz_l0mumu]=Truth_yzz_l0mumu[jzz_l0mumu];
            Truth_z_zz_l0mumu[jzz_l0mumu]=Truth_zzz_l0mumu[jzz_l0mumu];
        }


        EvtProcessor Processorzz_l0mumu("mm");
        Processorzz_l0mumu.initialization();
        Processorzz_l0mumu.setPFOinfo(PFO_E_zz_l0mumu,PFO_C_zz_l0mumu, PFO_x_zz_l0mumu, PFO_y_zz_l0mumu,PFO_z_zz_l0mumu,PFOsizezz_l0mumu);
        Processorzz_l0mumu.setMCtruthinfo(Truth_m_zz_l0mumu,Truth_PDG_zz_l0mumu,Truth_status_zz_l0mumu,Truth_x_zz_l0mumu,Truth_y_zz_l0mumu,Truth_z_zz_l0mumu,MCsizezz_l0mumu);
        Processorzz_l0mumu.GenLepSelection();
        Processorzz_l0mumu.PFOLepSelection();
        Processorzz_l0mumu.GoodLepSelection();
        Processorzz_l0mumu.findZleps();
        if(Processorzz_l0mumu.ZLep.size()<2) continue;
        Processorzz_l0mumu.GetZHMass();
        /*MH->Fill(Processor.rec_MH);
        MZ->Fill(Processor.rec_MZ);
        if ((Processorzz_l0mumu.rec_MZ<100)&&(Processorzz_l0mumu.rec_MZ>80)){
            invMassZ.setVal(Processor.rec_MZ);
            dataZ.add(invMassZ);
        }*/
        if(Processorzz_l0mumu.rec_MH>upmass) continue;
        if(Processorzz_l0mumu.rec_MH<downmass) continue;
        MH_zz_l0mumu.push_back(Processorzz_l0mumu.rec_MH);
        NUMzz_l0mumu++;
        

    }
    
    for (int p_zz_l0mumu=0; p_zz_l0mumu<NUMzz_l0mumu; p_zz_l0mumu++){
        invMass.setVal(MH_zz_l0mumu[p_zz_l0mumu]);
        weight.setVal(Lumi*xszz_l0mumu/DENzz_l0mumu);
        data->add(RooArgSet(invMass,weight));
    }
    while (std::getline(fileListsznu_l0mumu, linesznu_l0mumu)) {
        size_t startsznu_l0mumu = linesznu_l0mumu.find_first_not_of(" \t");
        size_t endsznu_l0mumu= linesznu_l0mumu.find_last_not_of(" \t");
        if (startsznu_l0mumu == std::string::npos || linesznu_l0mumu.empty()) {
            continue;
        }
        std::string fileNamesznu_l0mumu = linesznu_l0mumu.substr(startsznu_l0mumu, endsznu_l0mumu - startsznu_l0mumu + 1);

        Flowsznu_l0mumu->Add(fileNamesznu_l0mumu.c_str());
    }

    fileListsznu_l0mumu.close();
    TTreeReader R_Flowsznu_l0mumu(Flowsznu_l0mumu);
    TTreeReaderArray<Float_t> PFO_Esznu_l0mumu(R_Flowsznu_l0mumu,  "CyberPFO.energy");
    TTreeReaderArray<Float_t> PFO_Csznu_l0mumu(R_Flowsznu_l0mumu,  "CyberPFO.charge");
    TTreeReaderArray<Float_t> PFO_xsznu_l0mumu(R_Flowsznu_l0mumu,  "CyberPFO.momentum.x");
    TTreeReaderArray<Float_t> PFO_ysznu_l0mumu(R_Flowsznu_l0mumu,  "CyberPFO.momentum.y");
    TTreeReaderArray<Float_t> PFO_zsznu_l0mumu(R_Flowsznu_l0mumu,  "CyberPFO.momentum.z");

    TTreeReaderArray<double> Truth_msznu_l0mumu(R_Flowsznu_l0mumu,  "MCParticle.mass");
    TTreeReaderArray<int> Truth_PDGsznu_l0mumu(R_Flowsznu_l0mumu,   "MCParticle.PDG");
    TTreeReaderArray<int> Truth_statussznu_l0mumu(R_Flowsznu_l0mumu, "MCParticle.generatorStatus");
    TTreeReaderArray<Float_t> Truth_xsznu_l0mumu(R_Flowsznu_l0mumu,  "MCParticle.momentum.x");
    TTreeReaderArray<Float_t> Truth_ysznu_l0mumu(R_Flowsznu_l0mumu,  "MCParticle.momentum.y");
    TTreeReaderArray<Float_t> Truth_zsznu_l0mumu(R_Flowsznu_l0mumu,  "MCParticle.momentum.z");
    //RooRealVar invMassZ ("invMassZ", "invMassZ", 80, 100);
    //RooDataSet dataZ("dataZ","dataZ", invMassZ);
    //TH1D* MH; TH1D* MZ;
    //MH = new TH1D("H_recoil_Mass","H_recoil_Mass",80,100,160);
    //MZ = new TH1D("Z_Mass","Z_Mass",80,0,400);
    std::vector<float> MH_sznu_l0mumu;
    while (R_Flowsznu_l0mumu.Next()){
        DENsznu_l0mumu++;
        int PFOsizesznu_l0mumu=PFO_Esznu_l0mumu.GetSize();
        int MCsizesznu_l0mumu =Truth_msznu_l0mumu.GetSize();
        if (PFOsizesznu_l0mumu==0) continue;
        std::vector<Float_t> indices_PFO(PFOsizesznu_l0mumu);
        std::vector<Float_t> indices_Truth(MCsizesznu_l0mumu);

        std::vector<Float_t> PFO_E_sznu_l0mumu(PFOsizesznu_l0mumu);
        std::vector<Float_t> PFO_C_sznu_l0mumu(PFOsizesznu_l0mumu);
        std::vector<Float_t> PFO_x_sznu_l0mumu(PFOsizesznu_l0mumu);
        std::vector<Float_t> PFO_y_sznu_l0mumu(PFOsizesznu_l0mumu);
        std::vector<Float_t> PFO_z_sznu_l0mumu(PFOsizesznu_l0mumu);

        std::vector<double> Truth_m_sznu_l0mumu(MCsizesznu_l0mumu);
        std::vector<int> Truth_PDG_sznu_l0mumu(MCsizesznu_l0mumu);
        std::vector<int> Truth_status_sznu_l0mumu(MCsizesznu_l0mumu);
        std::vector<Float_t> Truth_x_sznu_l0mumu(MCsizesznu_l0mumu);
        std::vector<Float_t> Truth_y_sznu_l0mumu(MCsizesznu_l0mumu);
        std::vector<Float_t> Truth_z_sznu_l0mumu(MCsizesznu_l0mumu);

        for (int ksznu_l0mumu=0;ksznu_l0mumu<PFOsizesznu_l0mumu;ksznu_l0mumu++){
            PFO_E_sznu_l0mumu[ksznu_l0mumu]=PFO_Esznu_l0mumu[ksznu_l0mumu];
            PFO_C_sznu_l0mumu[ksznu_l0mumu]=PFO_Csznu_l0mumu[ksznu_l0mumu];
            PFO_x_sznu_l0mumu[ksznu_l0mumu]=PFO_xsznu_l0mumu[ksznu_l0mumu];
            PFO_y_sznu_l0mumu[ksznu_l0mumu]=PFO_ysznu_l0mumu[ksznu_l0mumu];
            PFO_z_sznu_l0mumu[ksznu_l0mumu]=PFO_zsznu_l0mumu[ksznu_l0mumu];
        }

        for (int jsznu_l0mumu=0;jsznu_l0mumu<MCsizesznu_l0mumu;jsznu_l0mumu++){
            Truth_m_sznu_l0mumu[jsznu_l0mumu]=Truth_msznu_l0mumu[jsznu_l0mumu];
            Truth_PDG_sznu_l0mumu[jsznu_l0mumu]=Truth_PDGsznu_l0mumu[jsznu_l0mumu];
            Truth_status_sznu_l0mumu[jsznu_l0mumu]=Truth_statussznu_l0mumu[jsznu_l0mumu];
            Truth_x_sznu_l0mumu[jsznu_l0mumu]=Truth_xsznu_l0mumu[jsznu_l0mumu];
            Truth_y_sznu_l0mumu[jsznu_l0mumu]=Truth_ysznu_l0mumu[jsznu_l0mumu];
            Truth_z_sznu_l0mumu[jsznu_l0mumu]=Truth_zsznu_l0mumu[jsznu_l0mumu];
        }


        EvtProcessor Processorsznu_l0mumu("mm");
        Processorsznu_l0mumu.initialization();
        Processorsznu_l0mumu.setPFOinfo(PFO_E_sznu_l0mumu,PFO_C_sznu_l0mumu, PFO_x_sznu_l0mumu, PFO_y_sznu_l0mumu,PFO_z_sznu_l0mumu,PFOsizesznu_l0mumu);
        Processorsznu_l0mumu.setMCtruthinfo(Truth_m_sznu_l0mumu,Truth_PDG_sznu_l0mumu,Truth_status_sznu_l0mumu,Truth_x_sznu_l0mumu,Truth_y_sznu_l0mumu,Truth_z_sznu_l0mumu,MCsizesznu_l0mumu);
        Processorsznu_l0mumu.GenLepSelection();
        Processorsznu_l0mumu.PFOLepSelection();
        Processorsznu_l0mumu.GoodLepSelection();
        Processorsznu_l0mumu.findZleps();
        if(Processorsznu_l0mumu.ZLep.size()<2) continue;
        Processorsznu_l0mumu.GetZHMass();
        /*MH->Fill(Processor.rec_MH);
        MZ->Fill(Processor.rec_MZ);
        if ((Processorsznu_l0mumu.rec_MZ<100)&&(Processorsznu_l0mumu.rec_MZ>80)){
            invMassZ.setVal(Processor.rec_MZ);
            dataZ.add(invMassZ);
        }*/
        if(Processorsznu_l0mumu.rec_MH>upmass) continue;
        if(Processorsznu_l0mumu.rec_MH<downmass) continue;
        MH_sznu_l0mumu.push_back(Processorsznu_l0mumu.rec_MH);
        NUMsznu_l0mumu++;
        

    }
    
    for (int p_sznu_l0mumu=0; p_sznu_l0mumu<NUMsznu_l0mumu; p_sznu_l0mumu++){
        invMass.setVal(MH_sznu_l0mumu[p_sznu_l0mumu]);
        weight.setVal(Lumi*xssznu_l0mumu/DENsznu_l0mumu);
        data->add(RooArgSet(invMass,weight));
    }
    while (std::getline(fileList4mu, line4mu)) {
        size_t start4mu = line4mu.find_first_not_of(" \t");
        size_t end4mu= line4mu.find_last_not_of(" \t");
        if (start4mu == std::string::npos || line4mu.empty()) {
            continue;
        }
        std::string fileName4mu = line4mu.substr(start4mu, end4mu - start4mu + 1);

        Flow4mu->Add(fileName4mu.c_str());
    }

    fileList4mu.close();
    TTreeReader R_Flow4mu(Flow4mu);
    TTreeReaderArray<Float_t> PFO_E4mu(R_Flow4mu,  "CyberPFO.energy");
    TTreeReaderArray<Float_t> PFO_C4mu(R_Flow4mu,  "CyberPFO.charge");
    TTreeReaderArray<Float_t> PFO_x4mu(R_Flow4mu,  "CyberPFO.momentum.x");
    TTreeReaderArray<Float_t> PFO_y4mu(R_Flow4mu,  "CyberPFO.momentum.y");
    TTreeReaderArray<Float_t> PFO_z4mu(R_Flow4mu,  "CyberPFO.momentum.z");

    TTreeReaderArray<double> Truth_m4mu(R_Flow4mu,  "MCParticle.mass");
    TTreeReaderArray<int> Truth_PDG4mu(R_Flow4mu,   "MCParticle.PDG");
    TTreeReaderArray<int> Truth_status4mu(R_Flow4mu, "MCParticle.generatorStatus");
    TTreeReaderArray<Float_t> Truth_x4mu(R_Flow4mu,  "MCParticle.momentum.x");
    TTreeReaderArray<Float_t> Truth_y4mu(R_Flow4mu,  "MCParticle.momentum.y");
    TTreeReaderArray<Float_t> Truth_z4mu(R_Flow4mu,  "MCParticle.momentum.z");
    RooRealVar invMass4mu ("invMass4mu", "invMass4mu", 100, 160);
    RooDataSet data4mu("data4mu","data4mu", invMass4mu);
    //RooRealVar invMassZ ("invMassZ", "invMassZ", 80, 100);
    //RooDataSet dataZ("dataZ","dataZ", invMassZ);
    //TH1D* MH; TH1D* MZ;
    //MH = new TH1D("H_recoil_Mass","H_recoil_Mass",80,100,160);
    //MZ = new TH1D("Z_Mass","Z_Mass",80,0,400);
    std::vector<float> MH_4mu;
    while (R_Flow4mu.Next()){
        DEN4mu++;
        int PFOsize4mu=PFO_E4mu.GetSize();
        int MCsize4mu =Truth_m4mu.GetSize();
        if (PFOsize4mu==0) continue;
        std::vector<Float_t> indices_PFO(PFOsize4mu);
        std::vector<Float_t> indices_Truth(MCsize4mu);

        std::vector<Float_t> PFO_E_4mu(PFOsize4mu);
        std::vector<Float_t> PFO_C_4mu(PFOsize4mu);
        std::vector<Float_t> PFO_x_4mu(PFOsize4mu);
        std::vector<Float_t> PFO_y_4mu(PFOsize4mu);
        std::vector<Float_t> PFO_z_4mu(PFOsize4mu);

        std::vector<double> Truth_m_4mu(MCsize4mu);
        std::vector<int> Truth_PDG_4mu(MCsize4mu);
        std::vector<int> Truth_status_4mu(MCsize4mu);
        std::vector<Float_t> Truth_x_4mu(MCsize4mu);
        std::vector<Float_t> Truth_y_4mu(MCsize4mu);
        std::vector<Float_t> Truth_z_4mu(MCsize4mu);

        for (int k4mu=0;k4mu<PFOsize4mu;k4mu++){
            PFO_E_4mu[k4mu]=PFO_E4mu[k4mu];
            PFO_C_4mu[k4mu]=PFO_C4mu[k4mu];
            PFO_x_4mu[k4mu]=PFO_x4mu[k4mu];
            PFO_y_4mu[k4mu]=PFO_y4mu[k4mu];
            PFO_z_4mu[k4mu]=PFO_z4mu[k4mu];
        }

        for (int j4mu=0;j4mu<MCsize4mu;j4mu++){
            Truth_m_4mu[j4mu]=Truth_m4mu[j4mu];
            Truth_PDG_4mu[j4mu]=Truth_PDG4mu[j4mu];
            Truth_status_4mu[j4mu]=Truth_status4mu[j4mu];
            Truth_x_4mu[j4mu]=Truth_x4mu[j4mu];
            Truth_y_4mu[j4mu]=Truth_y4mu[j4mu];
            Truth_z_4mu[j4mu]=Truth_z4mu[j4mu];
        }


        EvtProcessor Processor4mu("mm");
        Processor4mu.initialization();
        Processor4mu.setPFOinfo(PFO_E_4mu,PFO_C_4mu, PFO_x_4mu, PFO_y_4mu,PFO_z_4mu,PFOsize4mu);
        Processor4mu.setMCtruthinfo(Truth_m_4mu,Truth_PDG_4mu,Truth_status_4mu,Truth_x_4mu,Truth_y_4mu,Truth_z_4mu,MCsize4mu);
        Processor4mu.GenLepSelection();
        Processor4mu.PFOLepSelection();
        Processor4mu.GoodLepSelection();
        Processor4mu.findZleps();
        if(Processor4mu.ZLep.size()<2) continue;
        Processor4mu.GetZHMass();
        /*MH->Fill(Processor.rec_MH);
        MZ->Fill(Processor.rec_MZ);
        if ((Processor4mu.rec_MZ<100)&&(Processor4mu.rec_MZ>80)){
            invMassZ.setVal(Processor.rec_MZ);
            dataZ.add(invMassZ);
        }*/
        if(Processor4mu.rec_MH>upmass) continue;
        if(Processor4mu.rec_MH<downmass) continue;
        MH_4mu.push_back(Processor4mu.rec_MH);
        NUM4mu++;
        

    }
    
    for (int p_4mu=0; p_4mu<NUM4mu; p_4mu++){
        invMass.setVal(MH_4mu[p_4mu]);
        weight.setVal(Lumi*xs4mu/DEN4mu);
        data->add(RooArgSet(invMass,weight));
    }

    while (std::getline(fileListsl0mu_down, linesl0mu_down)) {
        size_t startsl0mu_down = linesl0mu_down.find_first_not_of(" \t");
        size_t endsl0mu_down= linesl0mu_down.find_last_not_of(" \t");
        if (startsl0mu_down == std::string::npos || linesl0mu_down.empty()) {
            continue;
        }
        std::string fileNamesl0mu_down = linesl0mu_down.substr(startsl0mu_down, endsl0mu_down - startsl0mu_down + 1);

        Flowsl0mu_down->Add(fileNamesl0mu_down.c_str());
    }

    fileListsl0mu_down.close();
    TTreeReader R_Flowsl0mu_down(Flowsl0mu_down);
    TTreeReaderArray<Float_t> PFO_Esl0mu_down(R_Flowsl0mu_down,  "CyberPFO.energy");
    TTreeReaderArray<Float_t> PFO_Csl0mu_down(R_Flowsl0mu_down,  "CyberPFO.charge");
    TTreeReaderArray<Float_t> PFO_xsl0mu_down(R_Flowsl0mu_down,  "CyberPFO.momentum.x");
    TTreeReaderArray<Float_t> PFO_ysl0mu_down(R_Flowsl0mu_down,  "CyberPFO.momentum.y");
    TTreeReaderArray<Float_t> PFO_zsl0mu_down(R_Flowsl0mu_down,  "CyberPFO.momentum.z");

    TTreeReaderArray<double> Truth_msl0mu_down(R_Flowsl0mu_down,  "MCParticle.mass");
    TTreeReaderArray<int> Truth_PDGsl0mu_down(R_Flowsl0mu_down,   "MCParticle.PDG");
    TTreeReaderArray<int> Truth_statussl0mu_down(R_Flowsl0mu_down, "MCParticle.generatorStatus");
    TTreeReaderArray<Float_t> Truth_xsl0mu_down(R_Flowsl0mu_down,  "MCParticle.momentum.x");
    TTreeReaderArray<Float_t> Truth_ysl0mu_down(R_Flowsl0mu_down,  "MCParticle.momentum.y");
    TTreeReaderArray<Float_t> Truth_zsl0mu_down(R_Flowsl0mu_down,  "MCParticle.momentum.z");
    RooRealVar invMasssl0mu_down ("invMasssl0mu_down", "invMasssl0mu_down", 100, 160);
    RooDataSet datasl0mu_down("datasl0mu_down","datasl0mu_down", invMasssl0mu_down);
    //RooRealVar invMassZ ("invMassZ", "invMassZ", 80, 100);
    //RooDataSet dataZ("dataZ","dataZ", invMassZ);
    //TH1D* MH; TH1D* MZ;
    //MH = new TH1D("H_recoil_Mass","H_recoil_Mass",80,100,160);
    //MZ = new TH1D("Z_Mass","Z_Mass",80,0,400);
    std::vector<float> MH_sl0mu_down;
    while (R_Flowsl0mu_down.Next()){
        DENsl0mu_down++;
        int PFOsizesl0mu_down=PFO_Esl0mu_down.GetSize();
        int MCsizesl0mu_down =Truth_msl0mu_down.GetSize();
        if (PFOsizesl0mu_down==0) continue;
        std::vector<Float_t> indices_PFO(PFOsizesl0mu_down);
        std::vector<Float_t> indices_Truth(MCsizesl0mu_down);

        std::vector<Float_t> PFO_E_sl0mu_down(PFOsizesl0mu_down);
        std::vector<Float_t> PFO_C_sl0mu_down(PFOsizesl0mu_down);
        std::vector<Float_t> PFO_x_sl0mu_down(PFOsizesl0mu_down);
        std::vector<Float_t> PFO_y_sl0mu_down(PFOsizesl0mu_down);
        std::vector<Float_t> PFO_z_sl0mu_down(PFOsizesl0mu_down);

        std::vector<double> Truth_m_sl0mu_down(MCsizesl0mu_down);
        std::vector<int> Truth_PDG_sl0mu_down(MCsizesl0mu_down);
        std::vector<int> Truth_status_sl0mu_down(MCsizesl0mu_down);
        std::vector<Float_t> Truth_x_sl0mu_down(MCsizesl0mu_down);
        std::vector<Float_t> Truth_y_sl0mu_down(MCsizesl0mu_down);
        std::vector<Float_t> Truth_z_sl0mu_down(MCsizesl0mu_down);

        for (int ksl0mu_down=0;ksl0mu_down<PFOsizesl0mu_down;ksl0mu_down++){
            PFO_E_sl0mu_down[ksl0mu_down]=PFO_Esl0mu_down[ksl0mu_down];
            PFO_C_sl0mu_down[ksl0mu_down]=PFO_Csl0mu_down[ksl0mu_down];
            PFO_x_sl0mu_down[ksl0mu_down]=PFO_xsl0mu_down[ksl0mu_down];
            PFO_y_sl0mu_down[ksl0mu_down]=PFO_ysl0mu_down[ksl0mu_down];
            PFO_z_sl0mu_down[ksl0mu_down]=PFO_zsl0mu_down[ksl0mu_down];
        }

        for (int jsl0mu_down=0;jsl0mu_down<MCsizesl0mu_down;jsl0mu_down++){
            Truth_m_sl0mu_down[jsl0mu_down]=Truth_msl0mu_down[jsl0mu_down];
            Truth_PDG_sl0mu_down[jsl0mu_down]=Truth_PDGsl0mu_down[jsl0mu_down];
            Truth_status_sl0mu_down[jsl0mu_down]=Truth_statussl0mu_down[jsl0mu_down];
            Truth_x_sl0mu_down[jsl0mu_down]=Truth_xsl0mu_down[jsl0mu_down];
            Truth_y_sl0mu_down[jsl0mu_down]=Truth_ysl0mu_down[jsl0mu_down];
            Truth_z_sl0mu_down[jsl0mu_down]=Truth_zsl0mu_down[jsl0mu_down];
        }


        EvtProcessor Processorsl0mu_down("mm");
        Processorsl0mu_down.initialization();
        Processorsl0mu_down.setPFOinfo(PFO_E_sl0mu_down,PFO_C_sl0mu_down, PFO_x_sl0mu_down, PFO_y_sl0mu_down,PFO_z_sl0mu_down,PFOsizesl0mu_down);
        Processorsl0mu_down.setMCtruthinfo(Truth_m_sl0mu_down,Truth_PDG_sl0mu_down,Truth_status_sl0mu_down,Truth_x_sl0mu_down,Truth_y_sl0mu_down,Truth_z_sl0mu_down,MCsizesl0mu_down);
        Processorsl0mu_down.GenLepSelection();
        Processorsl0mu_down.PFOLepSelection();
        Processorsl0mu_down.GoodLepSelection();
        Processorsl0mu_down.findZleps();
        if(Processorsl0mu_down.ZLep.size()<2) continue;
        Processorsl0mu_down.GetZHMass();
        /*MH->Fill(Processor.rec_MH);
        MZ->Fill(Processor.rec_MZ);
        if ((Processorsl0mu_down.rec_MZ<100)&&(Processorsl0mu_down.rec_MZ>80)){
            invMassZ.setVal(Processor.rec_MZ);
            dataZ.add(invMassZ);
        }*/
        if(Processorsl0mu_down.rec_MH>upmass) continue;
        if(Processorsl0mu_down.rec_MH<downmass) continue;
        MH_sl0mu_down.push_back(Processorsl0mu_down.rec_MH);
        NUMsl0mu_down++;
        

    }
    
    for (int p_sl0mu_down=0; p_sl0mu_down<NUMsl0mu_down; p_sl0mu_down++){
        invMass.setVal(MH_sl0mu_down[p_sl0mu_down]);
        weight.setVal(Lumi*xssl0mu_down/DENsl0mu_down);
        data->add(RooArgSet(invMass,weight));
    }

    while (std::getline(fileListsl0mu_up, linesl0mu_up)) {
        size_t startsl0mu_up = linesl0mu_up.find_first_not_of(" \t");
        size_t endsl0mu_up= linesl0mu_up.find_last_not_of(" \t");
        if (startsl0mu_up == std::string::npos || linesl0mu_up.empty()) {
            continue;
        }
        std::string fileNamesl0mu_up = linesl0mu_up.substr(startsl0mu_up, endsl0mu_up - startsl0mu_up + 1);

        Flowsl0mu_up->Add(fileNamesl0mu_up.c_str());
    }

    fileListsl0mu_up.close();
    TTreeReader R_Flowsl0mu_up(Flowsl0mu_up);
    TTreeReaderArray<Float_t> PFO_Esl0mu_up(R_Flowsl0mu_up,  "CyberPFO.energy");
    TTreeReaderArray<Float_t> PFO_Csl0mu_up(R_Flowsl0mu_up,  "CyberPFO.charge");
    TTreeReaderArray<Float_t> PFO_xsl0mu_up(R_Flowsl0mu_up,  "CyberPFO.momentum.x");
    TTreeReaderArray<Float_t> PFO_ysl0mu_up(R_Flowsl0mu_up,  "CyberPFO.momentum.y");
    TTreeReaderArray<Float_t> PFO_zsl0mu_up(R_Flowsl0mu_up,  "CyberPFO.momentum.z");

    TTreeReaderArray<double> Truth_msl0mu_up(R_Flowsl0mu_up,  "MCParticle.mass");
    TTreeReaderArray<int> Truth_PDGsl0mu_up(R_Flowsl0mu_up,   "MCParticle.PDG");
    TTreeReaderArray<int> Truth_statussl0mu_up(R_Flowsl0mu_up, "MCParticle.generatorStatus");
    TTreeReaderArray<Float_t> Truth_xsl0mu_up(R_Flowsl0mu_up,  "MCParticle.momentum.x");
    TTreeReaderArray<Float_t> Truth_ysl0mu_up(R_Flowsl0mu_up,  "MCParticle.momentum.y");
    TTreeReaderArray<Float_t> Truth_zsl0mu_up(R_Flowsl0mu_up,  "MCParticle.momentum.z");
    RooRealVar invMasssl0mu_up ("invMasssl0mu_up", "invMasssl0mu_up", 100, 160);
    RooDataSet datasl0mu_up("datasl0mu_up","datasl0mu_up", invMasssl0mu_up);
    //RooRealVar invMassZ ("invMassZ", "invMassZ", 80, 100);
    //RooDataSet dataZ("dataZ","dataZ", invMassZ);
    //TH1D* MH; TH1D* MZ;
    //MH = new TH1D("H_recoil_Mass","H_recoil_Mass",80,100,160);
    //MZ = new TH1D("Z_Mass","Z_Mass",80,0,400);
    std::vector<float> MH_sl0mu_up;
    while (R_Flowsl0mu_up.Next()){
        DENsl0mu_up++;
        int PFOsizesl0mu_up=PFO_Esl0mu_up.GetSize();
        int MCsizesl0mu_up =Truth_msl0mu_up.GetSize();
        if (PFOsizesl0mu_up==0) continue;
        std::vector<Float_t> indices_PFO(PFOsizesl0mu_up);
        std::vector<Float_t> indices_Truth(MCsizesl0mu_up);

        std::vector<Float_t> PFO_E_sl0mu_up(PFOsizesl0mu_up);
        std::vector<Float_t> PFO_C_sl0mu_up(PFOsizesl0mu_up);
        std::vector<Float_t> PFO_x_sl0mu_up(PFOsizesl0mu_up);
        std::vector<Float_t> PFO_y_sl0mu_up(PFOsizesl0mu_up);
        std::vector<Float_t> PFO_z_sl0mu_up(PFOsizesl0mu_up);

        std::vector<double> Truth_m_sl0mu_up(MCsizesl0mu_up);
        std::vector<int> Truth_PDG_sl0mu_up(MCsizesl0mu_up);
        std::vector<int> Truth_status_sl0mu_up(MCsizesl0mu_up);
        std::vector<Float_t> Truth_x_sl0mu_up(MCsizesl0mu_up);
        std::vector<Float_t> Truth_y_sl0mu_up(MCsizesl0mu_up);
        std::vector<Float_t> Truth_z_sl0mu_up(MCsizesl0mu_up);

        for (int ksl0mu_up=0;ksl0mu_up<PFOsizesl0mu_up;ksl0mu_up++){
            PFO_E_sl0mu_up[ksl0mu_up]=PFO_Esl0mu_up[ksl0mu_up];
            PFO_C_sl0mu_up[ksl0mu_up]=PFO_Csl0mu_up[ksl0mu_up];
            PFO_x_sl0mu_up[ksl0mu_up]=PFO_xsl0mu_up[ksl0mu_up];
            PFO_y_sl0mu_up[ksl0mu_up]=PFO_ysl0mu_up[ksl0mu_up];
            PFO_z_sl0mu_up[ksl0mu_up]=PFO_zsl0mu_up[ksl0mu_up];
        }

        for (int jsl0mu_up=0;jsl0mu_up<MCsizesl0mu_up;jsl0mu_up++){
            Truth_m_sl0mu_up[jsl0mu_up]=Truth_msl0mu_up[jsl0mu_up];
            Truth_PDG_sl0mu_up[jsl0mu_up]=Truth_PDGsl0mu_up[jsl0mu_up];
            Truth_status_sl0mu_up[jsl0mu_up]=Truth_statussl0mu_up[jsl0mu_up];
            Truth_x_sl0mu_up[jsl0mu_up]=Truth_xsl0mu_up[jsl0mu_up];
            Truth_y_sl0mu_up[jsl0mu_up]=Truth_ysl0mu_up[jsl0mu_up];
            Truth_z_sl0mu_up[jsl0mu_up]=Truth_zsl0mu_up[jsl0mu_up];
        }


        EvtProcessor Processorsl0mu_up("mm");
        Processorsl0mu_up.initialization();
        Processorsl0mu_up.setPFOinfo(PFO_E_sl0mu_up,PFO_C_sl0mu_up, PFO_x_sl0mu_up, PFO_y_sl0mu_up,PFO_z_sl0mu_up,PFOsizesl0mu_up);
        Processorsl0mu_up.setMCtruthinfo(Truth_m_sl0mu_up,Truth_PDG_sl0mu_up,Truth_status_sl0mu_up,Truth_x_sl0mu_up,Truth_y_sl0mu_up,Truth_z_sl0mu_up,MCsizesl0mu_up);
        Processorsl0mu_up.GenLepSelection();
        Processorsl0mu_up.PFOLepSelection();
        Processorsl0mu_up.GoodLepSelection();
        Processorsl0mu_up.findZleps();
        if(Processorsl0mu_up.ZLep.size()<2) continue;
        Processorsl0mu_up.GetZHMass();
        /*MH->Fill(Processor.rec_MH);
        MZ->Fill(Processor.rec_MZ);
        if ((Processorsl0mu_up.rec_MZ<100)&&(Processorsl0mu_up.rec_MZ>80)){
            invMassZ.setVal(Processor.rec_MZ);
            dataZ.add(invMassZ);
        }*/
        if(Processorsl0mu_up.rec_MH>upmass) continue;
        if(Processorsl0mu_up.rec_MH<downmass) continue;
        MH_sl0mu_up.push_back(Processorsl0mu_up.rec_MH);
        NUMsl0mu_up++;
        

    }
    
    for (int p_sl0mu_up=0; p_sl0mu_up<NUMsl0mu_up; p_sl0mu_up++){
        invMass.setVal(MH_sl0mu_up[p_sl0mu_up]);
        weight.setVal(Lumi*xssl0mu_up/DENsl0mu_up);
        data->add(RooArgSet(invMass,weight));
    }

    while (std::getline(fileListzzorww_l0, linezzorww_l0)) {
        size_t startzzorww_l0 = linezzorww_l0.find_first_not_of(" \t");
        size_t endzzorww_l0= linezzorww_l0.find_last_not_of(" \t");
        if (startzzorww_l0 == std::string::npos || linezzorww_l0.empty()) {
            continue;
        }
        std::string fileNamezzorww_l0 = linezzorww_l0.substr(startzzorww_l0, endzzorww_l0 - startzzorww_l0 + 1);

        Flowzzorww_l0->Add(fileNamezzorww_l0.c_str());
    }

    fileListzzorww_l0.close();
    TTreeReader R_Flowzzorww_l0(Flowzzorww_l0);
    TTreeReaderArray<Float_t> PFO_Ezzorww_l0(R_Flowzzorww_l0,  "CyberPFO.energy");
    TTreeReaderArray<Float_t> PFO_Czzorww_l0(R_Flowzzorww_l0,  "CyberPFO.charge");
    TTreeReaderArray<Float_t> PFO_xzzorww_l0(R_Flowzzorww_l0,  "CyberPFO.momentum.x");
    TTreeReaderArray<Float_t> PFO_yzzorww_l0(R_Flowzzorww_l0,  "CyberPFO.momentum.y");
    TTreeReaderArray<Float_t> PFO_zzzorww_l0(R_Flowzzorww_l0,  "CyberPFO.momentum.z");

    TTreeReaderArray<double> Truth_mzzorww_l0(R_Flowzzorww_l0,  "MCParticle.mass");
    TTreeReaderArray<int> Truth_PDGzzorww_l0(R_Flowzzorww_l0,   "MCParticle.PDG");
    TTreeReaderArray<int> Truth_statuszzorww_l0(R_Flowzzorww_l0, "MCParticle.generatorStatus");
    TTreeReaderArray<Float_t> Truth_xzzorww_l0(R_Flowzzorww_l0,  "MCParticle.momentum.x");
    TTreeReaderArray<Float_t> Truth_yzzorww_l0(R_Flowzzorww_l0,  "MCParticle.momentum.y");
    TTreeReaderArray<Float_t> Truth_zzzorww_l0(R_Flowzzorww_l0,  "MCParticle.momentum.z");
    RooRealVar invMasszzorww_l0 ("invMasszzorww_l0", "invMasszzorww_l0", 100, 160);
    RooDataSet datazzorww_l0("datazzorww_l0","datazzorww_l0", invMasszzorww_l0);
    //RooRealVar invMassZ ("invMassZ", "invMassZ", 80, 100);
    //RooDataSet dataZ("dataZ","dataZ", invMassZ);
    //TH1D* MH; TH1D* MZ;
    //MH = new TH1D("H_recoil_Mass","H_recoil_Mass",80,100,160);
    //MZ = new TH1D("Z_Mass","Z_Mass",80,0,400);
    std::vector<float> MH_zzorww_l0;
    while (R_Flowzzorww_l0.Next()){
        DENzzorww_l0++;
        int PFOsizezzorww_l0=PFO_Ezzorww_l0.GetSize();
        int MCsizezzorww_l0 =Truth_mzzorww_l0.GetSize();
        if (PFOsizezzorww_l0==0) continue;
        std::vector<Float_t> indices_PFO(PFOsizezzorww_l0);
        std::vector<Float_t> indices_Truth(MCsizezzorww_l0);

        std::vector<Float_t> PFO_E_zzorww_l0(PFOsizezzorww_l0);
        std::vector<Float_t> PFO_C_zzorww_l0(PFOsizezzorww_l0);
        std::vector<Float_t> PFO_x_zzorww_l0(PFOsizezzorww_l0);
        std::vector<Float_t> PFO_y_zzorww_l0(PFOsizezzorww_l0);
        std::vector<Float_t> PFO_z_zzorww_l0(PFOsizezzorww_l0);

        std::vector<double> Truth_m_zzorww_l0(MCsizezzorww_l0);
        std::vector<int> Truth_PDG_zzorww_l0(MCsizezzorww_l0);
        std::vector<int> Truth_status_zzorww_l0(MCsizezzorww_l0);
        std::vector<Float_t> Truth_x_zzorww_l0(MCsizezzorww_l0);
        std::vector<Float_t> Truth_y_zzorww_l0(MCsizezzorww_l0);
        std::vector<Float_t> Truth_z_zzorww_l0(MCsizezzorww_l0);

        for (int kzzorww_l0=0;kzzorww_l0<PFOsizezzorww_l0;kzzorww_l0++){
            PFO_E_zzorww_l0[kzzorww_l0]=PFO_Ezzorww_l0[kzzorww_l0];
            PFO_C_zzorww_l0[kzzorww_l0]=PFO_Czzorww_l0[kzzorww_l0];
            PFO_x_zzorww_l0[kzzorww_l0]=PFO_xzzorww_l0[kzzorww_l0];
            PFO_y_zzorww_l0[kzzorww_l0]=PFO_yzzorww_l0[kzzorww_l0];
            PFO_z_zzorww_l0[kzzorww_l0]=PFO_zzzorww_l0[kzzorww_l0];
        }

        for (int jzzorww_l0=0;jzzorww_l0<MCsizezzorww_l0;jzzorww_l0++){
            Truth_m_zzorww_l0[jzzorww_l0]=Truth_mzzorww_l0[jzzorww_l0];
            Truth_PDG_zzorww_l0[jzzorww_l0]=Truth_PDGzzorww_l0[jzzorww_l0];
            Truth_status_zzorww_l0[jzzorww_l0]=Truth_statuszzorww_l0[jzzorww_l0];
            Truth_x_zzorww_l0[jzzorww_l0]=Truth_xzzorww_l0[jzzorww_l0];
            Truth_y_zzorww_l0[jzzorww_l0]=Truth_yzzorww_l0[jzzorww_l0];
            Truth_z_zzorww_l0[jzzorww_l0]=Truth_zzzorww_l0[jzzorww_l0];
        }


        EvtProcessor Processorzzorww_l0("mm");
        Processorzzorww_l0.initialization();
        Processorzzorww_l0.setPFOinfo(PFO_E_zzorww_l0,PFO_C_zzorww_l0, PFO_x_zzorww_l0, PFO_y_zzorww_l0,PFO_z_zzorww_l0,PFOsizezzorww_l0);
        Processorzzorww_l0.setMCtruthinfo(Truth_m_zzorww_l0,Truth_PDG_zzorww_l0,Truth_status_zzorww_l0,Truth_x_zzorww_l0,Truth_y_zzorww_l0,Truth_z_zzorww_l0,MCsizezzorww_l0);
        Processorzzorww_l0.GenLepSelection();
        Processorzzorww_l0.PFOLepSelection();
        Processorzzorww_l0.GoodLepSelection();
        Processorzzorww_l0.findZleps();
        if(Processorzzorww_l0.ZLep.size()<2) continue;
        Processorzzorww_l0.GetZHMass();
        /*MH->Fill(Processor.rec_MH);
        MZ->Fill(Processor.rec_MZ);
        if ((Processorzzorww_l0.rec_MZ<100)&&(Processorzzorww_l0.rec_MZ>80)){
            invMassZ.setVal(Processor.rec_MZ);
            dataZ.add(invMassZ);
        }*/
        if(Processorzzorww_l0.rec_MH>upmass) continue;
        if(Processorzzorww_l0.rec_MH<downmass) continue;
        MH_zzorww_l0.push_back(Processorzzorww_l0.rec_MH);
        NUMzzorww_l0++;
        

    }
    
    for (int p_zzorww_l0=0; p_zzorww_l0<NUMzzorww_l0; p_zzorww_l0++){
        invMass.setVal(MH_zzorww_l0[p_zzorww_l0]);
        weight.setVal(Lumi*xszzorww_l0/DENzzorww_l0);
        data->add(RooArgSet(invMass,weight));
    }

    while (std::getline(fileListll, linell)) {
        size_t startll = linell.find_first_not_of(" \t");
        size_t endll= linell.find_last_not_of(" \t");
        if (startll == std::string::npos || linell.empty()) {
            continue;
        }
        std::string fileNamell = linell.substr(startll, endll - startll + 1);

        Flowll->Add(fileNamell.c_str());
    }

    fileListll.close();
    TTreeReader R_Flowll(Flowll);
    TTreeReaderArray<Float_t> PFO_Ell(R_Flowll,  "CyberPFO.energy");
    TTreeReaderArray<Float_t> PFO_Cll(R_Flowll,  "CyberPFO.charge");
    TTreeReaderArray<Float_t> PFO_xll(R_Flowll,  "CyberPFO.momentum.x");
    TTreeReaderArray<Float_t> PFO_yll(R_Flowll,  "CyberPFO.momentum.y");
    TTreeReaderArray<Float_t> PFO_zll(R_Flowll,  "CyberPFO.momentum.z");

    TTreeReaderArray<double> Truth_mll(R_Flowll,  "MCParticle.mass");
    TTreeReaderArray<int> Truth_PDGll(R_Flowll,   "MCParticle.PDG");
    TTreeReaderArray<int> Truth_statusll(R_Flowll, "MCParticle.generatorStatus");
    TTreeReaderArray<Float_t> Truth_xll(R_Flowll,  "MCParticle.momentum.x");
    TTreeReaderArray<Float_t> Truth_yll(R_Flowll,  "MCParticle.momentum.y");
    TTreeReaderArray<Float_t> Truth_zll(R_Flowll,  "MCParticle.momentum.z");
    RooRealVar invMassll ("invMassll", "invMassll", 100, 160);
    RooDataSet datall("datall","datall", invMassll);
    //RooRealVar invMassZ ("invMassZ", "invMassZ", 80, 100);
    //RooDataSet dataZ("dataZ","dataZ", invMassZ);
    //TH1D* MH; TH1D* MZ;
    //MH = new TH1D("H_recoil_Mass","H_recoil_Mass",80,100,160);
    //MZ = new TH1D("Z_Mass","Z_Mass",80,0,400);
    std::vector<float> MH_ll;
    while (R_Flowll.Next()){
        DENll++;
        int PFOsizell=PFO_Ell.GetSize();
        int MCsizell =Truth_mll.GetSize();
        if (PFOsizell==0) continue;
        std::vector<Float_t> indices_PFO(PFOsizell);
        std::vector<Float_t> indices_Truth(MCsizell);

        std::vector<Float_t> PFO_E_ll(PFOsizell);
        std::vector<Float_t> PFO_C_ll(PFOsizell);
        std::vector<Float_t> PFO_x_ll(PFOsizell);
        std::vector<Float_t> PFO_y_ll(PFOsizell);
        std::vector<Float_t> PFO_z_ll(PFOsizell);

        std::vector<double> Truth_m_ll(MCsizell);
        std::vector<int> Truth_PDG_ll(MCsizell);
        std::vector<int> Truth_status_ll(MCsizell);
        std::vector<Float_t> Truth_x_ll(MCsizell);
        std::vector<Float_t> Truth_y_ll(MCsizell);
        std::vector<Float_t> Truth_z_ll(MCsizell);

        for (int kll=0;kll<PFOsizell;kll++){
            PFO_E_ll[kll]=PFO_Ell[kll];
            PFO_C_ll[kll]=PFO_Cll[kll];
            PFO_x_ll[kll]=PFO_xll[kll];
            PFO_y_ll[kll]=PFO_yll[kll];
            PFO_z_ll[kll]=PFO_zll[kll];
        }

        for (int jll=0;jll<MCsizell;jll++){
            Truth_m_ll[jll]=Truth_mll[jll];
            Truth_PDG_ll[jll]=Truth_PDGll[jll];
            Truth_status_ll[jll]=Truth_statusll[jll];
            Truth_x_ll[jll]=Truth_xll[jll];
            Truth_y_ll[jll]=Truth_yll[jll];
            Truth_z_ll[jll]=Truth_zll[jll];
        }


        EvtProcessor Processorll("mm");
        Processorll.initialization();
        Processorll.setPFOinfo(PFO_E_ll,PFO_C_ll, PFO_x_ll, PFO_y_ll,PFO_z_ll,PFOsizell);
        Processorll.setMCtruthinfo(Truth_m_ll,Truth_PDG_ll,Truth_status_ll,Truth_x_ll,Truth_y_ll,Truth_z_ll,MCsizell);
        Processorll.GenLepSelection();
        Processorll.PFOLepSelection();
        Processorll.GoodLepSelection();
        Processorll.findZleps();
        if(Processorll.ZLep.size()<2) continue;
        Processorll.GetZHMass();
        /*MH->Fill(Processor.rec_MH);
        MZ->Fill(Processor.rec_MZ);
        if ((Processorll.rec_MZ<100)&&(Processorll.rec_MZ>80)){
            invMassZ.setVal(Processor.rec_MZ);
            dataZ.add(invMassZ);
        }*/
        if(Processorll.rec_MH>upmass) continue;
        if(Processorll.rec_MH<downmass) continue;
        MH_ll.push_back(Processorll.rec_MH);
        NUMll++;
        

    }
    
    for (int p_ll=0; p_ll<NUMll; p_ll++){
        invMass.setVal(MH_ll[p_ll]);
        weight.setVal(Lumi*xsll/DENll);
        //data->add(RooArgSet(invMass,weight));
    }

    while (std::getline(fileListsig, linesig)) {
        size_t startsig = linesig.find_first_not_of(" \t");
        size_t endsig= linesig.find_last_not_of(" \t");
        if (startsig == std::string::npos || linesig.empty()) {
            continue;
        }
        std::string fileNamesig = linesig.substr(startsig, endsig - startsig + 1);

        Flowsig->Add(fileNamesig.c_str());
    }

    fileListsig.close();
    TTreeReader R_Flowsig(Flowsig);
    TTreeReaderArray<Float_t> PFO_Esig(R_Flowsig,  "CyberPFO.energy");
    TTreeReaderArray<Float_t> PFO_Csig(R_Flowsig,  "CyberPFO.charge");
    TTreeReaderArray<Float_t> PFO_xsig(R_Flowsig,  "CyberPFO.momentum.x");
    TTreeReaderArray<Float_t> PFO_ysig(R_Flowsig,  "CyberPFO.momentum.y");
    TTreeReaderArray<Float_t> PFO_zsig(R_Flowsig,  "CyberPFO.momentum.z");

    TTreeReaderArray<double> Truth_msig(R_Flowsig,  "MCParticle.mass");
    TTreeReaderArray<int> Truth_PDGsig(R_Flowsig,   "MCParticle.PDG");
    TTreeReaderArray<int> Truth_statussig(R_Flowsig, "MCParticle.generatorStatus");
    TTreeReaderArray<Float_t> Truth_xsig(R_Flowsig,  "MCParticle.momentum.x");
    TTreeReaderArray<Float_t> Truth_ysig(R_Flowsig,  "MCParticle.momentum.y");
    TTreeReaderArray<Float_t> Truth_zsig(R_Flowsig,  "MCParticle.momentum.z");
    
    //RooRealVar invMassZ ("invMassZ", "invMassZ", 80, 100);
    //RooDataSet dataZ("dataZ","dataZ", invMassZ);
    //TH1D* MH; TH1D* MZ;
    //MH = new TH1D("H_recoil_Mass","H_recoil_Mass",80,100,160);
    //MZ = new TH1D("Z_Mass","Z_Mass",80,0,400);
    std::vector<float> MH_sig;
    while (R_Flowsig.Next()){
        DENsig++;
        int PFOsizesig=PFO_Esig.GetSize();
        int MCsizesig =Truth_msig.GetSize();
        if (PFOsizesig==0) continue;
        std::vector<Float_t> indices_PFO(PFOsizesig);
        std::vector<Float_t> indices_Truth(MCsizesig);

        std::vector<Float_t> PFO_E_sig(PFOsizesig);
        std::vector<Float_t> PFO_C_sig(PFOsizesig);
        std::vector<Float_t> PFO_x_sig(PFOsizesig);
        std::vector<Float_t> PFO_y_sig(PFOsizesig);
        std::vector<Float_t> PFO_z_sig(PFOsizesig);

        std::vector<double> Truth_m_sig(MCsizesig);
        std::vector<int> Truth_PDG_sig(MCsizesig);
        std::vector<int> Truth_status_sig(MCsizesig);
        std::vector<Float_t> Truth_x_sig(MCsizesig);
        std::vector<Float_t> Truth_y_sig(MCsizesig);
        std::vector<Float_t> Truth_z_sig(MCsizesig);

        for (int ksig=0;ksig<PFOsizesig;ksig++){
            PFO_E_sig[ksig]=PFO_Esig[ksig];
            PFO_C_sig[ksig]=PFO_Csig[ksig];
            PFO_x_sig[ksig]=PFO_xsig[ksig];
            PFO_y_sig[ksig]=PFO_ysig[ksig];
            PFO_z_sig[ksig]=PFO_zsig[ksig];
        }

        for (int jsig=0;jsig<MCsizesig;jsig++){
            Truth_m_sig[jsig]=Truth_msig[jsig];
            Truth_PDG_sig[jsig]=Truth_PDGsig[jsig];
            Truth_status_sig[jsig]=Truth_statussig[jsig];
            Truth_x_sig[jsig]=Truth_xsig[jsig];
            Truth_y_sig[jsig]=Truth_ysig[jsig];
            Truth_z_sig[jsig]=Truth_zsig[jsig];
        }


        EvtProcessor Processorsig("mm");
        Processorsig.initialization();
        Processorsig.setPFOinfo(PFO_E_sig,PFO_C_sig, PFO_x_sig, PFO_y_sig,PFO_z_sig,PFOsizesig);
        Processorsig.setMCtruthinfo(Truth_m_sig,Truth_PDG_sig,Truth_status_sig,Truth_x_sig,Truth_y_sig,Truth_z_sig,MCsizesig);
        Processorsig.GenLepSelection();
        Processorsig.PFOLepSelection();
        Processorsig.GoodLepSelection();
        Processorsig.findZleps();
        if(Processorsig.ZLep.size()<2) continue;
        Processorsig.GetZHMass();
        /*MH->Fill(Processor.rec_MH);
        MZ->Fill(Processor.rec_MZ);
        if ((Processorsig.rec_MZ<100)&&(Processorsig.rec_MZ>80)){
            invMassZ.setVal(Processor.rec_MZ);
            dataZ.add(invMassZ);
        }*/
        if(Processorsig.rec_MH>upmass) continue;
        if(Processorsig.rec_MH<downmass) continue;
        MH_sig.push_back(Processorsig.rec_MH);
        NUMsig++;
        

    }
    
    for (int p_sig=0; p_sig<NUMsig; p_sig++){
        invMass.setVal(MH_sig[p_sig]);
        weight.setVal(Lumi*xssig/DENsig);
        data->add(RooArgSet(invMass,weight));
    }
    //MH->SetDirectory(0);
    //MZ->SetDirectory(0);
    TCanvas* c = new TCanvas("c_roofit", "RooFit Gaussian Fit", 800, 600);
    //c->SetLogy();
    RooRealVar lambda("lambda","",-1.93879e-02,-1.93880e-02, -1.93878e-02);
    RooAbsPdf *bkgpdf = new RooExponential("exp_pdf", "Exponential PDF", invMass, lambda);
    // RooAbsPdf *pdf =new RooGaussian("pdf", "pdf", invMass, meanCB, sigmaCB);
    RooRealVar meanCB("meanCB","",125.,120.,130.);
    RooRealVar sigmaCB("sigmaCB","",0.5,0.,5.);
    RooRealVar alphaCBLo("alphaCBLo","",1.,0.,5.);
    RooRealVar alphaCBHi("alphaCBHi","",1.,0.,10.);
    RooRealVar nCBLo("nCBLo","",2.,0.,5.);
    RooRealVar nCBHi("nCBHi","",2.,0.,10.);
    RooAbsPdf *sigpdf =new HggTwoSidedCBPdf("sigpdf", "sigpdf", invMass, meanCB, sigmaCB, alphaCBLo, nCBLo, alphaCBHi, nCBHi);

    RooRealVar nsig("nsig", "Number of signal events", 500, 0, 300000);
    RooRealVar nbkg("nbkg", "Number of background events", 500, 0, 300000);

    RooAddPdf* model = new RooAddPdf("model", "Signal+Background", RooArgList(*sigpdf , *bkgpdf), RooArgList(nsig, nbkg));
    RooDataSet data_w("data_w","data_w",data,*data->get(),"",weight.GetName());
    RooFitResult *result=model->fitTo(data_w,RooFit::SumW2Error(kTRUE));

    RooPlot *frame=invMass.frame(50); 
    
    
    data_w.plotOn(frame, RooFit::DataError(RooAbsData::SumW2));
    model->plotOn(frame);
    model->plotOn(frame,
        RooFit::Components(*bkgpdf),
        RooFit::LineColor(kBlue),
        RooFit::LineStyle(kDashed));
    model->plotOn(frame,
        RooFit::Components(*sigpdf),
        RooFit::LineColor(kRed),
        RooFit::LineStyle(kDashed));
    frame->GetYaxis()->SetTitle("Events");
    frame->GetXaxis()->SetTitle("Recoil Mass [GeV]");
    frame->Draw();
    TLegend *leg = new TLegend(0.15, 0.7, 0.4, 0.85); // (x1,y1,x2,y2) 左上方坐标范围
    leg->SetBorderSize(0);                            // 去除边框
    leg->SetFillStyle(0);                             // 透明填充

    // 添加条目（参数：图形对象，标签，绘制样式）
    leg->AddEntry(frame->findObject("data"),   "Data",          "P");   // "P"表示点
    leg->AddEntry(frame->findObject("model"), "Total Fit",     "L");   // "L"表示线
    leg->AddEntry(frame->findObject("bkgpdf"), "Background",   "L");   // 背景虚线
    leg->AddEntry(frame->findObject("sigpdf"), "Signal",       "L");   // 信号虚线


    leg->Draw();
    TLatex *tex = new TLatex();
    FormatLatex(tex);
    tex->DrawLatexNDC(0.6, 0.8, "#it{#bf{CEPC}} Ref-TDR");
    tex->DrawLatexNDC(0.6, 0.74, "H recoil ");
    // tex->DrawLatexNDC(0.75, 0.66, (TString)"Output: "+ (TString)to_string(evt_h1)+(TString)" events");
    TPad *pad0 = (TPad*)c->GetPad(0);
    pad0->SetLeftMargin(0.15); 
    pad0->SetRightMargin(0.15);  
    pad0->SetBottomMargin(0.15); 
    pad0->SetTopMargin(0.1);  
    c->SaveAs("H_recoilMass_bkg.png");
    delete c;
    
    cout<<Lumi*xszz_l0mumu*(NUMzz_l0mumu/DENzz_l0mumu)<<endl;
    cout<<Lumi*xssig*(NUMsig/DENsig)<<endl;
}

