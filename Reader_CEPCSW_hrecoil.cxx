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

    RooRealVar x("x", "invMass", 120, 130);


    
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
void Reader_CEPCSW_hrecoil()
{   
    gROOT->ProcessLine(".L EvtProcessor.cxx");
    SetCEPCStyle();
    processing("Hrecoil");

}
void processing(TString outputTag){
    TChain *Flow = new TChain("events", "events");
    std::ifstream fileList("sample_input_mmHinv.txt");
    std::string line;
    while (std::getline(fileList, line)) {
        size_t start = line.find_first_not_of(" \t");
        size_t end = line.find_last_not_of(" \t");
        if (start == std::string::npos || line.empty()) {
            continue;
        }
        std::string fileName = line.substr(start, end - start + 1);

        Flow->Add(fileName.c_str());
    }

    fileList.close();
    TTreeReader R_Flow(Flow);
    TTreeReaderArray<Float_t> PFO_E(R_Flow,  "CyberPFO.energy");
    TTreeReaderArray<Float_t> PFO_C(R_Flow,  "CyberPFO.charge");
    TTreeReaderArray<Float_t> PFO_x(R_Flow,  "CyberPFO.momentum.x");
    TTreeReaderArray<Float_t> PFO_y(R_Flow,  "CyberPFO.momentum.y");
    TTreeReaderArray<Float_t> PFO_z(R_Flow,  "CyberPFO.momentum.z");

    TTreeReaderArray<double> Truth_m(R_Flow,  "MCParticle.mass");
    TTreeReaderArray<int> Truth_PDG(R_Flow,   "MCParticle.PDG");
    TTreeReaderArray<int> Truth_status(R_Flow, "MCParticle.generatorStatus");
    TTreeReaderArray<Float_t> Truth_x(R_Flow,  "MCParticle.momentum.x");
    TTreeReaderArray<Float_t> Truth_y(R_Flow,  "MCParticle.momentum.y");
    TTreeReaderArray<Float_t> Truth_z(R_Flow,  "MCParticle.momentum.z");
    RooRealVar invMass ("invMass", "invMass", 120, 130);
    RooDataSet data("data","data", invMass);
    RooRealVar invMassZ ("invMassZ", "invMassZ", 80, 100);
    RooDataSet dataZ("dataZ","dataZ", invMassZ);
    int TotalEvents=0;
    TH1D* MH; TH1D* MZ;
    MH = new TH1D("H_recoil_Mass","H_recoil_Mass",80,100,160);
    MZ = new TH1D("Z_Mass","Z_Mass",80,70,110);
    while (R_Flow.Next()){
        if(TotalEvents%100==0) std::cout<<"Processed "<<TotalEvents<<" Events!"<<std::endl;
        TotalEvents++;
        int PFOsize=PFO_E.GetSize();
        int MCsize =Truth_m.GetSize();
        if (PFOsize==0) continue;
        std::vector<Float_t> indices_PFO(PFOsize);
        std::vector<Float_t> indices_Truth(MCsize);

        std::vector<Float_t> PFO_E_(PFOsize);
        std::vector<Float_t> PFO_C_(PFOsize);
        std::vector<Float_t> PFO_x_(PFOsize);
        std::vector<Float_t> PFO_y_(PFOsize);
        std::vector<Float_t> PFO_z_(PFOsize);

        std::vector<double> Truth_m_(MCsize);
        std::vector<int> Truth_PDG_(MCsize);
        std::vector<int> Truth_status_(MCsize);
        std::vector<Float_t> Truth_x_(MCsize);
        std::vector<Float_t> Truth_y_(MCsize);
        std::vector<Float_t> Truth_z_(MCsize);

        for (int k=0;k<PFOsize;k++){
            PFO_E_[k]=PFO_E[k];
            PFO_C_[k]=PFO_C[k];
            PFO_x_[k]=PFO_x[k];
            PFO_y_[k]=PFO_y[k];
            PFO_z_[k]=PFO_z[k];
        }

        for (int j=0;j<MCsize;j++){
            Truth_m_[j]=Truth_m[j];
            Truth_PDG_[j]=Truth_PDG[j];
            Truth_status_[j]=Truth_status[j];
            Truth_x_[j]=Truth_x[j];
            Truth_y_[j]=Truth_y[j];
            Truth_z_[j]=Truth_z[j];
        }


        EvtProcessor Processor("mm");
        Processor.initialization();
        Processor.setPFOinfo(PFO_E_,PFO_C_, PFO_x_, PFO_y_,PFO_z_,PFOsize);
        Processor.setMCtruthinfo(Truth_m_,Truth_PDG_,Truth_status_,Truth_x_,Truth_y_,Truth_z_,MCsize);
        Processor.GenLepSelection();
        Processor.PFOLepSelection();
        Processor.GoodLepSelection();
        Processor.findZleps();
        if(Processor.ZLep.size()<2) continue;
        Processor.GetZHMass();
        MH->Fill(Processor.rec_MH);
        MZ->Fill(Processor.rec_MZ);
        if ((Processor.rec_MZ<100)&&(Processor.rec_MZ>80)){
            invMassZ.setVal(Processor.rec_MZ);
            dataZ.add(invMassZ);
        }
        if(Processor.rec_MH>130) continue;
        if(Processor.rec_MH<120) continue;
        invMass.setVal(Processor.rec_MH);
        data.add(invMass);

    }
    MH->SetDirectory(0);
    MZ->SetDirectory(0);
    TCanvas* c = new TCanvas("c_roofit", "RooFit Gaussian Fit", 800, 600);

    
    RooRealVar meanCB("meanCB","",125.,120.,130.);
    RooRealVar sigmaCB("sigmaCB","",0.5,0.,5.);
    RooRealVar alphaCBLo("alphaCBLo","",1.,0.,5.);
    RooRealVar alphaCBHi("alphaCBHi","",1.,0.,10.);
    RooRealVar nCBLo("nCBLo","",2.,0.,5.);
    RooRealVar nCBHi("nCBHi","",2.,0.,10.);

    RooAbsPdf *pdf =new HggTwoSidedCBPdf("pdf", "pdf", invMass, meanCB, sigmaCB, alphaCBLo, nCBLo, alphaCBHi, nCBHi);
    // RooAbsPdf *pdf =new RooGaussian("pdf", "pdf", invMass, meanCB, sigmaCB);

    RooFitResult *result=pdf->fitTo(data);

    RooPlot *frame=invMass.frame(50);
    
    
    data.plotOn(frame);
    pdf->plotOn(frame);
    frame->GetYaxis()->SetTitle("Events");
    frame->GetXaxis()->SetTitle("Recoil Mass [GeV]");
    frame->Draw();
    TLatex *tex = new TLatex();
    FormatLatex(tex);
    tex->DrawLatexNDC(0.6, 0.8, "#it{#bf{CEPC}} Ref-TDR");
    tex->DrawLatexNDC(0.6, 0.74, "ZH, H#rightarrow #nu#nu ");
    tex->DrawLatexNDC(0.6, 0.68, "Z#rightarrow #mu#mu");
    // tex->DrawLatexNDC(0.75, 0.66, (TString)"Output: "+ (TString)to_string(evt_h1)+(TString)" events");
    tex->DrawLatexNDC(0.17,0.58,Form("mean = %.3f #pm %.3f",meanCB.getVal(),meanCB.getError()));
    tex->DrawLatexNDC(0.17,0.51,Form("#sigma = %.3f #pm %.3f",sigmaCB.getVal(),sigmaCB.getError()));
    tex->DrawLatexNDC(0.17,0.44,Form("Resolution = %.2f%%",100*sigmaCB.getVal()/meanCB.getVal()));
    tex->DrawLatexNDC(0.17,0.37,Form("#chi^{2}/ndf = %.2f",frame->chiSquare(6)));
    TPad *pad0 = (TPad*)c->GetPad(0);
    pad0->SetLeftMargin(0.15); 
    pad0->SetRightMargin(0.15);  
    pad0->SetBottomMargin(0.15); 
    pad0->SetTopMargin(0.1);  
    c->SaveAs("H_recoilMass.png");
    delete c;

    TCanvas* c1 = new TCanvas("c_roofit", "RooFit Gaussian Fit", 800, 600);
    RooRealVar meanCBZ("meanCBZ","",90.,80.,100.);
    RooRealVar sigmaCBZ("sigmaCBZ","",0.5,0.,5.);
    RooRealVar alphaCBLoZ("alphaCBLoZ","",1.,0.,5.);
    RooRealVar alphaCBHiZ("alphaCBHiZ","",1.,0.,10.);
    RooRealVar nCBLoZ("nCBLoZ","",2.,0.,5.);
    RooRealVar nCBHiZ("nCBHiZ","",2.,0.,10.);

    RooAbsPdf *pdfZ =new HggTwoSidedCBPdf("pdfZ", "pdf", invMassZ, meanCBZ, sigmaCBZ, alphaCBLoZ, nCBLoZ, alphaCBHiZ, nCBHiZ);
    RooFitResult *result1=pdfZ->fitTo(dataZ);
    RooPlot *frame1=invMassZ.frame(50);
    dataZ.plotOn(frame1);
    pdfZ->plotOn(frame1);
    frame1->GetYaxis()->SetTitle("Events");
    frame1->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");
    frame1->Draw();

    TLatex *texz = new TLatex();
    FormatLatex(texz);
    texz->DrawLatexNDC(0.6, 0.8, "#it{#bf{CEPC}} Ref-TDR");
    texz->DrawLatexNDC(0.6, 0.74, "Z#rightarrow #mu#mu");
    // tex->DrawLatexNDC(0.75, 0.66, (TString)"Output: "+ (TString)to_string(evt_h1)+(TString)" events");
    texz->DrawLatexNDC(0.17,0.58,Form("mean = %.3f #pm %.3f",meanCBZ.getVal(),meanCBZ.getError()));
    texz->DrawLatexNDC(0.17,0.51,Form("#sigma = %.3f #pm %.3f",sigmaCBZ.getVal(),sigmaCBZ.getError()));
    texz->DrawLatexNDC(0.17,0.44,Form("Resolution = %.2f%%",100*sigmaCBZ.getVal()/meanCBZ.getVal()));
    texz->DrawLatexNDC(0.17,0.37,Form("#chi^{2}/ndf = %.2f",frame->chiSquare(6)));
    TPad *pad = (TPad*)c1->GetPad(0);
    pad->SetLeftMargin(0.15); 
    pad->SetRightMargin(0.15);  
    pad->SetBottomMargin(0.15); 
    pad->SetTopMargin(0.1);   
    c1->Update();
    c1->Draw();
    c1->SaveAs("ZMass.png");
}
