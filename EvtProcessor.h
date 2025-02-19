#ifndef EvtProcessor_h
#define EvtProcessor_h

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TLorentzVector.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"

class EvtProcessor{
    public:
      void initialization();
      void GenLepSelection();
      void Genmatching();
      void PFOLepSelection();
      void GoodLepSelection();
      bool MatchFlag(TLorentzVector A, TLorentzVector B);
      int Finalstate; //11 or 13
      int PFOsize;
      int MCsize;
      float rec_MZ; float rec_MH;
      float Getrecoilmass(TLorentzVector l1, TLorentzVector l2);
      float Getchisquare(TLorentzVector l1, TLorentzVector l2);
      void GetZHMass();
      std::vector<int> GenLepindex;
      std::vector<TLorentzVector> GenLep;
      std::vector<TLorentzVector> PFOLep;
      std::vector<int> PFOLepindex;
      std::vector<int> SelectedLepindex;
      std::vector<TLorentzVector> SelectedLep;
      void findZleps();
      std::vector<int> ZLepindex;
      std::vector<TLorentzVector> ZLep;
      void setPFOinfo(std::vector<Float_t> PFO_E_,std::vector<Float_t> PFO_C_,  std::vector<Float_t> PFO_x_,
                      std::vector<Float_t> PFO_y_,std::vector<Float_t> PFO_z_, int PFOsize_);
      void setMCtruthinfo(std::vector<double> Truth_m_, std::vector<int> Truth_PDG_, std::vector<int> Truth_status_,
                         std::vector<Float_t> Truth_x_, std::vector<Float_t> Truth_y_, std::vector<Float_t> Truth_z_, int MCsize_);
      EvtProcessor(TString fs);      
      float sigmaH = 0.276;
      float sigmaZ = 1.183;
    
      float MZ = 91.3;
      float MH = 125.38;
      std::vector<Float_t> PFO_E;
      std::vector<Float_t> PFO_C;
      std::vector<Float_t> PFO_x;
      std::vector<Float_t> PFO_y;
      std::vector<Float_t> PFO_z;
      std::vector<double> Truth_m;
      std::vector<int> Truth_PDG;
      std::vector<int> Truth_status;
      std::vector<Float_t> Truth_x;
      std::vector<Float_t> Truth_y;
      std::vector<Float_t> Truth_z;
};
void EvtProcessor::initialization(){
    PFOsize=0;MCsize=0;rec_MH=0;rec_MZ=0;
    PFO_C.clear();PFO_E.clear();PFO_x.clear();PFO_y.clear();PFO_z.clear();
    Truth_m.clear();Truth_PDG.clear();Truth_status.clear();Truth_x.clear();Truth_y.clear();Truth_z.clear();
    GenLepindex.clear();GenLep.clear();PFOLepindex.clear();SelectedLepindex.clear();PFOLep.clear();SelectedLep.clear();
    ZLep.clear();ZLepindex.clear();
}
EvtProcessor::EvtProcessor(TString fs){
    if (fs=="ee") Finalstate=11;
    if (fs=="mm") Finalstate=13;
    else Finalstate=11;

    if(Finalstate==11){
      float sigmaH = 0.269;
      float sigmaZ = 1.231;     
    }
}
void EvtProcessor::setPFOinfo(std::vector<Float_t> PFO_E_,std::vector<Float_t> PFO_C_,  std::vector<Float_t> PFO_x_,
                      std::vector<Float_t> PFO_y_,std::vector<Float_t> PFO_z_, int PFOsize_){
    PFOsize=PFOsize_;                    
    for (int k=0;k<PFOsize;k++){
        PFO_E.push_back(PFO_E_[k]);
        PFO_C.push_back(PFO_C_[k]);
        PFO_x.push_back(PFO_x_[k]);
        PFO_y.push_back(PFO_y_[k]);
        PFO_z.push_back(PFO_z_[k]);
    }
}
void EvtProcessor::setMCtruthinfo(std::vector<double> Truth_m_, std::vector<int> Truth_PDG_, std::vector<int> Truth_status_,
                         std::vector<Float_t> Truth_x_, std::vector<Float_t> Truth_y_, std::vector<Float_t> Truth_z_, int MCsize_){
    MCsize=MCsize_;
    for (int k=0;k<MCsize;k++){
        Truth_m.push_back(Truth_m_[k]);
        Truth_PDG.push_back(Truth_PDG_[k]);
        Truth_status.push_back(Truth_status_[k]);
        Truth_x.push_back(Truth_x_[k]);
        Truth_y.push_back(Truth_y_[k]);
        Truth_z.push_back(Truth_z_[k]);
    }
}
void EvtProcessor::GenLepSelection(){
    for(int k=0; k<MCsize; k++){
        if(Truth_status[k]!=1) continue;
        if(abs(Truth_PDG[k])==Finalstate) GenLepindex.push_back(k);
    }
    for (int j=0;j<GenLepindex.size();j++){
        TLorentzVector Lep;
        Lep.SetPxPyPzE(Truth_x[GenLepindex[j]],Truth_y[GenLepindex[j]],Truth_z[GenLepindex[j]],sqrt(Truth_x[GenLepindex[j]]*Truth_x[GenLepindex[j]]+Truth_y[GenLepindex[j]]*Truth_y[GenLepindex[j]]+Truth_z[GenLepindex[j]]*Truth_z[GenLepindex[j]]));
        GenLep.push_back(Lep);
    }
}
bool EvtProcessor::MatchFlag(TLorentzVector A, TLorentzVector B){
    bool matched=false;
    if (A.DeltaR(B, true)<0.1) matched=true;
    return matched;
}
void EvtProcessor::PFOLepSelection(){
    for(int m=0;m<PFOsize;m++){
        TLorentzVector Lep;
        bool Matched=false;
        Lep.SetPxPyPzE(PFO_x[m],PFO_y[m],PFO_z[m],PFO_E[m]);
        for(int n=0;n<GenLep.size();n++){
            if (MatchFlag(GenLep[n],Lep)) Matched=true;
        }
        if(Matched){
          PFOLepindex.push_back(m);
          PFOLep.push_back(Lep);
        } 
    }
}
float EvtProcessor::Getrecoilmass(TLorentzVector l1, TLorentzVector l2){
    return sqrt((240-(l1+l2).E())*(240-(l1+l2).E())-((l1+l2).Px()*(l1+l2).Px()+(l1+l2).Py()*(l1+l2).Py()+(l1+l2).Pz()*(l1+l2).Pz()));
}
void EvtProcessor::GoodLepSelection(){
    for (int i=0;i<PFOLep.size();i++){
        if(PFOLep[i].Pt()<5) continue;
        if(cos(PFOLep[i].Theta())>0.85) continue;
        SelectedLepindex.push_back(PFOLepindex[i]);
        SelectedLep.push_back(PFOLep[i]);
    }
}
float EvtProcessor::Getchisquare(TLorentzVector l1, TLorentzVector l2){
    float recoilmass = Getrecoilmass(l1, l2);
    return (((recoilmass-MH)*(recoilmass-MH))/(sigmaH*sigmaH)+(((l1+l2).M()-MZ)*((l1+l2).M()-MZ))/(sigmaZ*sigmaZ));
}
void EvtProcessor::findZleps(){
    if (SelectedLepindex.size()<2) return;
    bool foundZ =false;
    float minimalchisquare = 999999;
    int lep0id = 0;int lep1id = 1;
    for(int i=0;i<SelectedLepindex.size()-1;i++){
        for(int j=i+1;j<SelectedLepindex.size();j++){
            if((PFO_C[SelectedLepindex[i]]+PFO_C[SelectedLepindex[j]])!=0) continue;
            float chisquare = Getchisquare(SelectedLep[i],SelectedLep[j]);
            if(chisquare<minimalchisquare){
                lep0id = i; lep1id = j;
                minimalchisquare = chisquare;
                foundZ =true;
            }
        }
    }
    if (!foundZ) return;
    ZLepindex.push_back(SelectedLepindex[lep0id]);
    ZLepindex.push_back(SelectedLepindex[lep1id]);
    ZLep.push_back(SelectedLep[lep0id]);
    ZLep.push_back(SelectedLep[lep1id]);
}
void EvtProcessor::GetZHMass(){
    if (ZLep.size()<2) return;
    rec_MZ = (ZLep[0]+ZLep[1]).M();
    rec_MH = Getrecoilmass(ZLep[0],ZLep[1]);
}
#endif
