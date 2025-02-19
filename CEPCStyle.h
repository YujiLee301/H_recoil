
#ifndef CEPC_STYLE_H
#define CEPC_STYLE_H
using namespace ROOT;

#include "TStyle.h"
#include "TVirtualPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TH1.h"
#include "TGraph.h"
/* 
CEPC Style guide 
Version: 1.1
Author:  Zhang Kaili
Mail: zhangkl@ihep.ac.cn
Date: 2024.10.02

To use this script, include and use SetCEPCStyle() before you plot. If you use several styles in the same time, you'd better know what you are doing. (Though harmless……)
Also, FormatLatex(), FormatLegend() should be used.

Some general standards for this style:
0. Root 6 recommended. May meet unexpected behavior in root 5.
1. 800*800px, No global title and stats/fit box, no x axis error bar if the bin is uniform; 
2. sans-serif font with absolute font size font(43), Axis title size 36, label size 34, caption size 28.
3. Fix canvas box position and use left margin 0.16. So if for Y axis the label is long, change the Y axis title offset, do not edit the margin or label size.
4. If you need transparent color, you may set CEPCStyle->SetCanvasPreferGL(kTRUE); But make sure your root support it. This option is defaultly on for Mac and lxslc6. See more in https://root.cern.ch/doc/master/classTColor.html
5. After use this style, you still manually set the following things to make your plot better:
 .1 Axis title. 						Do not center the title.
	For Y axis in 1d plot, the format should be like Entries or A.U. / XX [GeV]. 
	For X axis, use small m, m_H or m_Z for a specific particle. Use big M, for the combination. M_{#mu#mu}^{Recoil}. a spce between varible and unit [GeV].
 .2 Legend position and content.		All align left.
	Notice Draw option p,l and f and set the correct name so that root can link them.
	For the order, data points, usually "CEPC simulation / CEPC Ref-TDR" is the first one.
 .3 Caption position and content. 		All align left.
	Use TLatex and tex->DrawLatexNDC(0.8, 0.9, "#bf{CEPC 2018}");
	For TLatex and TMathText, use universal conventions in the script.
	In CEPC white paper now we use bold #bf{CEPC 2018}. In CEPC CDR we use #bf{CEPC CDR} in the first line.
	Then make clear the situation. #int L and #sqrt{s} is often ommited. So use like CEPC-v4, 5.6 ab^{-1}, or 240 GeV. If the plot is normalized then do not show the integral luminosity.
	Then make clear the process of the plot. 

Other tips:
1. Make sure the error bar in your plot is for the total CEPC lumi or just your statstics, esp. for those plots with scale factor.  (RooAbsData::Poisson or RooAbsData::SumW2)
2. Sometimes the first X axis label and Y axis label would overlap. Try frame->GetYaxis()->ChangeLabel(1, -1, 0); //after ROOT 6
3. Use #kern[-0.8] to fine tune some space.
4. Enjoy you plots!


CEPC style Change Log:
2024.10.02 v1.0: Change CEPC CDR Style to CEPC Style;
CEPC CDR style Change Log:
2019.01.06 v1.1: Turn to .h file
2018.08.12 v1.0: Initial 
 */

void SetCEPCStyle() {
	TStyle *CEPCStyle = new TStyle("CEPCStyle","Style for CEPC");

	// For the canvas:
	CEPCStyle->SetCanvasBorderMode(0);
	CEPCStyle->SetCanvasColor(kWhite);
	CEPCStyle->SetCanvasDefH(800); //Height of canvas
	CEPCStyle->SetCanvasDefW(800); //Width of canvas
	CEPCStyle->SetCanvasDefX(0);   //POsition on screen
	CEPCStyle->SetCanvasDefY(0);

	// For the Pad:
	CEPCStyle->SetPadBorderMode(0);
	// CEPCStyle->SetPadBorderSize(Width_t size = 1);
	CEPCStyle->SetPadColor(kWhite);
	CEPCStyle->SetGridStyle(3);
	CEPCStyle->SetGridWidth(1);

	// For the frame:
	CEPCStyle->SetFrameBorderMode(0);
	CEPCStyle->SetFrameBorderSize(1);
	CEPCStyle->SetFrameFillColor(0);
	CEPCStyle->SetFrameFillStyle(0);
	CEPCStyle->SetFrameLineColor(1);
	CEPCStyle->SetFrameLineStyle(1);
	CEPCStyle->SetFrameLineWidth(2);

	// For the histo:
	// CEPCStyle->SetHistFillColor(1);
	// CEPCStyle->SetHistFillStyle(0);
	// CEPCStyle->SetHistLineColor(1);
	CEPCStyle->SetHistLineStyle(0);
	// CEPCStyle->SetHistLineWidth(2);
	// CEPCStyle->SetLegoInnerR(0.8);
	// CEPCStyle->SetNumberContours(Int_t number = 20);

	CEPCStyle->SetEndErrorSize(2);
	//CEPCStyle->SetErrorMarker(20);
	//CEPCStyle->SetErrorX(0.);

	CEPCStyle->SetMarkerStyle(20);

	//For the fit/function:
	CEPCStyle->SetOptFit(0);
	CEPCStyle->SetFitFormat("5.4g");
	CEPCStyle->SetFuncColor(2);
	CEPCStyle->SetFuncStyle(1);
	CEPCStyle->SetFuncWidth(2);

	//For the date:
	CEPCStyle->SetOptDate(0);
	// CEPCStyle->SetDateX(Float_t x = 0.01);
	// CEPCStyle->SetDateY(Float_t y = 0.01);

	// For the statistics box:
	CEPCStyle->SetOptFile(0);
	CEPCStyle->SetOptStat(0); 		// To display the mean and RMS:   SetOptStat("mr");
	CEPCStyle->SetStatColor(kWhite);
	CEPCStyle->SetStatTextColor(1);
	CEPCStyle->SetStatFormat("6.4g");
	CEPCStyle->SetStatBorderSize(1);
	CEPCStyle->SetStatH(0.1);
	CEPCStyle->SetStatW(0.15);
	// CEPCStyle->SetStatStyle(Style_t style = 1001);
	// CEPCStyle->SetStatX(Float_t x = 0);
	// CEPCStyle->SetStatY(Float_t y = 0);

	// Margins:
    // canvas->SetMargin(0.16, 0.04, 0.11, 0.02); // left, right, bottom, top

	// CEPCStyle->SetPadLeftMargin(0.16);
	CEPCStyle->SetPadLeftMargin(0.12);
	CEPCStyle->SetPadRightMargin(0.04);
	CEPCStyle->SetPadBottomMargin(0.11);
	CEPCStyle->SetPadTopMargin(0.02);

	// For the Global title:

	CEPCStyle->SetOptTitle(0);
	CEPCStyle->SetTitleFont(43);
	CEPCStyle->SetTitleColor(1);
	CEPCStyle->SetTitleTextColor(1);
	CEPCStyle->SetTitleFillColor(0);
	CEPCStyle->SetTitleFontSize(36);


	// For the axis titles:

	CEPCStyle->SetTitleColor(1, "XYZ");
	CEPCStyle->SetTitleFont(43, "XYZ");
	CEPCStyle->SetTitleSize(36, "XYZ");

	CEPCStyle->SetTitleXOffset(1.0);
	CEPCStyle->SetTitleYOffset(2.0);
	// CEPCStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

	// For the axis labels:

	CEPCStyle->SetLabelColor(1, "XYZ");
	CEPCStyle->SetLabelFont(43, "XYZ");
	CEPCStyle->SetLabelOffset(0.005, "XYZ");
	CEPCStyle->SetLabelSize(34, "XYZ");

	// For the axis:

	CEPCStyle->SetAxisColor(1, "XYZ");
	CEPCStyle->SetStripDecimals(kTRUE);
	CEPCStyle->SetTickLength(0.02, "XYZ");
	CEPCStyle->SetNdivisions(508, "XYZ");    
	CEPCStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
	CEPCStyle->SetPadTickY(1);
	// CEPCStyle->SetPadTickZ(1);

	// Change for log plots:
	CEPCStyle->SetOptLogx(0); 
	CEPCStyle->SetOptLogy(0); // Change to 1 if need log(y).
	CEPCStyle->SetOptLogz(0);

	// Postscript options:
	CEPCStyle->SetPaperSize(20.,20.);
	// CEPCStyle->SetLineScalePS(Float_t scale = 3);
	// CEPCStyle->SetLineStyleString(Int_t i, const char* text);
	// CEPCStyle->SetHeaderPS(const char* header);
	// CEPCStyle->SetTitlePS(const char* pstitle);

	// CEPCStyle->SetBarOffset(Float_t baroff = 0.5);
	// CEPCStyle->SetBarWidth(Float_t barwidth = 0.5);
	// CEPCStyle->SetPaintTextFormat(const char* format = "g");
	// CEPCStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
	// CEPCStyle->SetTimeOffset(Double_t toffset);
	// CEPCStyle->SetHistMinimumZero(kTRUE);

	CEPCStyle->SetHatchesLineWidth(5);
	CEPCStyle->SetHatchesSpacing(0.05);
	
	// CEPCStyle->SetCanvasPreferGL(kTRUE); // ok in Mac or lxslc6, not in local windows/linux.

	CEPCStyle->cd();



}

void FormatLatex(TLatex *ch1)
{
	ch1->SetNDC(kTRUE);
	ch1->SetTextFont(43);
	ch1->SetTextSize(28);
	// ch1->SetTextAlign(31);
	// ch1->DrawLatex(X_position, Y_position, text.c_str());
}
void FormatAxis(TAxis * axis, double offset)
{
   axis->SetNdivisions(508);

   axis->SetLabelFont(43);
   axis->SetLabelSize(34);
   axis->SetLabelOffset(0.005);

   axis->SetTitleFont(43);
   axis->SetTitleColor(1);
   axis->SetTitleSize(36);
   axis->SetTitleOffset(offset);
//    axis->SetTickLength(0.02);
//    axis->CenterTitle();
}
void FormatAxis(TAxis * axis, double offset, TString title)
{
   axis->SetNdivisions(508);

   axis->SetLabelFont(43);
   axis->SetLabelSize(34);
   axis->SetLabelOffset(0.005);
//    axis->SetLabelOffset(1);

   axis->SetTitleFont(43);
   axis->SetTitleColor(1);
   axis->SetTitleSize(36);
   axis->SetTitleOffset(offset);
   axis->SetTitle(title);
//    axis->SetTickLength(0.02);
//    axis->CenterTitle();
}
void FormatLegend(TLegend *legend)
{
    legend->SetBorderSize(0);
    legend->SetTextFont(43);
    legend->SetTextSize(28);
    legend->SetFillStyle(0);
    legend->SetFillColor(0);
    legend->SetLineColor(0);
}

void FormatH1(TH1* H1,int n){
	H1->SetTitle("");
	// if(n==1) H1->SetFillColorAlpha(kBlue,0.35);
	H1->SetFillColor(n);

        // if(n==1) H1->SetFillColor(kCyan);
        // if(n==1) H1->SetLineColor(kCyan);
		//         if(n==1) H1->SetFillColor(kBlack);
        // if(n==1) H1->SetLineColor(kBlack);
        if(n==1) H1->SetFillStyle(1001);


}

void FormatG1(TGraph* g,int n){
int color[7]={1,2,4,6,8,21,41};
int style[7]={20,21,22,23,24,25,26};
g->SetMarkerSize(1.1);
g->SetLineWidth(1);
g->SetMarkerColor(color[n]);
g->SetLineColor(color[n]);
g->SetMarkerStyle(style[n]);
}

#endif
