void drawIniConProfile(TString funcName="gaus2dwhotspot",double b = 4.0,double r=7.5) {
 
  gROOT->ProcessLine(Form(".L %s.cpp",funcName.Data()));
  gROOT->ProcessLine(".L loadStyle.C");
  loadStyle();

  TF2* fIniCon = new TF2(funcName.Data(),myinicon,-10.,10.,-10.,10.,5);
  fIniCon->SetTitle("");
  double par[] = {1.,0.,0.,(r-b/2.)/3.,r/3.};
  int npar = sizeof(par)/sizeof(double);
  for(int i=0;i<npar;++i) fIniCon->SetParameter(i,par[i]);

  TCanvas* c1 = new TCanvas("c1","Initial Energy Density Transverse Profile",10,10,600,600);
  c1->SetGridx();
  c1->SetGridy();
  fIniCon->SetNpx(70);
  fIniCon->SetNpy(70);
  fIniCon->SetLineColor(1);
  fIniCon->SetLineWidth(1);
  fIniCon->GetXaxis()->SetTitle("x [fm]");
  fIniCon->GetXaxis()->SetTitleOffset(2.0);
  fIniCon->GetXaxis()->CenterTitle(true);
  fIniCon->GetXaxis()->SetNdivisions(305);
  fIniCon->GetYaxis()->SetTitle("y [fm]");
  fIniCon->GetYaxis()->SetTitleOffset(2.0);
  fIniCon->GetYaxis()->CenterTitle(true);
  fIniCon->GetYaxis()->SetNdivisions(305);
  fIniCon->GetZaxis()->SetNdivisions(310);
  fIniCon->Draw("contz");

  TCanvas* c2 = new TCanvas("c2","Initial Energy Density Transverse Profile",600,10,600,600);
  fIniCon->Draw("surf1fb");

}

inline double myinicon(double *x,double *par) { return inicon(x,2,par); }
