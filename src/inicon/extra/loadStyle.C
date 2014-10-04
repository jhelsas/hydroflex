loadStyle() {
  gStyle->SetPalette(1);
  //  gStyle->SetCanvasPreferGL(true);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadTopMargin(0.15);
  gStyle->SetStatColor(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);

  gStyle->SetPalette(1, 0);
  const int NCont = 99; // should suffice, can be set to 255
  const int NRGBs = 5;
  double stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  double red[NRGBs] = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  double green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  double blue[NRGBs] = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  //
  // USE THIS FOR GRAYSCALE PLOTS
  //const int NCont = 25; // should suffice, can be set to 255
  //const int NRGBs = 7;
  //double stops[NRGBs] = { 0.00, 0.10, 0.25, 0.45, 0.60, 0.75, 1.00 };
  //double red[NRGBs] = { 1.00, 0.00, 0.00, 0.00, 0.97, 0.97, 0.10 };
  //double green[NRGBs] = { 1.00, 0.97, 0.30, 0.40, 0.97, 0.00, 0.00 };
  //double blue[NRGBs] = { 1.00, 0.97, 0.97, 0.00, 0.00, 0.00, 0.00 };
  //
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

}
