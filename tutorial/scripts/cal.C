void cal()
{
	// Constants
	const Int_t nBinTheta = 125;

	// Open file
	TFile* file = new TFile("tutorial_hist.root", "update");

	TH1D* hLeft  = (TH1D*) file -> Get("hLeft" );
	TH1D* hRight = (TH1D*) file -> Get("hRight");
	TH1D* hAsym  = new TH1D("hAsym", "hAsym", nBinTheta, 0, 25);

	// Asymmetry
	Double_t left, right;
	Double_t asym, sigma;
	for ( Int_t i = 0; i < nBinTheta; i++ )
	{
		left  = hLeft  -> GetBinContent(i+1);
		right = hRight -> GetBinContent(i+1);
		if ( left + right <= 0 ) continue;

		asym = (left - right) / (left + right);
		sigma = TMath::Sqrt( 4. * left * right / TMath::Power(left + right, 3) );

		hAsym -> SetBinContent(i+1, asym);
		hAsym -> SetBinError(i+1, sigma);
	}

	hAsym -> Write();

	// Bye bye
	file -> Close();
	delete file;

	return;
}
