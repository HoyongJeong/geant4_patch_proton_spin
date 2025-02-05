void drawTest()
{
	// Constants
	const Int_t nBinEnergy = 10000;
	const Int_t nBinTheta  = 125;
	const Int_t nBinPos    = 8000;
	const Double_t zDet    = 898.5; // in millimeter

	// Open file
	TFile* file = new TFile("tutorial_detect.root");

	// Get tree
	TTree* tree = (TTree*) file -> Get("hits");
	Int_t nEntries = tree -> GetEntries();
	cout << nEntries << " entries detected" << endl;

	// Branches
	Int_t trackID, parID, parCharge;
	Double_t kinEgyTar, posXDet, posYDet;
	tree -> SetBranchAddress("trackID"  , &trackID  );
	tree -> SetBranchAddress("parID"    , &parID    );
	tree -> SetBranchAddress("parCharge", &parCharge);
	tree -> SetBranchAddress("kinEgyTar", &kinEgyTar);
	tree -> SetBranchAddress("posXDet"  , &posXDet  );
	tree -> SetBranchAddress("posYDet"  , &posYDet  );

	// Hists
	TH1D* hEnergy = new TH1D("hEnergy", "hEnergy", nBinEnergy, 0, 250);
	TH1D* hAll    = new TH1D("hAll"   , "hAll"   , nBinTheta , 0 , 25);
	TH1D* hLeft   = new TH1D("hLeft"  , "hLeft"  , nBinTheta , 0 , 25);
	TH1D* hRight  = new TH1D("hRight" , "hRight" , nBinTheta , 0 , 25);
	TH2D* hMap    = new TH2D("hMap", "hMap", nBinPos, -400, 400, nBinPos, -400, 400);

	// Loop
	Double_t r, theta, phi;
	for (Int_t i = 0; i < nEntries; i++ )
	{
		tree -> GetEntry(i);

		// Cut
		if ( trackID != 1 ) continue; // primary protons only

		// Cal
		r = TMath::Sqrt(posXDet * posXDet + posYDet * posYDet);
		theta = TMath::ATan2(r, zDet) * TMath::RadToDeg();
		phi = TMath::ATan2(posYDet, posXDet) * TMath::RadToDeg();
		cout << phi << endl;

		// Fill hists
		hEnergy -> Fill(kinEgyTar);
		hAll    -> Fill(theta);
		if ( phi > - 90.0 && phi < 90.0 ) hLeft  -> Fill(theta);
		else                              hRight -> Fill(theta);
		hMap    -> Fill(posXDet, posYDet);

		if ( (i+1) % 100000 == 0 ) cout << (i+1) << "-th event processing..." << endl;
	}

	// Write
	TFile* file_out = new TFile("tutorial_hist.root", "recreate");
	hEnergy -> Write();
	hAll    -> Write();
	hLeft   -> Write();
	hRight  -> Write();
	hMap    -> Write();

	// Bye bye
	file_out -> Close();

	delete file_out;
	delete hEnergy;
	delete hAll;
	delete hLeft;
	delete hRight;
	delete hMap;
	delete file;

	return;
}
