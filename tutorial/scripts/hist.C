void hist(const char* fIn = "", const char* fOut = "")
{
	// Constants
	const Int_t nBinEnergy = 10000;
	const Int_t nBinTheta  = 125;
	const Int_t nBinPos    = 1000;
	const Double_t zDet    = 898.5; // in millimeter

	// Open file
	TFile* file_in = new TFile(fIn);
	if ( ! file_in -> IsOpen() )
	{
		cerr << "[Error] Something wrong with opening file" << fIn << endl;
		cerr << "Exiting..." << endl;
		return;
	}
	else
	{
		cout << "[Notice] Opening " << fIn << " as an input" << endl;
	}

	// Get tree
	TTree* tree = (TTree*) file_in -> Get("hits");
	Int_t nEntries = tree -> GetEntries();
	cout << nEntries << " entries detected" << endl;

	// Branches
	Int_t trackID, parID, parCharge;
	Double_t kinEgyTar, posXDet, posYDet;
	Int_t hadElastic, hIoni, CoulombScat, protonInelastic, etc;
	tree -> SetBranchAddress("trackID"        , &trackID        );
	tree -> SetBranchAddress("parID"          , &parID          );
	tree -> SetBranchAddress("parCharge"      , &parCharge      );
	tree -> SetBranchAddress("kinEgyTar"      , &kinEgyTar      );
	tree -> SetBranchAddress("posXDet"        , &posXDet        );
	tree -> SetBranchAddress("posYDet"        , &posYDet        );
	tree -> SetBranchAddress("hadElastic"     , &hadElastic     );
	tree -> SetBranchAddress("hIoni"          , &hIoni          );
	tree -> SetBranchAddress("CoulombScat"    , &CoulombScat    );
	tree -> SetBranchAddress("protonInelastic", &protonInelastic);
	tree -> SetBranchAddress("etc"            , &etc            );

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

		// Only primary protons which undergo hadElastic once
		if ( trackID != 1 || hadElastic != 1 || hIoni != 0 || CoulombScat != 0 || protonInelastic != 0 || etc != 0 ) continue;

		// Cal
		r = TMath::Sqrt(posXDet * posXDet + posYDet * posYDet);
		theta = TMath::ATan2(r, zDet) * TMath::RadToDeg();
		phi = TMath::ATan2(posYDet, posXDet) * TMath::RadToDeg();

		// Fill hists
		hEnergy -> Fill(kinEgyTar);
		hAll    -> Fill(theta);
		if ( phi > - 90.0 && phi < 90.0 ) hLeft  -> Fill(theta);
		else                              hRight -> Fill(theta);
		hMap    -> Fill(posXDet, posYDet);

		if ( (i+1) % 100000 == 0 ) cout << (i+1) << "-th event processing..." << endl;
	}

	// Write
	TFile* file_out = new TFile(fOut, "RECREATE");
	if ( ! file_out -> IsOpen() )
	{
		cerr << "[Error] Something wrong with opening file" << fOut << endl;
		cerr << "Exiting..." << endl;
		return;
	}
	else
	{
		cout << "[Notice] Opening " << fOut << " as an output" << endl;
	}
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
	delete file_in;

	return;
}
