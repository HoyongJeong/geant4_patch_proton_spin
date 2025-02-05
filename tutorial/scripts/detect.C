void detect(const char* fIn = "", const char* fOut = "")
{
	// Open file
	TFile* file_in = new TFile(fIn);
	if ( ! file_in -> IsOpen() )
	{
		cerr << "[Error] Something wrong with opening file" << fIn << endl;
		cerr << "Exiting..." << endl;
		return;
	}

	// Get tree
	TTree* tree = (TTree*) file_in -> Get("hits");
	cout << tree -> GetEntries() << " entries detected" << endl;

	// File to be written
	TFile* file_out = new TFile(fOut, "RECREATE");
	if ( ! file_out -> IsOpen() )
	{
		cerr << "[Error] Something wrong with opening file" << fOut << endl;
		cerr << "Exiting..." << endl;
		return;
	}

	// Detection
	TCut detection = "parCharge != 0 && posXDet != 0 && posYDet != 0";
	TTree* tree_detected = tree -> CopyTree(detection);
	cout << tree_detected -> GetEntries() << " events detected" << endl;

	// Bye bye
	file_out -> Write();
	file_out -> Close();
	file_in  -> Close();

	delete file_out;
	delete file_in;

	return;
}
