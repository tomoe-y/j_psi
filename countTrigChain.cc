void countTrigChain(){
	TChain *chain = new TChain("physics");
	//chain->Add("/gpfs/fs2001/junpei/L1TGCNtuple/data16_13TeV/data16_13TeV.periodL/user.junpei.00311481.physics_Main.merge.NTUP_MCP.1.f758_m1714.00-00-32_L1TGCNtuple_derivated.02-03-00.root");
	chain->Add("/home/toyamash/create_Ntuple/run/L1TGCNtuple.root");

	chain->SetBranchStatus("*", 0);
	chain->SetBranchStatus("trigger_info_chain", 1);

	std::vector<std::string> *trigger_info_chain = 0;

	chain->SetBranchAddress("trigger_info_chain", &trigger_info_chain);

	int entries = chain->GetEntries();
	int counts = 0;

	for (int i = 0; i < entries; i++){
		chain->GetEntry(i);
		for(int j = 0; j < trigger_info_chain->size(); j++){
			//if (trigger_info_chain->at(j) == "HLT_mu20_2mu0noL1_JpsimumuFS") counts++;
			if (trigger_info_chain->at(j) == "HLT_mu24_ivarmedium_mu4_probe_L1MU14FCH") counts++;
		}
	}
	std::cout << counts << std::endl;
}

