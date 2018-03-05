void ProofSparseStart()
{
    TProof::Open("");
    TChain* chain = new TChain("PhiTree/TrackTree");
    chain->AddFile("/Users/btran/research/RootFiles/Phi/Phi_TrackTree_PD1_v3_files1-80_022518_JL21.root");
    chain->SetProof();
    chain->Process("PhiTreeProofSparse.cc+");
}
