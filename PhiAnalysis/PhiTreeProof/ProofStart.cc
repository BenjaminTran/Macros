void ProofStart()
{
    TProof::Open("");
    //TChain* chain = new TChain("PhiTree/TrackTree");
    //chain->AddFile("/Users/btran/research/RootFiles/Phi/Phi_TrackTree_PD1_v3_complete_Mar042018.root");
    TChain* chain = new TChain("PhiGenMatch/SignalTree");
    chain->AddFile("/Users/btran/research/TestRootFiles/PhiGenMatch_v3_1.root");
    chain->SetProof();
    chain->Process("PhiTreeProof.cc+","");
}
