void ProofStart()
{
    TProof::Open("");
    TChain* chain = new TChain("PhiTree/TrackTree");
    chain->AddFile("/Users/btran/research/RootFiles/Phi/Phi_TrackTree_PD1_v3_files1-80_022518_JL21.root");
    chain->AddFile("/Users/btran/research/RootFiles/Phi/Phi_TrackTree_PD1_v3_files81-_022518_JL21.root");
    chain->SetProof();
    chain->Process("PhiTreeProof.cc+");
    //chain->Process("PhiTreeProof1.cc+");
    //chain->Process("PhiTreeProof2.cc+");
    //chain->Process("PhiTreeProof3.cc+");
    //chain->Process("PhiTreeProof4.cc+");
    //chain->Process("PhiTreeProof5.cc+");
    //chain->Process("PhiTreeProof6.cc+");
    //chain->Process("PhiTreeProof7.cc+");
    //chain->Process("PhiTreeProof8.cc+");
    //chain->Process("PhiTreeProof9.cc+");
    //chain->Process("PhiTreeProof10.cc+");
    //chain->Process("PhiTreeProof11.cc+");
    //chain->Process("PhiTreeProof12.cc+");
}
