#include "PhiTreeProof.h"

void PhiTreeProof::Init(TTree *tree)
{
    reader.SetTree(tree);
}

void PhiTreeProof::SlaveBegin(TTree *tree)
{
    h_pt = new TH1D("pt","pt",100,0,10);
    GetOutputList()->Add(h_pt);
}

bool PhiTreeProof::Process(Long64_t entry)
{
    reader.SetEntry(entry);

    for(int i=0; i<pt.GetSize(); ++i)
    {
        h_pt->Fill(pt[i]);
    }

    return true;

}

void PhiTreeProof::Terminate()
{
    h_pt->Draw();
}
