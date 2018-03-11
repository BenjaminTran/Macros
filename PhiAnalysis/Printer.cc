#include "TFile.h"
#include "TString.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TRegexp.h"

void Printer(TFile* file, TString reg, std::string Destination)
{
    TRegexp re(reg,kTRUE);
    for(const auto&& key : *file->GetListOfKeys())
    {
        TString st = key->GetName();
        if(st.Index(re) == kNPOS) continue;
        std::string obj_name = key->GetName();
        TCanvas* c;
        file->GetObject(obj_name.c_str(),c);
        std::string printedName = Destination + "/" + obj_name;
        std::cout << "working" << endl;
        c->Print(printedName.c_str());
        std::cout << "not working" << endl;
    }
}
