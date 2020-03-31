// Headers and Namespaces.
#include "Pythia8/Pythia.h" // Include Pythia headers.
using namespace Pythia8;
int main(int argc, char* argv[]) {
// Set up generation.
Pythia pythia;
pythia.readFile(argv[1]);
pythia.init();
Hist pT("invariant mass", 1000, 0., 500000);
//Hist pT("top transverse momentum", 100, 0., 500000);
//Hist eta("top pseudorapidity", 100, -5., 5.);
for (int iEvent = 0; iEvent < 10000; ++iEvent)
{
    pythia.next();
    double E = 0;
    double px = 0;
    double py = 0;
    double pz = 0;
    int iex = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
    {
        if (pythia.event[i].id() == 4000001 || pythia.event[i].id() == 4000002)
        {
            iex = i;
        }
    }
    for (int i = 0; i < pythia.event.size(); ++i)
    {
        if (pythia.event[i].mother1() == iex)
        {
            E += pythia.event[i].e();
            px += pythia.event[i].px();
            py += pythia.event[i].py();
            pz += pythia.event[i].pz();
        }
        //cout << "i = " << i << ", id = " << pythia.event[i].id() << endl;
    }
    double im = 0;
    im = (E*E - px*px - py*py - pz*pz);
    pT.fill( im );
    //cout << pythia.event[iTop].charge();
}
cout << pT;
pythia.stat();
return 0;
}
// End main program with error-free return.
