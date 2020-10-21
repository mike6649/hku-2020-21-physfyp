// Headers and Namespaces.
//http://home.thep.lu.se/~torbjorn/pythia83html/ParticleDataScheme.html#section0
#include "Pythia8/Pythia.h" // Include Pythia headers.
#include <iostream>
using namespace Pythia8; // Let Pythia8:: be implicit.

void switch_on_z_ee_only(Pythia& pythia) {
  pythia.readString("23:onMode = 0");
  pythia.readString("23:onIfAll = 11 -11"); // only Z -> ee processes
  return;
}

int main(int argc, char *argv[])
{                // Begin main program.
                 // Set up generation.
  Pythia pythia; // Declare Pythia object
  // switch_on_z_ee_only(pythia);
  int id;
  if (argv[1]!=NULL) {
    id = atoi(argv[1]);
  } else {
    id = 23;
    std::cout << "Usage: ./main02.exe [particleID]\n";
    std::cout << "Defaulting to 23 (Z0)\n\n";
  }
  ParticleDataEntry *entry = pythia.particleData.particleDataEntryPtr(id);
  std::cout << "Name: " << entry->name() << "\n";
  std::cout << "NUM OF DECAY CHANNELS: " << entry->sizeChannels() <<"\n";
  for (int i = 0; i < entry->sizeChannels(); i++)
  {
    DecayChannel &decayChannel = entry->channel(i);
    int numProducts = decayChannel.multiplicity();
    std::cout << "CHANNEL" << i << ", onMode: " << decayChannel.onMode() << "-- " << numProducts << " products.\n";
    for (int j = 0; j < numProducts; j++)
    {
      int productId = decayChannel.product(j);
      std::string name = pythia.particleData.findParticle(productId)->name();
      std::cout << productId << " " << name << endl;
    }
    std::cout << endl;
  }
  return 0;
}