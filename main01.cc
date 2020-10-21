// Headers and Namespaces.
#include "Pythia8/Pythia.h" // Include Pythia headers.
#include <iostream>
using namespace Pythia8;    // Let Pythia8:: be implicit.

int main(int argc, char *argv[])
{ // Begin main program.

  if (argv[1] == NULL)
  {
    std::cout << "Usage: ./main01.exe [pythia config file]\n";
    return 1;
  }

  Pythia pythia; // Declare Pythia object
  pythia.readFile(argv[1]);

  pythia.init(); // Initialize; incoming pp beams is default.

  pythia.next();

  pythia.stat();
  return 0;
} // End main program with error-free return.
