# Build and run pythia event generator:
- Copy over `PYTHIA_DIR/examples/Makefile.inc` into the main directory
- `make main01` then `./main01.exe`

# Measurement of the chance of misidentification of electrons in Accelerator event data
## Project Roadmap

- A brief rundown derivation of the steps in obtaining the charge flip rate ϵ<sub>i</sub>
  - Including the likelihood function 
  - Discuss about the related errors arose in the method used.

- Come up with a code to act as an analyzer.
  -	Work in Pythia & C++
  -	Input event data and output a weight factor of how likely misidentification occurs.

- Validation the code
  -	By using higher √s data
  -	Amend if necessary

# 2020 Oct 16 Meeting Notes
- TODO Generate the sample events before next week
- In Pythia, turn on detector simulation
  - Only turn on electroweak processes
  - pp -> Z -> ee
- The information about electron misidentification is stored in the sample
- After succeed
  - Use ROOT to perform further analysis

# 2020 Oct 21 Meeting Notes
- Done: Generate HepMC Events from Pythia
- TODO:
  - Install Delphes for detector simulation
    - Input .HepMC, output .ROOT
  - Pythia:
    - Further [specify allowed processes](http://home.thep.lu.se/~torbjorn/pythia82html/ElectroweakProcesses.html)
    - *Only* pure Z -> lepton lepton , no gamma
  - ROOT:
    - Get 4 momentum of electrons, to get the invariant mass of the Z Boson

# 2020 Oct 28 Meeting Notes
- No progress, Delphes rejected Pythia HepMC format
- TODO:
  - Further troubleshooting together

# 2020 Nov 04 Meeting Notes 
- Done
  - Using Madgraph & Delphes, generate ROOT file
  - linking Delphes libraries and headers with ROOT
- TODO:
  - Use ROOT TLorentzVector::M (?) to find the invariant mass
