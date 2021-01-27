# Build
1. Install ROOT and HepMC3
2. From here either use MadGraph or Pythia for MC event generation
  - Pythia
    - [Install Pythia and Delphes](https://cp3.irmp.ucl.ac.be/projects/delphes/wiki/WorkBook/Pythia8)
    - `. $DELPHES/DephesPythia8 $DELPHES/cards/delphes_card_ATLAS.tcl configZtoee.cmnd delphes_output.root`
  - MadGraph
    - Install MadGraph
    - Generate the process `p p > Z , Z > e+ e-` and output to HepMC format
    - `. $DELPHES/DephesHepMC $DELPHES/cards/delphes_card_ATLAS.tcl delphes_output.root output.hepmc`
4. `root -l`
5. `gROOT->ProcessLine(".include /path/to/delphes/");`
`gROOT->ProcessLine(".include /path/to/delphes/external/");`
6. `gSystem->Load("/path/to/delphes/libDelphes");`
7. `.x diElectronMass.C`

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

# 2021 Jan 27 Meeting Notes
- Done
  - Calculate mis-id rate with `misIdCount2.C`
- TODO:
  - Implement likelihood method

# 2020 Nov 11 Meeting Notes
- Done
  - Generate ROOT file (10000 events)
  - Graph the Z invariant Mass with `diElectronMass.C`
  - *0* mis-identified electrons (expected ~10)
  - Verify eta distribution cuts off at 2.5
- TODO:
  - Try Delphes again with 100k events
  - Check parents of all electrons are from Z
  - Discuss further if *0* misId rate again

# 2020 Nov 04 Meeting Notes 
- Done
  - Using Madgraph & Delphes, generate ROOT file
  - linking Delphes libraries and headers with ROOT
- TODO:
  - Use ROOT TLorentzVector::M (?) to find the invariant mass

# 2020 Oct 28 Meeting Notes
- No progress, Delphes rejected Pythia HepMC format
- TODO:
  - Further troubleshooting together
  
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

# 2020 Oct 16 Meeting Notes
- TODO Generate the sample events before next week
- In Pythia, turn on detector simulation
  - Only turn on electroweak processes
  - pp -> Z -> ee
- The information about electron misidentification is stored in the sample
- After succeed
  - Use ROOT to perform further analysis