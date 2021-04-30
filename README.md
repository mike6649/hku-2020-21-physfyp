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
  
# 2020 Nov 18
- No progress, 0 misId rate.
- TODO:
  - Get ATLAS GEANT4 sample and retry

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
