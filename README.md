Roadmap for the project on:
# Measurement of the chance of misidentification of electrons in Accelerator event data

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
