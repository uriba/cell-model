from sympy import *
from math import log

#Constants:
RibosomeMass = 2500000 #BIND
RNAPolymeraseMass = 500000 #BIND
AverageAAMass = 110
AverageNucleotideMass = 325
ProteinFractionOfRibosome = 1.0/3
RNAFractionOfRibosome = 2.0/3
RibosomalProteinAA = (RibosomeMass * 1/3) / AverageAAMass

#Translation rate in aa per second:
TranslationRate = 20

#Transcription rate in nucletotides per second
TranscriptionRate = 80

#Cellular components are measured as fractions of the biomass.
Ribosomes = Symbol('Ribosomes',positive=True)
RibosomalProtein = Ribosomes*ProteinFractionOfRibosome
RibosomalRNA = Ribosomes*RNAFractionOfRibosome

RNAPolymerase = Symbol('RNAPolymerase',positive=True)

#Growth rate is in the natural base, e.
GrowthRate = Symbol('GrowthRate')

#Total protein is the sum of all protein components:
Protein=RibosomalProtein+RNAPolymerase

#Total RNA is ribosomal RNA:
RNA=RibosomalRNA

#AA per ribosome in proteome
PerRibosomeAA = (Protein/RibosomalProtein) * RibosomalProteinAA

#A constant to connect the mass of RNA to the mass of RNAPolymerase.
PerRNAPolymeraseRNAMass = RNAPolymeraseMass * RNA / RNAPolymerase
PerRNAPNucleotides = PerRNAPolymeraseRNAMass / AverageNucleotideMass

#Protein growth rate:
ProteinDoubling = PerRibosomeAA/TranslationRate

#RNA growth rate:
RNADoubling = PerRNAPNucleotides/TranscriptionRate

#Under balanced growth the two are equal:
Growth = Eq(RNADoubling,ProteinDoubling)

print Growth
[sol] = solve([Growth,Eq(Ribosomes+RNAPolymerase,1)])
results = {}
results['RNAPolymerase'] = N(sol[RNAPolymerase])
results['Ribosomes'] = N(sol[Ribosomes])
results['RibosomalProtein'] = ProteinFractionOfRibosome*results['Ribosomes']
results['RibosomalRNA'] = RNAFractionOfRibosome*results['Ribosomes']
results['Protein'] = results['RibosomalProtein'] + results['RNAPolymerase']
results['RNA'] = results['RibosomalRNA']
results['PerRibosomeAA'] = (results['Protein']/results['RibosomalProtein']) * RibosomalProteinAA
results['PerRNAPolymeraseRNAMass'] = RNAPolymeraseMass * results['RNA'] / results['RNAPolymerase']
results['PerRNAPNucleotides'] = results['PerRNAPolymeraseRNAMass'] / AverageNucleotideMass
results['ProteinDoubling'] = results['PerRibosomeAA']/TranslationRate
results['RNADoubling'] = results['PerRNAPNucleotides']/TranscriptionRate

print results
