from sympy import *
from math import log

#Constants:
#Weights:
RibosomeMass              = 2700000.0 #BIND100118
RNAPolymeraseMass         =  500000.0 #BIND104925
#DNAPolymeraseMass         =  790000.0 #BIND104931
tRNAMass                  = 25000.0 #BIND101177
AverageAAMass             = 110.0
AverageNucleotideMass     = 325.0

#Ratios
ProteinFractionOfRibosome = 1.0/3
RNAFractionOfRibosome     = 2.0/3
RibosomalProteinAA        = (RibosomeMass * 1/3) / AverageAAMass
tRNAPerRibosome           = 5.0 #Rough minimal value, assuming each ribosome has 3 bound tRNA's, one being charged, and another one floating around.

#Rates
#Translation rate in aa per second:
TranslationRate = 20.0

#Transcription rate in nucletotides per second
TranscriptionRate = 80.0

#DNA replication rate in nucletotides per second
#DNAReplicationRate = 600.0 #BIND109251 

#Cellular components are measured as fractions of the biomass.
Ribosomes = Symbol('Ribosomes',positive=True)
RibosomalProtein = Ribosomes*ProteinFractionOfRibosome
RibosomalRNA = Ribosomes*RNAFractionOfRibosome

tRNAs = (Ribosomes/RibosomeMass)*tRNAPerRibosome*tRNAMass

RNAPolymerase = Symbol('RNAPolymerase',positive=True)
#DNAPolymerase = Symbol('DNAPolymerase',positive=True)

#Growth rate is in the natural base, e.
GrowthRate = Symbol('GrowthRate')

#Total protein is the sum of all protein components:
Protein=RibosomalProtein+RNAPolymerase #+DNAPolymerase

#Total RNA is ribosomal RNA:
RNA=RibosomalRNA + tRNAs

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
[sol] = solve([Growth,Eq(Protein+RNA,1)])
results = {}
results['RNAPolymerase'] = N(sol[RNAPolymerase])
results['Ribosomes'] = N(sol[Ribosomes])
results['tRNAs'] = (N(sol[Ribosomes])/RibosomeMass)*tRNAPerRibosome*tRNAMass
results['RibosomalProtein'] = ProteinFractionOfRibosome*results['Ribosomes']
results['RibosomalRNA'] = RNAFractionOfRibosome*results['Ribosomes']
results['Protein'] = results['RibosomalProtein'] + results['RNAPolymerase']
results['RNA'] = results['RibosomalRNA'] + results['tRNAs']
results['PerRibosomeAA'] = (results['Protein']/results['RibosomalProtein']) * RibosomalProteinAA
results['PerRNAPolymeraseRNAMass'] = RNAPolymeraseMass * results['RNA'] / results['RNAPolymerase']
results['PerRNAPNucleotides'] = results['PerRNAPolymeraseRNAMass'] / AverageNucleotideMass
results['ProteinDoubling'] = results['PerRibosomeAA']/TranslationRate
results['RNADoubling'] = results['PerRNAPNucleotides']/TranscriptionRate

for key,value in results.iteritems():
    print "%s\t:%f" % (key,value)
