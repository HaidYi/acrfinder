#!/usr/bin/python3
'''
	*****************************************************************************************************
	Purpose:	Uses candidate Acr/Aca loci and CDD's that infer pathogenicity to attempt to filter false positive Acr/Aca locus.

	Author:
		Javi Gomez - https://github.com/rtomyj wrote the version 1.0.
		Haidong Yi - https://github.com/haidyi revised the codes, remove some bugs, wrote the version 2.0. 
	Project:
		This work is advised by Dr. Yanbin Yin at UNL - yyin@unl.edu
	*****************************************************************************************************
'''
from collections import defaultdict

from sys import path as sys_path
sys_path.append('dependencies/PyGornism/')
from organism import Organism



'''
	Purpose:
		Uses the candidate Acr/Aca loci as a starting point to obtain neighboring proteins.
		Goes up stream from first Acr/Aca protein in locus with a value of PROTEIN_UP_DOWN.
		Goes down stream from last Acr/Aca protei in locus with a value of PROTEIN_UP_DOWN.
		Creates a dict containing locus number as a value and locus proteins of neighbors.
		Creates an faa string of the locus proteins of neighbors.
	Arguments:
		candidateAcrs - list of lists containing proteins in the inner list. Each protein list is an Acr/Aca locus.
		ORGANISM_SUBJECT - Organism object that contains information of the organism being parsed.
		PROTEIN_UP_DOWN - Number of proteins get upstream and downstream.
	Returns:
		neighborsFaaStr - faa representation of the newly acquired neighbors
		NEIGHBORHOOD_NUM_maps_NEIGHBORHOOD_WP - dict containing newly acquired neighbors (Protein objects) maped by the locus number
'''
def get_acr_neighbors(candidateAcrs, ORGANISM_SUBJECT, PROTEIN_UP_DOWN):
	neighborsFaaStr = ""
	NEIGHBORHOOD_NUM_maps_NEIGHBORHOOD_WP = defaultdict(list)

	'''
		Goes through the candidate Acr/Aca loci.
		Finds the neighboring proteins surrounding the locus.
		Obtains the FAA string representation of each neighbors.
	'''
	from formated_output import get_faa
	for neighborhoodNum, locus in enumerate(candidateAcrs):
		startProtein, endProtein = locus[0].wp, locus[len(locus) - 1].wp	# first and last protein of locus
		downstream = ORGANISM_SUBJECT.get_downstream_neighbors(PROTEIN_UP_DOWN, startProtein, inclusive=False)	# gets downstream neighbors
		upstream = ORGANISM_SUBJECT.get_upstream_neighbors(PROTEIN_UP_DOWN, endProtein)	# gets upstream neighbors


		'''
			Combines downstream, upstream and the remaining Acr/Aca proteins (neglecting start/end)
		'''
		neighborhood = downstream[:]
		neighborhood.extend(locus[1: len(locus) -1])
		neighborhood.extend(upstream[:])

		for protein in neighborhood:	# adds neighbor list to neighbor dict
			NEIGHBORHOOD_NUM_maps_NEIGHBORHOOD_WP[neighborhoodNum].append(protein.wp)

		neighborsFaaStr += '# Neighborhood for locus {0} starting with {1} and ending with {2}\n'.format(neighborhoodNum, startProtein, endProtein)
		neighborsFaaStr += get_faa(neighborhood)	# continues building faa string
		neighborsFaaStr += '\n\n'

	return neighborsFaaStr, NEIGHBORHOOD_NUM_maps_NEIGHBORHOOD_WP



'''
	Purpose:
		Uses rpsblast+ using CDD's to filter out candidate Acr/Aca loci
	Arguments:
		candidateAcrs - list of lists containing proteins in the inner list. Each protein list is an Acr/Aca locus
		ORGANISM_SUBJECT - Organism object that contains information of the organism being parsed
		NEIGHBORHOOD_FAA_PATH - Path to file containing FAA's of neighborhoods
		CDD_RESULTS_PATH - Where to store psiblast results
		CDD_DB_PATH - Path to CDD's being used
		PROTEIN_UP_DOWN - Number of proteins get upstream and downstream
		MIN_NUM_PROTEINS_MATCH_CDD - Minimum number of proteins that had CDD hits needed to keep a locus
	Returns:
		result from parse_cdd_results()
'''
def use_cdd(candidateAcrs, ORGANISM_SUBJECT, NEIGHBORHOOD_FAA_PATH, CDD_RESULTS_PATH, CDD_DB_PATH, PROTEIN_UP_DOWN, MIN_NUM_PROTEINS_MATCH_CDD):
	with open(NEIGHBORHOOD_FAA_PATH, 'w') as handle:
		neighborsFaaStr, NEIGHBORHOOD_NUM_maps_NEIGHBORHOOD_WP = get_acr_neighbors(candidateAcrs, ORGANISM_SUBJECT, PROTEIN_UP_DOWN)
		handle.write(neighborsFaaStr)

	from subprocess import call as execute
	execute(['rpsblast+', '-query', NEIGHBORHOOD_FAA_PATH, '-db', CDD_DB_PATH, '-evalue', '.01', '-outfmt', '7', '-out', CDD_RESULTS_PATH])

	return parse_cdd_results(CDD_RESULTS_PATH, NEIGHBORHOOD_NUM_maps_NEIGHBORHOOD_WP, MIN_NUM_PROTEINS_MATCH_CDD)




'''
	Purpose:
		Goes through the CDD results and determines which Acr/Aca regions to keep
	Arguments:
		CDD_RESULTS_PATH - Where to store psiblast results
		NEIGHBORHOOD_NUM_maps_NEIGHBORHOOD_WP - dict containing newly acquired neighbors (Protein objects) maped by the locus number
		MIN_NUM_PROTEINS_MATCH_CDD - Minimum number of proteins that had CDD hits needed to keep a locus
	Returns:
		goodNeighborhoods - list containing the locus number of candidate Acr/Aca regions that passed CDD parsing
'''
def parse_cdd_results(CDD_RESULTS_PATH, NEIGHBORHOOD_NUM_maps_NEIGHBORHOOD_WP, MIN_NUM_PROTEINS_MATCH_CDD):
	'''
		If MIN_NUM_PROTEINS_MATCH_CDD == 0 it means we shouldn't even run psiblast
	'''
	if MIN_NUM_PROTEINS_MATCH_CDD == 0:
		goodNeighborhoods = list(NEIGHBORHOOD_NUM_maps_NEIGHBORHOOD_WP.keys())
		return goodNeighborhoods

	uniqueWPHits = set()
	'''
		Traverses through psiblast results
	'''
	with open(CDD_RESULTS_PATH, 'r', 512) as handle:
		for line in handle:
			if not line.startswith('#'):
				cols = line.rstrip().split('\t')
				proteinInfo = cols[0].split('|')
				wp = proteinInfo[1]

				uniqueWPHits.add(wp)

	NEIGHBORHOOD_NUM_maps_CDD_HITS = defaultdict(lambda: 0)
	for wp in uniqueWPHits:
		for neighborhoodNum, neighborhoodWP in NEIGHBORHOOD_NUM_maps_NEIGHBORHOOD_WP.items():
			if wp in neighborhoodWP:
				NEIGHBORHOOD_NUM_maps_CDD_HITS[neighborhoodNum] += 1


	goodNeighborhoods = list()
	for neighborhoodNum in NEIGHBORHOOD_NUM_maps_CDD_HITS:
		if NEIGHBORHOOD_NUM_maps_CDD_HITS[neighborhoodNum] >= MIN_NUM_PROTEINS_MATCH_CDD:
			goodNeighborhoods.append(neighborhoodNum)

	return goodNeighborhoods



'''
	Purpose:
		Uses the indexes of good neighborhoods to get a list of Acr/Aca proteins at that index.
	Arguments:
		candidateAcrs - list of lists containing proteins in the inner list. Each protein list is an Acr/Aca locus
		goodNeighborhoods - list containing the locus number of candidate Acr/Aca regions that passed CDD parsing
	Returns:
		goodAcrs - list of Acr/Aca proteins that passed CDD filter
'''
def get_good_candidates(candidateAcrs, goodNeighborhoods):
	goodAcrs = list()
	for n in goodNeighborhoods:
		goodAcrs.append(candidateAcrs[n])

	return goodAcrs
