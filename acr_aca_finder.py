#!/usr/bin/python3
'''
	*********************************************************************************************************
	Purpose:
		Uses three other scripts to find and then limit Acr/Aca proteins.
			Steps:
				1)	Candidate proteins are identified
				2)	Candidate proteins are compared against CDD's that imply Anti-CRISPR funtionality
				3)	Candiate proteins are compared against known areas of pathogenicity using two databases (IslandViewer/PHASTER)
				4)	Only proteins that passed 2 & 3 are put into a seperate list. These loci are considered good Acr/Aca candidates
				5)	Good candidates are written to a file
		Note: - To get more candidate Acrs, we do not recommend users use step 2) and 3)
					- Step 3) is not a default option.
	Author:
		Javi Gomez - https://github.com/rtomyj wrote the version 1.0.
		Haidong Yi - https://github.com/haidyi revised the codes, remove some bugs, wrote the version 2.0. 
	Project:
		This work is advised by Dr. Yanbin Yin at UNL - yyin@unl.edu
    **********************************************************************************************************
'''


from os import path as os_path
from os import devnull as devnull
from time import sleep
from collections import defaultdict

from find_candidate_acr_aca import first_filter, second_and_third_filter, fourth_filter, get_acr_loci, get_candidate_acr_loci, print_acrs, cleanup_acr_id_files, acr_homolog
from parse_acr_aca_with_cdd import use_cdd
from parse_acr_aca_with_db import use_gi_db_on_acr, use_pai_db_on_acr
from command_options import define_acr_finder_options, parse_acr_aca_id_options, parse_cdd_options, parse_db_options, create_sub_directories

from sys import path as sys_path
from subprocess import call as execute
sys_path.append('dependencies/PyGornism/')
from organism import Organism
from Bio import SeqIO


'''
	Summary:
		Selects only specified Acr/Aca loci. The index number of the desired loci are found in uniqueHits.
	Params:
		candidateAcrs - list of lists holding protein neighborhoods that passed filters 1, 2, 3 and 4.
		uniqueHits - set. Contains indexes of Acr/Aca loci that are going to be kept from the orignal set.
	Returns:
		finalAcrs - list of list with only the candidate Acr/Aca's that had an index number found in uniqueHits.
		CANDIDATE_INDEX_maps_FINAL_ACRS - dict that contains the contents of finalAcrs but also preserves the index of each locus from the original candidate Acrs list
'''
def finalizeLoci(candidateAcrs, uniqueHits, ORGANISM_SUBJECT, WP_ID_maps_Aca_HOMOLOG, WP_ID_maps_CDD_META, WP_ID_maps_Acr_HOMOLOG, GCF, OUTPUT_DIR):
	finalResultsFile = OUTPUT_DIR + GCF + '_guilt_by_association.out'	# file to put all selected Acr/Aca loci
	finalAcrs = []	# list containing all loci with an index number found in uniqeHits
	CANDIDATE_INDEX_maps_FINAL_ACRS = {}	# dict that maps loci index to a locus

	'''
		Traverses through all indices wanted to get corresponding locus
	'''
	for position in uniqueHits:
		finalAcrs.append(candidateAcrs[position])
		CANDIDATE_INDEX_maps_FINAL_ACRS[position] = candidateAcrs[position]

	'''
		Writes desired loci to file
	'''
	if len(finalAcrs) > 0:
		with open(finalResultsFile, 'w', 512) as handle:
			handle.write(get_acr_loci(finalAcrs, ORGANISM_SUBJECT, WP_ID_maps_Aca_HOMOLOG, WP_ID_maps_CDD_META, WP_ID_maps_Acr_HOMOLOG))
	else:
		print('\033[92m')	# color output
		print('\nNo guilt-by-association result was found, terminating...\n')
		print('\033[0m')	# makes output normal color
		exit(0)
	return finalAcrs, CANDIDATE_INDEX_maps_FINAL_ACRS, finalResultsFile





'''
	Summary:
		Uses find_candidate_acr_aca.py to find candidate Acr/Aca loci.
	Args:
		AA_THRESHOLD - amino acid threshold
		DISTANCE_THRESHOLD - max allowable distance between two adjacent proteins
		MIN_PROTEINS_IN_LOCUS - minimum allowable proteins to have per locus
		KNOWN_ACA_DATABASE - path to known aca database (.faa)
		KNOWN_ACR_DATABASE - path to known acr database (.faa)
		GFF_FILE - path to gff file of organism to analyze
		FAA_FILE - path to faa file of organism to analyze
		INTERMEDIATES - path to directory used to store intermediate files
		OUTPUT_DIR - path of directory to store all output

	Returns:
		candidateAcrs - candidate Acr/Aca loci found
		WP_ID_maps_Aca_HOMOLOG - dict. Keys are unique protein ID's that had a hit when using diamond with Aca DB. The value is the name of the Aca Hit
		ORGANISM_SUBJECT - object representing ORGANISM_SUBJECT that Acr data is for
		GCF - ID of the organism being used
'''
def identify_acr_aca(AA_THRESHOLD, DISTANCE_THRESHOLD, MIN_PROTEINS_IN_LOCUS, KNOWN_ACA_DATABASE, KNOWN_ACR_DATABASE, GFF_FILE, FAA_FILE, INTERMEDIATES, OUTPUT_DIR, isProdigalUsed):
	print('Parsing GFF and FAA file\n...')

	'''
		Creates organism object using gff and faa files
	'''
	ORGANISM_SUBJECT = Organism([GFF_FILE, FAA_FILE], isProdigalUsed, bufferSize = 30720, twoFileParse=True)    # creates Organism object used to parse gff file
	GCF = ORGANISM_SUBJECT.GCF  # obtaions GCF ID that corresponds to subject
	
	# print(GCF + '\n')
	# print(ORGANISM_SUBJECT.get_ncid_contents().items())
	# print('Done\n\n')

	'''
		Parses ORGANISM_SUBJECT file (gff) to find candidate Acr/Aca protein. An Acr locus is picked according to filters
		Filters being used are:
			1) All genes in locus are on the same strand.
			2) All encoding proteins are less than <AA_THRESHOLD> amino acids (aa) long. Default = 150.
			3) The distance between two proteins is less than <DISTANCE_THRESHOLD> long. Default = 250.
			4) At least one protein in loci is of HTH domain.
	'''


	print('Finding candidate Acr/Aca regions with the following conditions: Protein length = {0}aa, intergenic distance = {1}bb\nUsing Aca database found here -> {2}\n...'.format(AA_THRESHOLD, DISTANCE_THRESHOLD, KNOWN_ACA_DATABASE))


	'''
	Applies filter 1 and then filter 2 and 3  then finally filter 4 to the contents of the results of the first, sencond and third filters.
	obtains the candidate Acr's and the proteins that have an HTH domain
	'''
	WP_ID_maps_Aca_HOMOLOG, WP_ID_maps_Acr_HOMOLOG, candidateAcrs = fourth_filter(second_and_third_filter(first_filter(ORGANISM_SUBJECT, MIN_PROTEINS_IN_LOCUS),
		GCF, AA_THRESHOLD, DISTANCE_THRESHOLD, MIN_PROTEINS_IN_LOCUS), GCF, KNOWN_ACA_DATABASE, KNOWN_ACR_DATABASE, OUTPUT_DIR, isProdigalUsed)

	print('Filter 1-4 is Done.\n\n')

	'''
		If there are no Acr/Aca loci then there is no need to continue.
		Print message and end program
	'''
	if len(candidateAcrs) == 0:
		print('No candidate Acr/Aca regions found')
		exit(0)

	'''
	Writes acr/aca candidate file to file
	'''
	with open(INTERMEDIATES + GCF + '_candidate_acr_aca.txt', 'w', 3072) as handle:
		handle.write(get_candidate_acr_loci(candidateAcrs, ORGANISM_SUBJECT, WP_ID_maps_Aca_HOMOLOG, WP_ID_maps_Acr_HOMOLOG))
		handle.flush()

	print('\n\n')
	print_acrs(candidateAcrs, 'Candidate Acr/Aca regions:')

	return candidateAcrs, WP_ID_maps_Aca_HOMOLOG, ORGANISM_SUBJECT, GCF, WP_ID_maps_Acr_HOMOLOG





'''
	Summary:
		Uses CDD's to limit the number of Acr/Aca regions.
		Produces an FAA file that contains all the Acr/Aca regions and their neighbors +- PROTEIN_UP_DOWN.
		Uses newly produced FAA file with the CDD database. Forms a list of indices that correspond to a locus that has the CDD hit.
	Params:
		candidateAcrs - candidate Acr/Aca loci found in previous steps
		ORGANISM_SUBJECT - path to organism file for subject
		PROTEIN_UP_DOWN - number of surrounding proteins to use with CDD's
		MIN_NUM_PROTEINS_MATCH_CDD
		INTERMEDIATES - path to directory used to store intermediate files
		OUTPUT_DIR - path of directory to store all output
		GCF - ID of the organism being used
	Returns:
		neighborhoodsFromCDD - List of indices that correspond to a locus that has the CDD test.
'''
def limit_with_cdd(candidateAcrs, ORGANISM_SUBJECT, PROTEIN_UP_DOWN, MIN_NUM_PROTEINS_MATCH_CDD, INTERMEDIATES, OUTPUT_DIR, GCF, isProdigalUsed):
	'''
		Use CDD to parse
	'''
	print('Limiting Acr/Aca with CDD (psiblast+) using the following conditions: up/downstream neighbors = {0}, min number of CDD hits per locus to keep = {1}'.format(PROTEIN_UP_DOWN, MIN_NUM_PROTEINS_MATCH_CDD))

	print('...')

	NEIGHBORHOOD_FAA_PATH = INTERMEDIATES + GCF + '_candidate_acr_aca_neighborhood.faa'	# file name for temp FAA file
	CDD_DB_PATH = 'dependencies/cdds/mge'	# Path to CDD's
	CDD_RESULTS_PATH = INTERMEDIATES + GCF + '_candidate_acr_aca_cdd_results.txt'	# where to put output of psiblast+

	'''
		Abstracts CDD processing.
		Obtains indices of Acr/Aca we can keep.
	'''
	neighborhoodsFromCDD, psiblastHitFromCDD = use_cdd(candidateAcrs, ORGANISM_SUBJECT, NEIGHBORHOOD_FAA_PATH, CDD_RESULTS_PATH, CDD_DB_PATH, PROTEIN_UP_DOWN, MIN_NUM_PROTEINS_MATCH_CDD, isProdigalUsed)


	print("Number of Acr/Aca regions perserved after CDD: {0} out of {1} original loci".format(len(neighborhoodsFromCDD), len(candidateAcrs)))

	print('Keeping ' + str(sorted(neighborhoodsFromCDD)) + ' candidte Acr/Aca regions')

	print('Done\n\n')

	return neighborhoodsFromCDD, psiblastHitFromCDD





'''
	Background:
		Uses two databases that contain regions in organisms that have statistically proven areas of pathogenicty.
			- Databases are PHASTER and IslandViewer.
		The candidate Acr/Aca loci are compared with the regions in these databases.
		Since the results are dependant on whether the user wanted a lax or a strict comparison there can be very different results.
	Summary:
		Compares candidate Acr/Aca regions with regions found in the PHASTER and IslandViewer database.
		Both databases are used in parallel
		Saves the index of a locus if comparison to a region in either database is favorable.
	Params:
		candidateAcrs - candidate Acr/Aca loci found in previous steps
		USE_GI_DB - boolean, whether the GI (IslandViewer) database should be used
		USE_PAI_DB - boolean, wheter the PAI (PHASTER) database should be used
		GI_DB_FILE - path to file containing contents of GI DB
		PAI_DB_FILE - path to file containing contents of the PAI DB
		DB_LAX - boolean, whether comparison of Acr/Aca and DB entries should be lax
		DB_STRICT - boolean, whether comparison of Acr/Aca and DB entries should be lax
	Returns:
		neighborhoodsFromDB - contains all the indices of the loci that had a hit with one or both the DB's
'''
def limit_with_db(candidateAcrs, USE_GI_DB, USE_PAI_DB, GI_DB_FILE, PAI_DB_FILE, DB_LAX, DB_STRICT):
	from multiprocessing import Process, Queue

	'''
		Uses DB to parse
	'''
	queue = Queue()	# used with multiprocessing so two process share the same object (set)
	neighborhoodsFromDB = set()	# unique instances of matched Acr/Aca loci
	queue.put(neighborhoodsFromDB)

	GI_PROCESS, PAI_PROCESS = None, None
	if USE_GI_DB:	# ISLANDVIEWER
		print('Using IslandViewer to limit Acr/Aca regions using the following conditions: ', end='')

		if DB_STRICT:
			print('Strict (all proteins in locus must have same IslandViewer hit)')
		else:
			print('Lax (at least one protein in locus musth have a IslandViewer hit)')

		GI_PROCESS = Process(target = use_gi_db_on_acr, args = (candidateAcrs, DB_STRICT, DB_LAX, GI_DB_FILE, queue))
		GI_PROCESS.start()


	if USE_PAI_DB:	# PHASTER
		print('Using PHASTER to limit Acr/Aca regions using the following conditions: ', end='')

		if DB_STRICT:
			print('Strict (all proteins in locus must have same PHASTER hit)')
		else:
			print('Lax (at least one protein in locus musth have a PHASTER hit)')

		PAI_PROCESS = Process(target = use_pai_db_on_acr, args = (candidateAcrs, DB_STRICT, DB_LAX, PAI_DB_FILE, queue))
		PAI_PROCESS.start()


	print('...')

	'''
		Joins processes
	'''
	if GI_PROCESS != None:
		GI_PROCESS.join()
	if PAI_PROCESS != None:
		PAI_PROCESS.join()


	neighborhoodsFromDB = queue.get()	# obtains the indexes of the loci that had hits (we want to keep)
	print('Number of Acr/Aca regions perserved after IslandViewer/PHASTER: {0} out of {1} original loci'.format(len(neighborhoodsFromDB), len(candidateAcrs)))
	print('Keeping ' + str(sorted(neighborhoodsFromDB)) + ' candidate Acr/Aca regions')

	print('Done\n\n')

	return neighborhoodsFromDB





'''
	Summary:
		Gets all the command line arguments user has used/wanted from executing this script (whether directly or indirectly)
	Params:
		parser - OptionParser that contians defined user options
	returns:
		User args
'''
def get_options(parser = None, fna_faaNeeded=True):

	'''
		Defines parser with the arguments needed to run acr_aca_finder.
	'''
	if parser == None:
		from optparse import OptionParser
		parser = OptionParser()
		define_acr_finder_options(parser)

	options = parser.parse_args()[0]

	AA_THRESHOLD, DISTANCE_THRESHOLD, MIN_PROTEINS_IN_LOCUS, KNOWN_ACA_DATABASE, KNOWN_ACR_DATABASE, OUTPUT_DIR, GFF_FILE, FAA_FILE = parse_acr_aca_id_options(options, fna_faaNeeded=fna_faaNeeded)

	PROTEIN_UP_DOWN, MIN_NUM_PROTEINS_MATCH_CDD = parse_cdd_options(options)

	GI_DB_FILE, PAI_DB_FILE, USE_GI_DB, USE_PAI_DB, DB_STRICT, DB_LAX = parse_db_options(options)


	'''
		Tries to find a GCF ID for naming purposes
	'''
	from sys import path as sys_path
	sys_path.append('dependencies/PyGornism/')
	from regex import Regex

	match = Regex.GCF_REGEX.search(FAA_FILE)

	if match != None:
		GCF = match.group(1)
	else:
		GCF = 'Undefined_Organism'

	'''
		Creates a temp directory if no dir is specified
	'''
	if OUTPUT_DIR == "./":
		from datetime import datetime as date
		from re import compile
		from subprocess import call as execute

		DATE_TIME_REGEX = compile(r'([0-9]+)-([0-9]+)-([0-9]+) ([0-9]+):([0-9]+):([0-9]+)\.([0-9]+)')
		match = DATE_TIME_REGEX.search(str(date.now()))
		OUTPUT_DIR = '_'.join([GCF, match.group(1), match.group(2), match.group(3), match.group(4), match.group(
			5), match.group(7)]) + '/'  # dir path to store all output and depenencies created or used by script
		print('CRISPRCasFinder and Acr/Aca Finder results will be found here -> ' + OUTPUT_DIR)

		execute(['mkdir', OUTPUT_DIR])


	INTERMEDIATES, SUBJECTS_DIR = create_sub_directories(OUTPUT_DIR)  # creates directory to store intermediate files and dir to store organism files that are going to be used through execution


	return AA_THRESHOLD, DISTANCE_THRESHOLD, MIN_PROTEINS_IN_LOCUS, KNOWN_ACA_DATABASE, KNOWN_ACR_DATABASE, OUTPUT_DIR, GFF_FILE, FAA_FILE, PROTEIN_UP_DOWN, MIN_NUM_PROTEINS_MATCH_CDD, GI_DB_FILE, PAI_DB_FILE, USE_GI_DB, USE_PAI_DB, DB_STRICT, DB_LAX, INTERMEDIATES, SUBJECTS_DIR, GCF





'''
	Summary:
		Executes all methods etc needed for the proper execution of acr_aca_finder.py (this script)
	Params:
		User Args
	Returns:
		Nothing
'''
def acr_aca_run(AA_THRESHOLD, DISTANCE_THRESHOLD, MIN_PROTEINS_IN_LOCUS, KNOWN_ACA_DATABASE, KNOWN_ACR_DATABASE, OUTPUT_DIR, GFF_FILE, FAA_FILE, PROTEIN_UP_DOWN, MIN_NUM_PROTEINS_MATCH_CDD, GI_DB_FILE, PAI_DB_FILE, USE_GI_DB, USE_PAI_DB, DB_STRICT, DB_LAX, INTERMEDIATES, SUBJECTS_DIR, GCF, isProdigalUsed):
	'''
		Both a valid GFF and FAA file are needed.
		This will make sure they both exist before commencing.
	'''
	if GFF_FILE == "" or FAA_FILE == "":
		print('Both a GFF and an FAA file need to be specified')
		print('Use --help for usage')
		exit(-1)

		'''
			Checks gff and faa file to see if they exist.
			If any of them don't exist, exit with -1
		'''
		if not os_path.isfile(GFF_FILE):
			print('gff file doesn\'t exist')
			exit(-1)
		if not os_path.isfile(FAA_FILE):
			print('faa file doesn\'t exist')
			exit(-1)

	# KEEP_INTERMEDIATES = True

	'''
	Search Acrs with homolog methods: diamond (diamond against the whole faa)
	'''
	DIAMOND_DATA_BASE = INTERMEDIATES + GCF + '_acr_diamond_database'
	DIAMOND_ACR_QUERY = KNOWN_ACR_DATABASE
	DIAMOND_ACRHOMOLOG_FILE = INTERMEDIATES + GCF + '_acr_homolog_result.txt'
	with open(devnull, 'w') as DEV_NULL:
		execute(['diamond', 'makedb', '--in', FAA_FILE, '-d', DIAMOND_DATA_BASE], stdout=DEV_NULL)
	with open(devnull, 'w') as DEV_NULL:
		execute(['diamond', 'blastp', '-q', DIAMOND_ACR_QUERY, '--db', DIAMOND_DATA_BASE, '-e', '.01', '-f', '6', 'qseqid', 'sseqid', 'pident', 'slen', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', '-o', DIAMOND_ACRHOMOLOG_FILE], stdout=DEV_NULL)

	#acr_hit_record, _ = acr_homolog(FAA_FILE, DIAMOND_ACRHOMOLOG_FILE, INTERMEDIATES, GCF, isProdigalUsed)

	candidateAcrs, WP_ID_maps_Aca_HOMOLOG, ORGANISM_SUBJECT, GCF, WP_ID_maps_Acr_HOMOLOG = identify_acr_aca(AA_THRESHOLD, DISTANCE_THRESHOLD, MIN_PROTEINS_IN_LOCUS, KNOWN_ACA_DATABASE, KNOWN_ACR_DATABASE, GFF_FILE, FAA_FILE, INTERMEDIATES, OUTPUT_DIR, isProdigalUsed)
	neighborhoodsFromCDD, WP_ID_maps_CDD_META = limit_with_cdd(candidateAcrs, ORGANISM_SUBJECT, PROTEIN_UP_DOWN, MIN_NUM_PROTEINS_MATCH_CDD, INTERMEDIATES, OUTPUT_DIR, GCF, isProdigalUsed)


	if USE_GI_DB or USE_PAI_DB:
		neighborhoodsFromDB = limit_with_db(candidateAcrs, USE_GI_DB, USE_PAI_DB, GI_DB_FILE, PAI_DB_FILE, DB_LAX, DB_STRICT)
	else:
		neighborhoodsFromDB = set()
		print('IslandViewer/PHASTER option not selected. \nSkipping...\n\n')


	uniqueHits = set(neighborhoodsFromCDD).union(neighborhoodsFromDB)
	print('Number of Acr/Aca regions perserved after combining CDD/PHASTER/IslandViewer results: {0} out of {1} original locus'.format(len(uniqueHits), len(candidateAcrs)))
	print('Keeping ' + str(sorted(uniqueHits)) + ' candidate Acr/Aca regions')

	uniqueHits, CANDIDATE_INDEX_maps_FINAL_ACRS, finalResultsFile = finalizeLoci(candidateAcrs, uniqueHits, ORGANISM_SUBJECT, WP_ID_maps_Aca_HOMOLOG, WP_ID_maps_CDD_META, WP_ID_maps_Acr_HOMOLOG, GCF, OUTPUT_DIR)

	#finalHomologFile = finalizeHomolog(acr_hit_record, uniqueHits, OUTPUT_DIR, GCF, isProdigalUsed)
	#print('\033[31m')  # highlight the homology based method
	#print('\nAcr Homolog final results can be found here -> {0}'.format(os_path.abspath(finalHomologFile)))
	#print('\033[0m\n\n')
	print_acrs(CANDIDATE_INDEX_maps_FINAL_ACRS, 'Final Acr/Aca regions:', isDict=True)

	if OUTPUT_DIR == "":
		OUTPUT_DIR = '.'
	print('Results can be found here -> {0}'.format(os_path.abspath(OUTPUT_DIR)))

	return finalResultsFile

