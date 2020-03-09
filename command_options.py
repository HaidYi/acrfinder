#!/usr/bin/python3
'''
	*****************************************************************************************************
	Purpose:
		Contains methods used to define and parse command line arguments for Acr/Aca Finder.

	Author:
		Javi Gomez - https://github.com/rtomyj wrote the version 1.0.
		Haidong Yi - https://github.com/haidyi revised the codes, remove some bugs, wrote the version 2.0. 
	Project:
		This work is advised by Dr. Yanbin Yin at UNL - yyin@unl.edu
	*****************************************************************************************************
'''
from os import path as os_path
from subprocess import call as call



'''
	Functionality:
		Validates paths for directories and files that user submits.
		If the path is a dir and it doesn't exist then it will be created.
		If path isn't a dir:
			If path is empty or if path is to a file that doesn't exist, message is displayed and program terminates
	Args:
		PATH - absolute/relative path to resource
		ARG - short description of argument being checked, appended to err message
		IS_DIR - if PATH is for a dir
		EXT - whether program should stop exection if paths validity is bad
	Returns:
		invalid - whether PATH is valid or not
'''
def validate_path(PATH: str, ARG:str, IS_DIR:bool=False, EXT:bool=True):
	invalid = False
	'''
		Checks to see if dir exists if not it is created.
	'''
	if IS_DIR and not os_path.isdir(PATH):
		call(['mkdir', PATH])

	if not IS_DIR:
		if PATH == "":
			print(ARG + ' path is an empty string but required.')
			invalid = True
		elif not os_path.isfile(PATH):
			print('Path is invalid/DNE @ ' + PATH + ' for dependency ' + ARG)
			invalid = True

	if invalid and EXT:
		exit(-1)
	return invalid



'''
	Functionality:
		Defines options for IO.
	Args:
		parser - will contain all user options
	Returns:
		Nothing
'''
def define_io_options(parser, INPUT_DIR_HELP:str= 'directory containing data to be parsed by script', OUTPUT_DIR_HELP:str= 'directory used to store output'):
	parser.add_option('-o', '--outDir', action = 'store', dest = 'outDir', help = OUTPUT_DIR_HELP, default = './')



'''
	Functionality:
		Defines additional options available when running acr_aca_cri_runner.py
		Mainly it will allow a user to submit an FNA file for processing.
	Args:
		parser - will contain all user options
	Returns:
		Nothing
'''
def define_master_options(parser):
	parser.add_option('-n', '--inFNA', action = 'store', dest = 'fna', help = 'input fna file', default = "")
	parser.add_option('-z', '--genomeType', action = 'store', dest = 'genome', help = 'Virus = V, Archaea = A, Bacteria = B', default = "V")



'''
	Functionality:
		Defines all acceptable options a user can use to modify how Acr/Aca proteins are identified.
	Args:
		parser - will contain all user options
	Returns:
		Nothing
'''
def define_acr_aca_id_options(parser):
	parser.add_option('-m', '--aaThresh', action = 'store', dest = 'aaThresh', help = 'max amino acid length of protein. Default = 200', default = '200')
	parser.add_option('-d', '--distThresh', action = 'store', dest = 'distThresh', help = 'max intergenic distance between proteins. Default = 150', default = '150')
	parser.add_option('-r', '--minProteins', action = 'store', dest = 'minProteins', help = 'Minimum number of proteins a locus must have in order to keep as candidate. Default = 2', default = '2')
	parser.add_option('-t', '--aca', action = 'store', dest = 'acaDB', help = 'Known Aca file (.faa) to diamond candidate aca in candidate Acr-Aca loci', default = 'dependencies/diamond_query/401-aca.faa')
	parser.add_option('-u', '--acr', action='store', dest='acrDB', help='Known Acr file (.faa) to diamond the homolog of Acr', default = 'dependencies/diamond_query/known-acr.faa')
	
	parser.add_option('--escape_file', action='store', dest='escapeDB', help='', default = 'dependencies/escape_list_HTH_MGE')
	parser.add_option('--cdd_db', action='store', dest='cddDB', help='', default = 'dependencies/cdd/Cdd')
	parser.add_option('--blsType', action='store', dest='blastType', choices=['blastp', 'rpsblast'], default='blastp', help='which blast mode to use')

	parser.add_option('--identity', action= 'store', type=str, dest='Identity', default='30', help='diamond aca identity')
	parser.add_option('--coverage', action= 'store', type=str, dest='Coverage', default='0.8', help='diamond aca coverage')
	parser.add_option('--e_value', action= 'store', type=str, dest='E_Value', default='0.01', help='diamond aca e_value')

	parser.add_option('--blast_slack', action='store', dest='blsSlack', default=5000, type=int, help='how far an Acr/Aca locus is allowed to be from a blastn hit to be considered high confidence')
	parser.add_option('--no_dmd_ss', action='store_true', dest='nodmdSS', default=False, help='whether to use sensive mode of diamond.')

	parser.add_option('-f', '--inGFF', action = 'store', dest = 'gff', help = 'input gff file', default = '')
	parser.add_option('-a', '--inFAA', action = 'store', dest = 'faa', help = 'input faa file', default = '')



'''
	Functionality:
		Defines arguments that modify how CDD's are used to identify potential Acr/Aca regions.
	Args:
		parser - will contain all user options
	Returns:
		Nothing
'''
def define_cdd_options(parser):
	parser.add_option('-e', '--proteinUpDown', action = 'store', dest = 'proteinUpDown', help = 'How many proteins upstream and downstream to gather to use with CDD\'s. Default = 5', default = '10')
	parser.add_option('-c', '--minCDDProteins', action = 'store', dest = 'minCDDProteins', help = 'Minimum number of proteins that should have a CDD match in order to include Acr/Aca locus. Default = 2', default = '1')



'''
	Functionality:
		Defines options that modify how PHASTER/IslandViewr databases are used to identify potential Acr/Aca regions.
	Args:
		parser - will contain all user options
	Returns:
		Nothing
'''
def define_db_options(parser):
	parser.add_option('-g', '--gi', action = 'store', dest = 'gi', help = 'Use GI DB to parse subject, t=true, f=false. Default = f', default = 'f')
	parser.add_option('-p', '--pai', action = 'store', dest = 'pai', help = 'Use PAI DB to parse subject, t=true, f=false. Default = f', default = 'f')
	parser.add_option('-s', '--strict', action = 'store', dest = 'strict', help = 'Only print out regions that are strictly within known GI/PAI, t=true, f=false. Default = f', default = 'f')
	parser.add_option('-l', '--lax', action = 'store', dest = 'lax', help = 'Prints regions that touch a GI/PAI, t=true, f=false. Default = f', default = 't')



'''
	Functionality:
		Defines additional options available when running acr_aca_cri_runner.py
		Mainly it will allow a user to submit an FNA file for processing.
	Args:
		parser - will contain all user options
	Returns:
		Nothing
'''
def define_crispr_cas_finder_options(parser):
	parser.add_option('-y', '--arrayEvidence', action='store', dest='arrayEvidence', help='The min evidence level of crisprs to fastify (min = 1, max = 4). Default = 3', default='3')



'''
	Functionality:
		Calls all methods that define options that are needed/available for Acr/Aca identification and the execution of acr_aca_finder.py
	Args:
		parser - will contain all user options
	Returns:
		Nothing
'''
def define_acr_finder_options(parser):
	define_acr_aca_id_options(parser)
	define_cdd_options(parser)
	define_db_options(parser)
	define_io_options(parser)



'''
	Functionality:
		Calls all methods that define options that are needed/available for Acr/Aca identification as well as all options needed to run CRISPRCasFinder (options needed to run acr_aca_cri_runner.py)
	Args:
		parser - will contain all user options
	Returns:
		Nothing
'''
def define_acr_aca_cri_runner_options(parser):
	define_master_options(parser)
	define_acr_finder_options(parser)
	define_crispr_cas_finder_options(parser)



'''
	Functionality:
		Processes arguments that contian IO options user wanted.
		Creates Output directory.
	Args:
		options - object containing all options available with the values the user wanted for each option.
	Returns:
		OUTPUT_DIR - path to the output directory
'''
def parse_io_options(options):
	'''
		Corrects directory name to include terminating '/'
		Creates a new output dir if the one specified by user doesn't exist
	'''
	OUTPUT_DIR = options.outDir
	if not OUTPUT_DIR.endswith('/'):
		OUTPUT_DIR += '/'
		validate_path(OUTPUT_DIR, 'Output Dir', IS_DIR=True)

	return OUTPUT_DIR



'''
	Functionality:
		Processes arguments specific to the script acr_aca_cri_runner.py
	Args:
		options - object containing all options available with the values the user wanted for each option.
	Returns:
		FNA - path to FNA file user wants parsed
'''
def parse_master_options(options):
	validate_path(options.fna, 'FNA file')
	genome = options.genome[0].upper()
	return options.fna, genome



''''
	Functionality:
		Processes all argumnets needed to identify Acr/Aca protiens and run acr_aca_finder.py
		Establishes paths to an FAA and GFF file.
		Establishes path to HTH HMM DB.
		Sets parameters for amino acids length, distance between proteins, etc needed for identifying proteins.
	Args:
		options - object containing all options available with the values the user wanted for each option.
	Returns:
		AA_THRESHOLD - int, max number of amino acids (how big) an Acr/Aca protein is allowed to be
		DISTANCE_THRESHOLD - int, max distnace between two adjacent proteins
		MIN_PROTEINS_IN_LOCUS - int, min number of proteins in an Acr/Aca locus
		KNOWN_ACA_DATABASE - string, path to known aca database (.faa)
		KNOWN_ACR_DATABASE - string, path to known acr database (.faa)
		OUTPUT_DIR - string, path to a directory to store all output
		GFF_FILE - string, filename of GFF to parse
		FAA_FILE - string, filename of FAA to parse
'''
def parse_acr_aca_id_options(options, fna_faaNeeded=True):
	GFF_FILE, FAA_FILE = options.gff, options.faa	# path to gff and faa files
	validate_path(GFF_FILE, 'GFF file', EXT=fna_faaNeeded), validate_path(FAA_FILE, 'FAA file', EXT=fna_faaNeeded)	# validates GFF and FAA

	KNOWN_ACA_DATABASE = options.acaDB
	KNOWN_ACR_DATABASE = options.acrDB
	validate_path(KNOWN_ACA_DATABASE, 'KNOWN ACA DB')	# validates users custom ACA DB
	validate_path(KNOWN_ACR_DATABASE, 'KNOWN ACR DB') # validates users custom ACR DB

	'''
		Changes threshold for amino acids, intergenic distance and minimum number of proteins for a locus to user specified values.
	'''
	if not options.aaThresh.isdigit():
		print('Amino acid threshold (-m/--aaThresh) must be an integer >= 0')
		exit(-1)
	else:
		AA_THRESHOLD = int(options.aaThresh)

	if not options.distThresh.isdigit():
		print('Intergenic distance threshold (-d/--distThresh) must be an integer >= 0')
		exit(-1)
	else:
		DISTANCE_THRESHOLD = int(options.distThresh)

	if not options.minProteins.isdigit():
		print('Intergenic distance threshold (-r/--minProteins) must be an integer >= 0')
		exit(-1)
	else:
		MIN_PROTEINS_IN_LOCUS = int(options.minProteins)


	OUTPUT_DIR = parse_io_options(options)	# gets the output dir

	BLAST_SLACK = options.blsSlack # blast slack paramter
	NO_DIAMOND_SS = options.nodmdSS # flag whether use sensitive mode
	BLAST_TYPE = options.blastType # which blast comp to use, 'blastp' or 'rpsblast+'

	ESCAPE_DBFILE = options.escapeDB # escape file path
	CDD_DBFILE = options.cddDB # cdd database file path

	# diamond parameters
	IDENTITY = options.Identity  # the identity parameters of diamond aca
	COVERAGE = options.Coverage # the coverage parameter of diamond aca
	E_VALUE = options.E_Value # the e-value parameter of diamond aca

	return AA_THRESHOLD, DISTANCE_THRESHOLD, MIN_PROTEINS_IN_LOCUS, KNOWN_ACA_DATABASE, KNOWN_ACR_DATABASE, OUTPUT_DIR, BLAST_SLACK, NO_DIAMOND_SS, ESCAPE_DBFILE, CDD_DBFILE, BLAST_TYPE, IDENTITY, COVERAGE, E_VALUE, GFF_FILE, FAA_FILE



'''
	Functionality:
		Parses the options that pertain to Acr/Aca discrimination using CDD
	Args:
		options - object containing all options available with the values the user wanted for each option.
	Returns:
		PROTEIN_UP_DOWN - int, number of proteins needed to make/create a protein neighborhood used for parsing
		MIN_NUM_PROTEINS_MATCH_CDD - int, min number of proteins that must have a CDD hit from neighborhood
'''
def parse_cdd_options(options):
	'''
		Checks user options for number of proteins to use (upstream and downstream) for neighborhoods adn the minimun number of proteins that must have a CDD hit in order to keep locus.
		Makes sure both numbers are integers. Terminates script with error if it's not.
	'''
	if options.proteinUpDown != "":
		if not options.proteinUpDown.isdigit():
			print('CDD neighborhood (-e/--proteinUpDown) must be an integer >= 0')
			exit(-1)
		else:
			PROTEIN_UP_DOWN = int(options.proteinUpDown)

	if options.minCDDProteins != "":
		if not options.minCDDProteins.isdigit():
			print('Minumum number of proteins with a CDD hit (-c/--minCDDProteins) must be an integer >= 0')
			exit(-1)
		else:
			MIN_NUM_PROTEINS_MATCH_CDD = int(options.minCDDProteins)

	return PROTEIN_UP_DOWN, MIN_NUM_PROTEINS_MATCH_CDD



'''
	Functionality:
		Parses only the options that pertain to Acr/Aca discrimination using PHASTER/IslandViewer
	Args:
		options - object containing all options available with the values the user wanted for each option.
	Returns:
		GI_DB_FILE - string, path to GI DB file
		PAI_DB_FILE - string, path to PAI DB file
		USE_GI_DB - boolean, if GI DB should be used
		USE_PAI_DB - boolean, if PAI DB should be used
		DB_STRICT - boolean, if strict search should be used with GI/PAI DB
		DB_LAX - boolean, if lax search should be used with GI/PAI DB
'''
def parse_db_options(options):
	'''
		DB DISCRIMINATION PARAMS
	'''
	GI_DB_FILE, PAI_DB_FILE = 'dependencies/all_gis_islandviewer_iv4.txt', 'dependencies/z_DNA_fragment_DB.header'  # database files for GI and PAI

	gi, pai, strict, lax = options.gi[0].upper(), options.pai[0].upper(), options.strict[0].upper(), options.lax[0].upper()
	if strict == 'T': lax = 'F'	# lax is the default, if user wants strict mode then it overrides the default lax

	print('Using GI DB:\t\t{0}\nUsing PAI DB:\t\t{1}\nUsing Strict Mode:\t{2}\nUsing Lax Mode:\t\t{3}\n'.format(gi, pai, strict, lax))


	'''
		Looks to see if PHASTER/IslandViewer DB should be used
	'''
	if gi == 'T':	USE_GI_DB = True
	else:			USE_GI_DB = False

	if pai == 'T':	USE_PAI_DB = True
	else:			USE_PAI_DB = False


	'''
		Looks to see if strict or lax algorithm should be used
	'''
	if strict == 'T':
		DB_STRICT = True
		DB_LAX = False
	elif lax == 'T':
		DB_LAX = True
		DB_STRICT = False


	'''
		If user doesn't want either strict or lax output and error message will display and the program will exit
	'''
	if DB_STRICT == False and DB_LAX == False:
		print('Both strict and lax output options cannot be false')
		print('Use --help for details')
		exit(-1)


	return GI_DB_FILE, PAI_DB_FILE, USE_GI_DB, USE_PAI_DB, DB_STRICT, DB_LAX



'''
	Functionality:
		Parses user options that pertain to CRISPRCasFinder.
	Args:
		options - object containing all options available with the values the user wanted for each option.
	Returns:
		EVIDENCE_LEVEL - int range->[1,4] the evidence level desired of each crispr array.
'''
def parse_crispr_cas_finder_options(options):
	if options.arrayEvidence != "" and options.arrayEvidence.isdigit():
		EVIDENCE_LEVEL = int(options.arrayEvidence)
		if EVIDENCE_LEVEL > 4 or EVIDENCE_LEVEL < 1:
			print('Invalid evidence level value.\nUse -h for usage.')
			exit(-1)

	else:
		print('Invalid evidence level value.\nUse -h for usage.')
		exit(-1)

	return EVIDENCE_LEVEL



'''
	Functionality:
		creates directories to store files produced/being parsed
	Args:
		OUTPUT_DIR - directory to hold all output
	Returns:
		INTERMEDIATES - string, path to directory that stores intermediate files
		SUBJECTS_DIR - string, path to directory that stores organism/sequence data that will be parsed
'''
def create_sub_directories(OUTPUT_DIR:str):
	INTERMEDIATES = OUTPUT_DIR + 'intermediates/'
	# dir path to store all genome data user submitted with query or data produced using prodigal
	SUBJECTS_DIR = OUTPUT_DIR + 'subjects/'

	if not os_path.isdir(INTERMEDIATES):
		call(['mkdir', INTERMEDIATES])

	if not os_path.isdir(SUBJECTS_DIR):
		call(['mkdir', SUBJECTS_DIR])

	return INTERMEDIATES, SUBJECTS_DIR
