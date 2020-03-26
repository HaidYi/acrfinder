#!/usr/bin/python3
'''
	*****************************************************************************************************
	Purpose:
		Executes CRISPRCasFinder and acr_aca_finder.py to find Cas/CRISPR data as well as Acr/Aca proteins.
		If CRISPRCasFinder finds no CRISPR-Cas info then acr_aca_finder.py will not execute.
		If Acr/Aca proteins are found, CRISPRCasFinder is used to classify loci ranging from low, med, and high confidence
		based on whether spacer self-target occurs and the distance between target position and Anti-Crispr protein loci.
	Author:
		Javi Gomez - https://github.com/rtomyj wrote the version 1.0.
		Haidong Yi - https://github.com/haidyi revised the codes, remove some bugs, wrote the version 2.0. 
	Project:
		This work is advised by Dr. Yanbin Yin at UNL - yyin@unl.edu
	*****************************************************************************************************
'''
import subprocess
from subprocess import call as execute
from os import devnull as devnull
from os import path as os_path
from optparse import OptionParser
from collections import defaultdict

from acr_aca import getLocus, getLocusStartAndEnd
from acr_aca_finder import get_options as aaf_get_options, acr_aca_run as aa_runner
from command_options import define_acr_aca_cri_runner_options, parse_master_options, parse_crispr_cas_finder_options, parse_io_options, create_sub_directories
from find_candidate_acr_aca import acr_homolog

'''
	Functionality:
		Uses prodigal to generate a GFF and FAA file.
	Args:
		FNA_FILE - FNA file of organism, needed by prodigal
		DIR - directory used to store output
		GCF - the GCF ID of the organism
	Returns:
		GFF_FILE - string, path to GFF file created by prodigal
		FAA_FILE - string, path to FAA file created by prodigal
'''
def create_gff_faa_with_prodigal(FNA_FILE, DIR, GCF):
	print('Using prodigal to generate gff and faa file...')
	GFF_FILE = DIR + GCF + '.gff'
	FAA_FILE = DIR + GCF + '.faa'

	execute(['prodigal', '-f', 'gff', '-i', FNA_FILE, '-o', GFF_FILE, '-q', '-a', FAA_FILE ])
	execute(['cp', FNA_FILE, SUBJECTS_DIR])

	print('Done\n\n')


	return GFF_FILE, FAA_FILE


def generate_filter_file(acr_aca_file, rpsdb_file, INTERMEDIATES):
    sequence = ''
    with open(acr_aca_file, 'r') as acr_file:
        for locus in getLocus(acr_file):
            proteins = locus.split('\n')

            if proteins[0].startswith('#'):
                proteins = proteins[1:len(proteins)]

            for protein in proteins:
                item_list = protein.split('\t')
                if item_list[-4] == 'Acr':
                    sequence += '>'+item_list[6]+'\n'+item_list[-1]+'\n'
    
    fasta_file = INTERMEDIATES + 'final_acr.fasta'
    with open(fasta_file, 'w') as handle:
        handle.write(sequence)
    # run rpsblast here
    rpsblast_file = INTERMEDIATES + 'final_rpsblast.out'
    execute(['rpsblast+', '-query', fasta_file, '-db', rpsdb_file, '-evalue', '.001', '-outfmt', '6', '-out', rpsblast_file])

    return fasta_file, rpsblast_file

'''
	Functionality:
		Looks at the result file 
'''
def cdd_loci_filter(proteins, escape_set, fasta_file, rpsblast_file):
	#make rps result set
	with open(rpsblast_file, "r") as rps:
		rps_dic={}
		for line in rps.readlines():
			line = line.rstrip().split()
			wp = line[0]
			cddNUM = line[1].replace("CDD:", "")
			if wp not in rps_dic.keys():
				rps_dic.setdefault(wp, [cddNUM])
			if wp in rps_dic.keys():
				rps_dic[wp].append(cddNUM)
	rps_result = set()
	for key in rps_dic.keys():
		if any(v for v in rps_dic[key] if v in escape_set) is False:
			rps_result.add(key)

	# filter the protein using the rpsblast results
	proteinID_list = []
	for protein in proteins:
		item_list = protein.split('\t')
		if item_list[-4] == 'Acr':
			proteinID_list.append(item_list[6])
	if any(v for v in proteinID_list if v in rps_result) is False:
		return False
	else:
		return True



'''
	Functionality:
		Looks at blast result file to create a list of start/end positions (with slack) and the results info (such as array start/end, spacer start/end and CRISPR Cas system).
		Reads Acr/Aca results file to get the start and end of all loci.
		Compares the start and end of blast results (with slack) and Acr/Aca loci to see if any Acr/Aca locus should be considered High confidence.
		Inserts three new columns into the existing Acr/Aca file;
			classification - low, med, hight
			blast info for blast hits (with slack) that are NEAR Acr/Aca loci.
			blast info for blast hits (with slack) that are NOT NEAR Acr/Aca loci.
	Args:
		BLAST_FILE - file containing blast results between masked.fna and generated FASTA file of spacers
		ACR_ACA_FILE - file containing Acr/Aca results, will be modified to include Acr/Aca classification
		OUTPUT_DIR - directory used to store output files
	Returns:
		Nothing
'''

def classify_acr_aca(BLAST_FILE, CRISPR_NAME_maps_SEQ_NAME, ACR_ACA_FILE, OUTPUT_DIR, INTERMEDIATES, BLAST_SLACK, ESCAPE_DBFILE, CDD_DBFILE, GENOME_TYPE):
	BLAST_HIT_SLACK = BLAST_SLACK	# how far an Acr/Aca locus is allowed to be from a blastn hit to be considered high confidence
	lookForHigh = False	# whether there is a possibility of a high confidence classifcation

	'''
		Defines the default confidence level for the genome.
		Low confidence is the default classification if there were no self targeting blast hits using the masked FNA file and the spacer FASTA file.
		Medium confidence is the default classification if 
			- there was at least one self targeting blast hit.
			- Anti-CRISPR accession and target accession are different.
			- Anti-CRISPR and target locate in the same accession but their distance is more than 5000 bp.
		High confidence is the default classification if
			- Anti-CRISPR and target locate in teh same accession and their distance is less than 5000 bp.
	'''
	if BLAST_FILE == None or os_path.getsize(BLAST_FILE) == 0:
		defaultClassification = 'Low Confidence'	# no self targeted hits
	else:
		defaultClassification = 'Medium Confidence'	# at least one self targeted hit

		START_END_maps_BLAST_HIT = defaultdict(list)	# dict to hold all blast hit info for given region
		'''
			Reads blastn output and builds a dict of regions (start/end) and the corresponding blast info
		'''
		with open(BLAST_FILE) as blastHandle:
			for line in blastHandle:
				cols = line.rstrip().split('\t')
				target_accession, hitStart, hitEnd = cols[1], cols[8], cols[9]
				spacers_info = cols[0].split('|')
				locus_accession = CRISPR_NAME_maps_SEQ_NAME[ '_'.join(spacers_info[0].split('_')[0:-1]) ]  # fetch the locus accession of spacers.
				
				'''
					Get the self-target info 
				'''
				spacer_nc = 'Spacer Accession={0}'.format(locus_accession)
				spacer_pos = 'Spacer_Pos={0}-{1}'.format(spacers_info[4].split('=')[1], spacers_info[5].split('=')[1])
				cas_type = spacers_info[-1]
				
				target_nc = 'Target Accession={0}'.format(target_accession)
				target_pos = 'Target_Pos={0}-{1}'.format(hitStart, hitEnd)

				START_END_maps_BLAST_HIT[(hitStart, hitEnd, target_accession)].append('|'.join([spacer_nc, spacer_pos, cas_type, target_nc, target_pos]))

				lookForHigh = True

	fasta_file, rpsblast_file = generate_filter_file(ACR_ACA_FILE, CDD_DBFILE, INTERMEDIATES)
	# obtain the escape list
	escape_set = set()
	for line in open(ESCAPE_DBFILE).readlines():
		line=line.rstrip()
		escape_set.add(line)

	with open(ACR_ACA_FILE) as acrFile:
		newLoci = []
		header = ''

		'''
			Cycles through all Acr/Aca loci that were generated using acr_aca_finder.py
		'''
		for locus in getLocus(acrFile):
			proteins = locus.split('\n')	# gets all protein info from a locus
			highConfidence = False	# whether locus is high confidence
			closeBlastHits = ""	# contains blast hits that are close to Acr/Aca region
			farBlastHits = ""	# contains blast hits that aren't close to Acr/Aca region

			'''
				Gets the header and removes it from proteins list.
				There will only be one header to remove.
			'''
			if proteins[0].startswith('#'):
				header = proteins[0]
				proteins = proteins[1:len(proteins)]


			locusStart, locusEnd = getLocusStartAndEnd(proteins)	# start and ending position of current Acr/Aca locus

			# add cdd filter there
			is_cddfilter = (GENOME_TYPE != 'V')
			if is_cddfilter and cdd_loci_filter(proteins, escape_set, fasta_file, rpsblast_file):
				continue


			if lookForHigh:  # if self targeting blast hit
				for startEnd in START_END_maps_BLAST_HIT.keys():
					# blastHitStart, blastHitEnd, target_accession = startEnd.split('-')
					blastHitStart, blastHitEnd, target_accession = startEnd
					blastHitStart, blastHitEnd = int(blastHitStart), int(blastHitEnd)	# start/end position +/- slack of current blast hit info
					blastHitStart, blastHitEnd = blastHitEnd - BLAST_HIT_SLACK, blastHitStart + BLAST_HIT_SLACK

					if blastHitStart < 0:  	# corrects over subtraction
						blastHitStart = 0
					'''
						If Acr/Aca locus is within the blast hit (with slack) then the blast hit info is added to closeBlastHits.
						Otherwise it is added to farBlastHits.
					'''
					if proteins[0].split('\t')[2] == target_accession: # if target accession equals acr accession
						if max(locusStart, blastHitStart) < min(locusEnd, blastHitEnd):
							highConfidence = True

							for infoList in START_END_maps_BLAST_HIT[startEnd]:
								if closeBlastHits != "":
									closeBlastHits += ','  # split different self-target with ','
								closeBlastHits += infoList

						else:
							for infoList in START_END_maps_BLAST_HIT[startEnd]:
								if farBlastHits != "":
									farBlastHits += ','  # split different self-target with ','
								farBlastHits += infoList
					else:
						for infoList in START_END_maps_BLAST_HIT[startEnd]:
							if farBlastHits != "":
								farBlastHits += ','
							farBlastHits += infoList

			if farBlastHits == "":
				farBlastHits = '---'

			newProteins = []
			for protein in proteins:
				if highConfidence:
					newProteins.append('\t'.join(['High Confidence', protein, closeBlastHits, farBlastHits]))
				else:
					newProteins.append('\t'.join([defaultClassification, protein, '---', farBlastHits]))

			newLoci.append('\n'.join(newProteins))

	with open(ACR_ACA_FILE, 'w') as handle:
		header = '#Classification\t' + header.lstrip('#') + '\tSelf Target w/in {0} BP\tSelf Target Outside {0} BP\n\n'.format(BLAST_HIT_SLACK)
		handle.write(header)
		handle.write('\n\n'.join(newLoci))

	'''
		Functionality: 
	'''
def acr_aca_homolog(FAA_FILE, GFF_FILE, FNA_FILE, MIN_PROTEINS_IN_LOCUS, AA_THRESHOLD, DISTANCE_THRESHOLD, KNOWN_ACR_DATABASE, INTERMEDIATES, GCF, OUTPUT_DIR, isProdigalUsed):
	DIAMOND_DATA_BASE = INTERMEDIATES + GCF + '_acr_diamond_database'
	DIAMOND_ACR_QUERY = KNOWN_ACR_DATABASE
	DIAMOND_ACRHOMOLOG_FILE = INTERMEDIATES + GCF + '_acr_homolog_result.txt'
	with open(devnull, 'w') as DEV_NULL:
		execute(['diamond', 'makedb', '--in', FAA_FILE, '-d', DIAMOND_DATA_BASE], stdout=DEV_NULL, stderr=subprocess.STDOUT)
	with open(devnull, 'w') as DEV_NULL:
		execute(['diamond', 'blastp', '-q', DIAMOND_ACR_QUERY, '--db', DIAMOND_DATA_BASE, '-e', '.01', '-f', '6', 'qseqid', 'sseqid', 'pident', 'slen', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', '-o', DIAMOND_ACRHOMOLOG_FILE], stdout=DEV_NULL, stderr=subprocess.STDOUT)

	_, HOMOLOG_FINAL_RESULT_FILE = acr_homolog(FAA_FILE, GFF_FILE, FNA_FILE, MIN_PROTEINS_IN_LOCUS, AA_THRESHOLD, DISTANCE_THRESHOLD, DIAMOND_ACRHOMOLOG_FILE, OUTPUT_DIR, isProdigalUsed)

	print('\033[31m')  # highlight the homology based method with red color
	if HOMOLOG_FINAL_RESULT_FILE != None:
		print('The result of using homolog-based methods can be found -> {0}\n'.format(os_path.abspath(HOMOLOG_FINAL_RESULT_FILE)))
	else:
		print("No homology-based result was found. Try to use guilt-by-association method in the next step.")
	print('\033[0m')


'''
	********************************************************************************************************
	**************************************	MAIN FUNCTIONALITY	********************************************
	*******************************************************************************************************
'''
GENOME_TYPE_maps_NAME = {'A': 'Archaea', 'B': 'Bacteria', 'V': 'Virus'}

'''
	Defines options needed for script execution
'''
parser = OptionParser()
define_acr_aca_cri_runner_options(parser)

'''
	Gets user options
'''
options = parser.parse_args()[0]
FNA_FILE, GENOME_TYPE = parse_master_options(options)
EVIDENCE_LEVEL = parse_crispr_cas_finder_options(options)
AA_THRESHOLD, DISTANCE_THRESHOLD, MIN_PROTEINS_IN_LOCUS, KNOWN_ACA_DATABASE, KNOWN_ACR_DATABASE, OUTPUT_DIR, BLAST_SLACK, NO_DIAMOND_SS, ESCAPE_DBFILE, CDD_DBFILE, BLAST_TYPE, IDENTITY, COVERAGE, E_VALUE, THREADS_NUM, GFF_FILE, FAA_FILE, PROTEIN_UP_DOWN, MIN_NUM_PROTEINS_MATCH_CDD, GI_DB_FILE, PAI_DB_FILE, USE_GI_DB, USE_PAI_DB, DB_STRICT, DB_LAX, INTERMEDIATES, SUBJECTS_DIR, GCF = aaf_get_options(parser, fna_faaNeeded=False)
isProdigalUsed = False  # Flag whether prodigal is used to generate .gff and .faa files.
CRISPR_CAS_FINDER_EXECUTABLE, CRISPR_CAS_FINDER_SO, NUM_CPUS = 'dependencies/CRISPRCasFinder/CRISPRCasFinder.pl', 'dependencies/CRISPRCasFinder/sel392v2.so', '4'  # CRISPRCasFinder files

print('Treating organism as {0}\n\n\n'.format(GENOME_TYPE_maps_NAME[GENOME_TYPE]))


'''
	Validating fna/faa input files
'''
if GFF_FILE == "" or FAA_FILE == "":  # creates new gff and faa file if user didn't specify either one
	GFF_FILE, FAA_FILE = create_gff_faa_with_prodigal(FNA_FILE, SUBJECTS_DIR, GCF)
	isProdigalUsed = True
else:
	execute(['cp', FNA_FILE, GFF_FILE, FAA_FILE, SUBJECTS_DIR])	# copies organism files to subjects dir
	gffFileParts, faaFileParts = GFF_FILE.split('/'), FAA_FILE.split('/')
	GFF_FILE, FAA_FILE = SUBJECTS_DIR + gffFileParts[len(gffFileParts) - 1], SUBJECTS_DIR + faaFileParts[len(faaFileParts) - 1]	# gets path to GFF and FAA file in subjects dir
	fnaFileParts = FNA_FILE.split('/')
	FNA_FILE = SUBJECTS_DIR + fnaFileParts[len(gffFileParts) - 1]	# gets path to FNA file in subjects dir



from crispr_cas_runner import crispr_cas_runner as cc_runner
'''
	Runs CRISPRCasFinder.
	If there are no CRISPR arrays found then the program will only use homology-based method.
	If there are CRISPR arrays the program contiunes and will find arrays with an evidence level equal to or greater than EVIDENCE_LEVEL.
	BLAST_FILE will be null/None when there are no arrays with the wanted evidence level.
'''
acr_aca_homolog(FAA_FILE, GFF_FILE, FNA_FILE, MIN_PROTEINS_IN_LOCUS, AA_THRESHOLD, DISTANCE_THRESHOLD, KNOWN_ACR_DATABASE, INTERMEDIATES, GCF, OUTPUT_DIR, isProdigalUsed)  # No crispr-cas system was found above the given evidence level, terminating...
if GENOME_TYPE != 'V':
	BLAST_FILE, CRISPR_NAME_maps_SEQ_NAME = cc_runner(CRISPR_CAS_FINDER_EXECUTABLE, CRISPR_CAS_FINDER_SO, FNA_FILE, OUTPUT_DIR, INTERMEDIATES, GENOME_TYPE, EVIDENCE_LEVEL) # Get crispr spacers and spacers blast hit file
	if BLAST_FILE == None and CRISPR_NAME_maps_SEQ_NAME == None:
		print('Acr-Aca finder program was finished using homolog methods, terminating...\n')
		exit(0)
	else:
		ACR_ACA_FILE = aa_runner(AA_THRESHOLD, DISTANCE_THRESHOLD, MIN_PROTEINS_IN_LOCUS, KNOWN_ACA_DATABASE, KNOWN_ACR_DATABASE, OUTPUT_DIR, NO_DIAMOND_SS, BLAST_TYPE, IDENTITY, COVERAGE, E_VALUE, THREADS_NUM, GFF_FILE, FAA_FILE, FNA_FILE, PROTEIN_UP_DOWN, MIN_NUM_PROTEINS_MATCH_CDD, GI_DB_FILE, PAI_DB_FILE, USE_GI_DB, USE_PAI_DB, DB_STRICT, DB_LAX, INTERMEDIATES, SUBJECTS_DIR, GCF, isProdigalUsed)
		classify_acr_aca(BLAST_FILE, CRISPR_NAME_maps_SEQ_NAME, ACR_ACA_FILE, OUTPUT_DIR, INTERMEDIATES, BLAST_SLACK, ESCAPE_DBFILE, CDD_DBFILE, GENOME_TYPE)	# classifies Acr/Aca proteins
else: # if organism is virus, no need to run CRISPRCas-Finder
	ACR_ACA_FILE = aa_runner(AA_THRESHOLD, DISTANCE_THRESHOLD, MIN_PROTEINS_IN_LOCUS, KNOWN_ACA_DATABASE, KNOWN_ACR_DATABASE, OUTPUT_DIR, NO_DIAMOND_SS, BLAST_TYPE, IDENTITY, COVERAGE, E_VALUE, THREADS_NUM, GFF_FILE, FAA_FILE, FNA_FILE, PROTEIN_UP_DOWN, MIN_NUM_PROTEINS_MATCH_CDD, GI_DB_FILE, PAI_DB_FILE, USE_GI_DB, USE_PAI_DB, DB_STRICT, DB_LAX, INTERMEDIATES, SUBJECTS_DIR, GCF, isProdigalUsed)

	classify_acr_aca(None, None, ACR_ACA_FILE, OUTPUT_DIR, INTERMEDIATES, BLAST_SLACK, ESCAPE_DBFILE, CDD_DBFILE, GENOME_TYPE)	# classifies Acr/Aca proteins

