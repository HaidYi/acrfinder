#!/usr/bin/python3
'''
	*****************************************************************************************************
	Purpose:
		Handles execution of CRISPRCasFinder.
		Parses CRISPRCasFinder output and creates a fasta file with CRISPR spacers found.
		Creates a file that blocks (masks) all spacers from the parsed organism.
		uses blastn to find occurences of all CRISPR spacers within the organism using the masked organism file.

	Author:
		Javi Gomez - https://github.com/rtomyj wrote the version 1.0.
		Haidong Yi - https://github.com/haidyi revised the codes, remove some bugs, wrote the version 2.0. 
	Project:
		This work is advised by Dr. Yanbin Yin at UNL - yyin@unl.edu
	*****************************************************************************************************
'''
from subprocess import call as execute

'''
	Functionality:
		Uses makeblastdb to create a database using the masked FNA file.
		Uses blastn along with the database created and fasta file of spacers to find other occurences of the spacers within the same organism.
	Args:
		SELECT_SPACER_FASTA - fasta file containing all CRISPR spacers found on organism
		MASKED_FNA - file of organism that has masked regions where CRISPR spacers were found
		DIR - directory to store files
	Returns:
		BLAST_FILE - string, path to file containing blast results
'''
def blast_spacers(SELECT_SPACER_FASTA: str, MASKED_FNA: str, DIR: str):
	BLAST_FILE = DIR + 'blast_out.txt'	# file to store blastn results
	REQURIED_IDENTITY = '95.0'	# percent identity needed in order to consider blast result

	print('\n\nUsing blastn on CRISPR arrays...')
	execute(['makeblastdb', '-in', MASKED_FNA, '-dbtype', 'nucl'])	# makes DB for blast using the masked FNA file created in previous steps

	with open(BLAST_FILE, 'w') as handle:
		execute(['blastn', '-db', MASKED_FNA, '-query', SELECT_SPACER_FASTA, '-outfmt', '6', '-perc_identity', REQURIED_IDENTITY], stdout=handle)	# runs blastn

	print('\nDone\n\n')

	return BLAST_FILE

'''
	Functionality:
		Runs CRISPRCasFinder.
		Creates a fasta file containing all CRISPR spacers found.
		Creates a masked fna file that takes the orginal FNA file and puts non sensical but valid characters in the locations where spacers were found.
			This will prevent spacers from matching themselves but will preserve the length of the sequence and location of all other sequences.
		Runs blastn.
	Args:
		CRISPR_CAS_FINDER_EXECUTABLE - file used to execute CRISPRCasFinder program
		CRISPR_CAS_FINDER_SO - file needed by CRISPRCasFinder
		FNA_FILE - path to FNA file of organism to study
		OUTPUT_DIR - directory to store any output
		INTERMEDIATES - directory to store intermediate files
		EVIDENCE_LEVEL - evidence level desired of the spacers
	Returns:
		BLAST_FILE - string, path to file containing blast results
'''

def crispr_cas_runner(CRISPR_CAS_FINDER_EXECUTABLE: str, CRISPR_CAS_FINDER_SO: str, FNA_FILE: str, OUTPUT_DIR: str, INTERMEDIATES: str, GENOME_TYPE: str, EVIDENCE_LEVEL:int=3):
	from mask_fna_with_spacers import mask_fna_with_spacers as mask_fna
	from fastafy_select_spacers import fastafy_select_spacers as fastafy

	CRISPR_CAS_OUTPUT = OUTPUT_DIR + 'CRISPRCas_OUTPUT'	# directory to store CRISPRCasFinder output

	'''
		Executes a slightly different variant of the CRISPRCasFinder command depending if genome is Archaea or Bacteria
	'''
	if GENOME_TYPE == 'A':   # 'A' means Archaea
		execute(['perl', CRISPR_CAS_FINDER_EXECUTABLE, '-in', FNA_FILE, '-cas', '-so', CRISPR_CAS_FINDER_SO, '-out', CRISPR_CAS_OUTPUT, '-q', '-ArchaCas'])	# runs CRISPRCasFinder
	else:
		execute(['perl', CRISPR_CAS_FINDER_EXECUTABLE, '-in', FNA_FILE, '-cas', '-so', CRISPR_CAS_FINDER_SO, '-out', CRISPR_CAS_OUTPUT, '-q'])	# runs CRISPRCasFinder

	print('\n\n\nCRISPRCasFinder has finished executing...\n\n')

	SELECT_SPACER_FASTA = fastafy(CRISPR_CAS_OUTPUT, INTERMEDIATES, EVIDENCE_LEVEL)	# creates fasta of spacers
	if SELECT_SPACER_FASTA == None:
		return None, None	# if there were no spacers found then running blastn is impossible
	MASKED_FNA, CRISPR_NAME_maps_SEQ_NAME = mask_fna(INTERMEDIATES, FNA_FILE, SELECT_SPACER_FASTA)	# creates masked fna file using original FNA and spacers

	return blast_spacers(SELECT_SPACER_FASTA, MASKED_FNA, INTERMEDIATES), CRISPR_NAME_maps_SEQ_NAME	# returns blast file

