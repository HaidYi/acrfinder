#!/usr/bin/python3
'''
	*****************************************************************************************************
	Purpose:
		Creates mask FNA file.
		This file is a modification of the original FNA file of the organism being studied however;
			All arrays that have a spacer with desired evidence level will be 'blanked out'
			Sorounding base pairs of array will also be blanked out.
			Other sequence data is left as is including the positions.
	Author:
		Javi Gomez - https://github.com/rtomyj wrote the version 1.0.
		Haidong Yi - https://github.com/haidyi revised the codes, remove some bugs, wrote the version 2.0. 
	Project:
		This work is advised by Dr. Yanbin Yin at UNL - yyin@unl.edu
	*****************************************************************************************************
'''

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import SingleLetterAlphabet
import re

'''
	Functionality:
		Creates a list of start and end positions for all CRISPR arrays that have a spacer in the fasta spacer file
	Args:
		SELECT_SPACER_FASTA- fasta file containing spacers with desired confidence level
	Returns:
		arrayStartEndList - list, contains all start/end positions for all arrays
'''
def get_start_end_list(SELECT_SPACER_FASTA: str):
	arrayStartEndSet = set()
	with open(SELECT_SPACER_FASTA) as handle:
		for line in handle:
			if line.startswith('>'):
				spacerInfo = line.split('|')	# gets spacer array start and end
				arrayStart = spacerInfo[1].split("=")[1]
				arrayEnd = spacerInfo[2].split("=")[1]

				# reg = re.compile(r'(_)(?!.*\1)')
				# array_nc = re.sub(reg, '.', spacerInfo[0].lstrip('>'))
				array_nc = '_'.join(spacerInfo[0].lstrip('>').split('_')[0:-1])

				arrayStartEnd = tuple( [int(arrayStart), int(arrayEnd), array_nc] )
				arrayStartEndSet.add(arrayStartEnd)

	return list(arrayStartEndSet)



# '''
# 	GENERATOR

# 	Functionality:
# 		Reads an FNA file handle.
# 		Cycles through and gets the header and sequence info the sequence is describing.
# 	Args:
# 		FNA_FILE - FNA file of organism
# 	YIELDS:
# 		header and sequence info for portion of FNA file
# '''
# def fna_to_string(FNA_FILE: str):
# 	with open(FNA_FILE) as handle:
# 		header, sequence = "", ""
# 		for line in handle:
# 			if line.startswith('>'):
# 				if header != "":
# 					yield header, sequence
# 				header = line
# 				sequence = ""

# 			else:
# 				sequence += line.rstrip()

# 		if header != "" and sequence != "":
# 			yield header, sequence



# '''
# 	Functionality:
# 		Reads an FNA file.
# 		Puts all header info in a dict along with the start/end pos it is describing.
# 		Accumalates sequence info in a string.
# 	Args:
# 		FNA_FILE - FNA file of organism
# 	Returns:
# 		START_END_maps_HEADER - dict, contains the start/end position for all sequence info found in FNA as well as the header info for that specific sequence
# 		sequence - string, complete sequence read from FNA file
# '''
# def fna_to_dict(FNA_FILE: str):
# 	START_END_maps_HEADER = {}
# 	sequence = ""
# 	start = 0

# 	for head, seq in fna_to_string(FNA_FILE):
# 		end = start + len(seq)
# 		START_END_maps_HEADER[tuple([start, end])] = head
# 		sequence += seq
# 		start = end
# 	print(START_END_maps_HEADER); print('\n')
# 	print(sequence)
# 	return START_END_maps_HEADER, sequence

'''
	Functionality:
		Reads an FNA file using the SeqIO module in Biopython
		Puts all header (sequence id) in a dict along with the sequence info.
'''
def fna_to_dict(FNA_FILE: str):
	sequence_record = SeqIO.to_dict(SeqIO.parse(FNA_FILE, 'fasta'))
	return sequence_record



'''
	Functionality:

	Args:
		INTERMEDIATES - dir to store intermediate files
		FNA_FILE - FNA file of organism
		SELECT_SPACER_FASTA - fasta file containing spacer info
		MASKED_FNA - modified FNA file of organism that masks spacers
	Returns:
		MASKED_FNA - modified FNA file of organism that masks spacers
'''
def mask_fna_with_spacers(INTERMEDIATES: str, FNA_FILE: str, SELECT_SPACER_FASTA: str, MASKED_FNA:str='masked.fna'):
	from subprocess import call as execute

	MASKED_ORGANISM_DB = INTERMEDIATES + 'masked_db/'	# file to store masked FNA file
	execute(['mkdir', MASKED_ORGANISM_DB])
	MASKED_FNA = MASKED_ORGANISM_DB + MASKED_FNA


	NEIGHBORING_NUCLEOTIDES = 500	# number of BP (+/-) to also mask around the array(s)

	arrayStartEndList = get_start_end_list(SELECT_SPACER_FASTA)	# list containing the start and end of all arrays
	sequence_dict = fna_to_dict(FNA_FILE)	# header along with its start/end and the full sequence of the organism
	crisprName_maps_seqName = dict()
	for key in sequence_dict.keys():
		nc_name = key.split('.')
		if len(nc_name) > 1:
			crisprName_maps_seqName[''.join(nc_name[0:-1])] = key
		else:
			crisprName_maps_seqName[nc_name[0]] = key
	# crisprName_maps_seqName = {''.join(key.split('.')[0:-1]):key for key in sequence_dict.keys()}


	'''
		Masks CRISPR spacers from the organisms FNA file.
	'''
	for arrayStartEnd in arrayStartEndList:
		start = int(arrayStartEnd[0]) - 1
		end = int(arrayStartEnd[1])
		crispr_name = arrayStartEnd[2]  # fetch the nc # of the array
		nc_id = crisprName_maps_seqName[crispr_name]

		start, end = start - NEIGHBORING_NUCLEOTIDES, end + NEIGHBORING_NUCLEOTIDES
		sequence = str(sequence_dict[nc_id].seq) # obtain the sequence: str from the sequence dict

		if start < 0:	# corrects neg number
			start = 0
		if end > len(sequence):	# corrects if end happens to be bigger than the length of the sequence
			end = len(sequence)

		blank = 'N'
		blank = blank * ((end-start))	# multiplies masking to cover the array plus the neighboring BP's

		sequence = sequence[0:start] + blank + sequence[end:len(sequence)]	# masks sequences
		sequence_dict[nc_id].seq = Seq(sequence, SingleLetterAlphabet())

	from sys import path as sys_path
	sys_path.append('dependencies/PyGornism/')
	# from regex import string_with_limited_width

	'''
		Writes new FNA file with masking to disk.
	'''
	with open(MASKED_FNA, 'w') as handle:
		for seq_record in sequence_dict:
			SeqIO.write(sequence_dict[seq_record], handle, 'fasta')

	return MASKED_FNA, crisprName_maps_seqName


'''
The following code is used to test the functionality of masking the sequence
'''
#  mask_fna_with_spacers('INTERMEDIATES/', './GCF_003836565.1/subjects/GCF_003836565.1_ASM383656v1_genomic.fna', './GCF_003836565.1/intermediates/spacers_with_desired_evidence.fna')