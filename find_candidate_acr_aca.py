#!/usr/bin/python3
'''
	*********************************************************************************************************

	Purpose:
		Parses ORGANISM_SUBJECT file (gff) to find candidate Acr/Aca protein. An Acr locus is picked according to filters
		Filters being used are:
			1) All genes in locus are on the same strand.
			2) All encoding proteins are less than <AA_THRESHOLD> amino acids (aa) long. Default = 150.
			3) The distance between two proteins is less than <DISTANCE_THRESHOLD> long. Default = 250.
			4) At least one protein in loci is of HTH domain.

		Additionally, A locus must have the following attributes:
			1) Must contain a minimum number of proteins. Default = 2
			2) All proteins must be adjacent, that is no other proteins must exist between the two consecutive proteins in a locus.
	Author:
		Javi Gomez - https://github.com/rtomyj wrote the version 1.0.
		Haidong Yi - https://github.com/haidyi revised the codes, remove some bugs, wrote the version 2.0. 
	Project:
		This work is advised by Dr. Yanbin Yin at UNL - yyin@unl.edu
    **********************************************************************************************************
'''
from sys import path as sys_path
from os import devnull as devnull
from subprocess import call as execute
from collections import defaultdict

sys_path.append('dependencies/PyGornism/')
from organism import Organism
from formated_output import write_faa
from regex import Regex
from Bio import SeqIO
import math
import heapq



'''
	Functionality:
		Applies first filter to ORGANISM_SUBJECT (extracts neighborhoods of genes where each neighborhood is on the same strand).
		Creates a list of these neighborhoods.
		Adds neighborhoods to another list, candidateAcrs

		Eg:
			ProteinA = +, ProteinB = +, ProteinC = -, ProteinD = +
			neigborhoods = [ProteinA, ProteinB], [ProteinC], [ProteinD]
	Note:
		A neighborhood must have at least two proteins to be added to the candidateAcrs list
	Args:
		ORGANISM_SUBJECT - organism object (from organism.py) that contains data such as protein and ncid information for the organism being parsed.
	Returns:
		candidateAcrs - list of lists holding like strand neighborhoods of Protein objects
'''
def first_filter(ORGANISM_SUBJECT, MIN_PROTEINS_IN_LOCUS = 2):
	candidateAcrs = list()	# holds nighborhoods

	'''
		Loops through all ncids in the ORGANISM_SUBJECT and forms neighborhoods where each item in neigbhor hood is on same strand.
	'''
	for ncid, proteinList in ORGANISM_SUBJECT.get_ncid_contents().items():
		strand = ""
		strandList = list()	# empty neighborhood
		for protein in proteinList:
			'''
				Since the current neighborhood strand doesn't match the next protein a new neighborhood list needs to be created.
				Before that a check is done to make sure neighborhood is holds at least two elements,
					if so the protein list gets added to candidateAcrs
			'''
			if strand != protein.strand:
				strand = protein.strand
				if len(strandList) >= MIN_PROTEINS_IN_LOCUS:
					candidateAcrs.append(strandList)

				strandList = list()

			'''
				Appends protein to newly created list or previously populated list
			'''
			strandList.append(protein)


		'''
			Adds last remaining strand list after loop terminates if it has 2 or more proteins.
		'''
		if len(strandList) >= MIN_PROTEINS_IN_LOCUS:
			candidateAcrs.append(strandList)

	return candidateAcrs



'''
	Functionality:
		Applies second (protein amino acid threshold) and third (intergenic region bp threshold) filters to candidateAcrs proteins.
		Creates a new list of neighborhoods where each neighborhood is also a list of proteins that passed filter 2 and 3
	Args:
		candidateAcrs - list of lists holding like strand neighborhoods of Protein objects with first filter applied
		AA_THRESHOLD - Max length a protein is allowed to be to remain a candidate Acr/Aca protein.
		DISTANCE_THRESHOLD - Max allowable distance between two adjacent proteins.
	Returns:
		candidateAcrs_2 - list of lists holding like strand neighborhoods of Protein objects with the second and third filter applied to it
'''
def second_and_third_filter(candidateAcrs, OUTPUT_DIR, GCF, AA_THRESHOLD = 150, DISTANCE_THRESHOLD = 250, MIN_PROTEINS_IN_LOCUS = 2):
	candidateAcrs_2 = list()
	# INTERMEDIATES = OUTPUT_DIR + 'intermediates/'
	# flattenCandidates = [acr for loci in candidateAcrs for acr in loci]
	# CANDIDATES_FAA_FILE = INTERMEDIATES + GCF + '_candidate_acr_aca_filter_1.faa'
	# write_faa(flattenCandidates, CANDIDATES_FAA_FILE)
	'''
		Loops through neighborhoods of proteins, each neighborhood a list
	'''
	for proteinList in candidateAcrs:
		aaList = list()	# contains neighborhoods of proteins with amino acid length less than AA_THRESHOLD

		'''
			Loops through individual Protein objects within a specific neighborhood
		'''
		for protein in proteinList:
			aaLength = int( (protein.end - protein.start + 1) / 3 )	# amino acid length of protein

			if aaLength <= AA_THRESHOLD:
				'''
					Will attempt to add protein to aaList
					Three things can happen:
						1) aaLength isn't less than AA_THRESHOLD -> do nothing
						2) aaLength is less than AA_THRESHOLD and (protein position of current protein is exactly one more than protein previouls parsed and distance between the two proteins is less than DISTANCE_THRESHOLD) -> add current protein to aaList
						3) aaLength is less than AA_THRESHOLD and not (protein position of current protein is exactly one more than protein previouls parsed and distance between the two proteins is less than DISTANCE_THRESHOLD) -> check current aaList - if len >=2 add to candidateAcrs_2, create new aaList, add current protein to new aaList.
				'''
				if len(aaList) == 0:
					aaList.append(protein)
				else:
					'''
						Checks to see position of current protein is right after previous protein and checks filter 3
					'''
					if aaList[len(aaList) -1].position == (protein.position - 1) and (protein.start - aaList[len(aaList) -1].end <= DISTANCE_THRESHOLD):
						aaList.append(protein)
					else:
						if len(aaList) >= MIN_PROTEINS_IN_LOCUS:
							candidateAcrs_2.append(aaList)
						aaList=list()
						aaList.append(protein)

		'''
			Adds last aaList to candidateAcrs_2 if condition holds
		'''
		if len(aaList) >= MIN_PROTEINS_IN_LOCUS:
			candidateAcrs_2.append(aaList)

	return candidateAcrs_2



'''
	Functionality:
		Applies fourth filter (only locus that contain at least one protein with an HTH domain are kept)
		Creates a new list of neighborhoods where each neighborhood is also a list of proteins that passed filter 4
	Args:
		candidateAcrs - list of lists holding protein neighborhoods that passed filters 1, 2, and 3.
		GCF - ID of organism used as a prefix when naming/creating files.
		KNOWN_ACA_DATABASE - path to known aca database (.faa)
		KNOWN_ACR_DATABASE -path to known acr database (.faa)
		OUTPUT_DIR - Directory used to store results/created files.
	Returns:
		WP_ID_maps_HTH_DOMAIN - dict. Keys are unique protein ID's that had a hit when using hmmscan with HTH DB. The value is the name of the HTH domain.
		candidateAcrs_filter4 - list of lists holding like strand neighborhoods of Protein objects with the second and third filter applied to it
'''
def fourth_filter(candidateAcrs, GCF, KNOWN_ACA_DATABASE, KNOWN_ACR_DATABASE, OUTPUT_DIR, isProdigalUsed):
	'''
		Creates an faa file of the remaining proteins of organism
	'''
	INTERMEDIATES = OUTPUT_DIR + 'intermediates/'	# directory to store intermediate (not primarily important) files

	'''
		If there are no Acr/Aca loci then there is no need to continue.
		Print message and end program
	'''
	if len(candidateAcrs) == 0:
		print('No candidate Acr/Aca regions found after 1,2,3 filter.')
		exit(0)

	flattenCandidates = [acr for loci in candidateAcrs for acr in loci]	# converts list of lists into a 1D list
	CANDIDATES_FAA_FILE = INTERMEDIATES + GCF + '_candidate_acr_aca.faa'
	write_faa(flattenCandidates, CANDIDATES_FAA_FILE)

	'''
	  Uses diamond on newly created faa file to find proteins that has homelog in the query: 401 Acas.
		Uses user uploaded .faa to make the database
		Sends stdout to /dev/null.
	'''
	DIAMOND_DATA_BASE = INTERMEDIATES + GCF + '_candidate_acr_aca_diamond_database'
	DIAMOND_QUERY = KNOWN_ACA_DATABASE
	DIAMOND_OUTPUT_FILE = INTERMEDIATES + GCF + '_candidate_acr_aca_diamond_result.txt'

	DIAMOND_ACR_QUERY = KNOWN_ACR_DATABASE
	DIAMOND_ACRHOMOLOG_FILE = INTERMEDIATES + GCF + '_candidate_acr_homolog_result.txt'

	with open(devnull, 'w') as DEV_NULL:
		execute(['diamond', 'makedb', '--in', CANDIDATES_FAA_FILE, '-d', DIAMOND_DATA_BASE], stdout=DEV_NULL)
	
	with open(devnull, 'w') as DEV_NULL:
		execute(['diamond', 'blastp', '-q', DIAMOND_ACR_QUERY, '--db', DIAMOND_DATA_BASE, '-e', '.01', '-f', '6', 
		         'qseqid', 'sseqid', 'pident', 'slen', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', '-o', DIAMOND_ACRHOMOLOG_FILE], stdout=DEV_NULL)

	with open(devnull, 'w') as DEV_NULL:
		execute(['diamond', 'blastp', '-q', DIAMOND_QUERY, '--db', DIAMOND_DATA_BASE, '-e', '.01', '--id', '40', '--query-cover', '0.8', '-f', '6', 
		         'qseqid', 'sseqid', 'pident', 'slen', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', '-o', DIAMOND_OUTPUT_FILE], stdout=DEV_NULL)

	'''
	    Parses DIAMOND_OUTPUT_FILE created by diamond blastp
			Populates dict of dicts with useful diamond output
	'''
	WP_ID_maps_Aca_HOMOLOG = defaultdict(list)	# dict of list, holds info of all wp's with Aca Hit (Considering the difference using prodigal.)
	WP_ID_maps_Acr_HOMOLOG = defaultdict(list)  # dict of list, wp --> Acr Homologs (Considering the difference using prodigal.)
	
	with open(DIAMOND_OUTPUT_FILE, 'r', 512) as handle:
		for line in handle:
			cols = line.rstrip().split('\t')
			if int(cols[3]) > 50 and int(cols[3]) < 200:
				aca_hit, regionInfo = cols[0], cols[1].split('|')
				start, end = regionInfo[2], regionInfo[3]
				evalue = float(cols[11])

				if isProdigalUsed:
					wp = '-'.join(regionInfo[0:2])
				else:
					wp = regionInfo[1]

				WP_ID_maps_Aca_HOMOLOG[wp].append( {'aca_hit': aca_hit, 'start': start, 'end': end, 'evalue': evalue} )
	
	with open(DIAMOND_ACRHOMOLOG_FILE, 'r', 512) as handle:
		for line in handle:
			cols = line.rstrip().split('\t')
			if int(cols[3]) < 200:
				acr, regionInfo, pident = cols[0], cols[1].split('|'), cols[2]
				evalue = float(cols[11])
				
				if isProdigalUsed:
					wp = '-'.join(regionInfo[0:2])
				else:
					wp = regionInfo[1]
				
				WP_ID_maps_Acr_HOMOLOG[wp].append( {'acr_hit': '|'.join([acr, pident]), 'evalue': evalue} )
	
	# Fetch the most significant output from the results of diamond
	for wp_id in WP_ID_maps_Aca_HOMOLOG.keys():
		WP_ID_maps_Aca_HOMOLOG[wp_id] = heapq.nsmallest(1, WP_ID_maps_Aca_HOMOLOG[wp_id], key=lambda s: s['evalue'])[0]
	
	for wp_id in WP_ID_maps_Acr_HOMOLOG.keys():
		WP_ID_maps_Acr_HOMOLOG[wp_id] = heapq.nsmallest(1, WP_ID_maps_Acr_HOMOLOG[wp_id], key=lambda s: s['evalue'])[0]
		

	'''
		Creates new candidate list that contains only loci where at least one protein is of aca hit.
	'''
	candidateAcrs_filter4 = list()
	for locus in candidateAcrs:
		add = False
		for protein in locus:
			'''
				Attempts to find protein ID (wp) in the dict holding wp's that were found to be of aca hit.
			'''
			if protein.id in WP_ID_maps_Aca_HOMOLOG.keys():
				add = True
				break

		if add:
			candidateAcrs_filter4.append(locus)

	return WP_ID_maps_Aca_HOMOLOG, WP_ID_maps_Acr_HOMOLOG, candidateAcrs_filter4



'''
	Functionality:
		Creates a string representation of all Acr loci found.
		Each locus is seperated by a line break (empty line).
		All columns are tab delimited
	Args:
		candidateAcrs - list of lists holding nieghborhood of Protein objects that passed all filters
		ORGANISM_SUBJECT - object representing ORGANISM_SUBJECT that Acr data is for.
		header - optional, what to print before the actual data is displayed.

	Returns:
		output - string representing all Acr loci.
'''
def get_acr_loci(candidateAcrs, ORGANISM_SUBJECT, WP_ID_maps_Aca_HOMOLOG, WP_ID_maps_CDD_META, WP_ID_maps_Acr_HOMOLOG, header='#GCF\tPosition\tNC ID\tStart\tEnd\tStrand\tProtein ID\taa Length\tAcr/Aca\tCDD MetaData\tAcr_Hit|pident\n'):
	output = header
	gcf = ORGANISM_SUBJECT.GCF

	'''
		Formatted output.
		Prints 'Acr' if protein has no Aca hit, otherwise it prints out the name of the domain
	'''
	for locus in candidateAcrs:
		for protein in locus:
			'''
				Searches WP_ID_maps_Aca_HOMOLOG. Only wp/proteins that had an Aca hit will be keys in the dict: WP_ID_maps_Aca_HOMOLOG
				Searches WP_ID_maps_Acr_HOMOLOG. Only wp/proteins that had an Acr hit will be keys in the dict: WP_ID_maps_Acr_HOMOLOG
			'''
			if protein.id in WP_ID_maps_Acr_HOMOLOG.keys():
				acr_hit = WP_ID_maps_Acr_HOMOLOG[protein.id]['acr_hit']
			else:
				acr_hit = '-'
			
			if protein.id in WP_ID_maps_CDD_META.keys():
				cddMetaData = WP_ID_maps_CDD_META[protein.id]['sid'] + '|' + WP_ID_maps_CDD_META[protein.id]['evalue']
				# cddMetaData = WP_ID_maps_CDD_META[protein.wp]['qid'] + '|' + WP_ID_maps_CDD_META[protein.wp]['evalue'] + '-' + WP_ID_maps_CDD_META[protein.wp]['sid']
			else:
				cddMetaData = '-'

			if protein.id in WP_ID_maps_Aca_HOMOLOG.keys():
				aca_hit = WP_ID_maps_Aca_HOMOLOG[protein.id]['aca_hit']
				output += '\t'.join([gcf, str(protein.position), protein.nc, str(protein.start), str(protein.end), protein.strand, protein.wp, str( int(((protein.end - protein.start + 1) / 3)) ), aca_hit, cddMetaData, acr_hit])
			else:
				output += '\t'.join([gcf, str(protein.position), protein.nc, str(protein.start), str(protein.end), protein.strand, protein.wp, str( int(((protein.end - protein.start + 1) / 3)) ), 'Acr', cddMetaData, acr_hit])


			output += '\n'

		output += '\n'

	return output


'''
	Functionality:
		Creates a string representation of all Acr loci found.
		Each locus is seperated by a line break (empty line).
		All columns are tab delimited
	Args:
		candidateAcrs - list of lists holding nieghborhood of Protein objects that passed all filters
		ORGANISM_SUBJECT - object representing ORGANISM_SUBJECT that Acr data is for.
		header - optional, what to print before the actual data is displayed.

	Returns:
		output - string representing all Acr loci (without cdd meta).
'''
def get_candidate_acr_loci(candidateAcrs, ORGANISM_SUBJECT, WP_ID_maps_Aca_HOMOLOG, WP_ID_maps_Acr_HOMOLOG, header='#GCF\tPosition\tNC ID\tStart\tEnd\tStrand\tProtein ID\taa Length\tDomain\tAcr_Hit|pident\n'):
	output = header
	gcf = ORGANISM_SUBJECT.GCF

	'''
		Formatted output.
		Prints 'Acr' if protein has no Aca hit, otherwise it prints out the name of the domain
	'''
	for locus in candidateAcrs:
		for protein in locus:
			'''
				Searches WP_ID_maps_Aca_HOMOLOG. Only wp/proteins that had an Aca hit will be keys in the dict: WP_ID_maps_Aca_HOMOLOG
				Searches WP_ID_maps_Acr_HOMOLOG. Only wp/proteins that had an Acr hit will be keys in the dict: WP_ID_maps_Acr_HOMOLOG
			'''
			if protein.id in WP_ID_maps_Acr_HOMOLOG.keys():
				acr_hit = WP_ID_maps_Acr_HOMOLOG[protein.id]['acr_hit']
			else:
				acr_hit = '-'
			if protein.id in WP_ID_maps_Aca_HOMOLOG.keys():
				aca_hit = WP_ID_maps_Aca_HOMOLOG[protein.id]['aca_hit']
				output += '\t'.join([gcf, str(protein.position), protein.nc, str(protein.start), str(protein.end), protein.strand, protein.wp, str( int(((protein.end - protein.start + 1) / 3)) ), aca_hit, acr_hit])
			else:
				output += '\t'.join([gcf, str(protein.position), protein.nc, str(protein.start), str(protein.end), protein.strand, protein.wp, str( int(((protein.end - protein.start + 1) / 3)) ), 'Acr', acr_hit])

			output += '\n'

		output += '\n'

	return output



'''
	Functionality:
		Prints a short summary of the Acr/Aca loci contained in list.
		Can also print Acr/Aca loci and preserve the original position of all locus from the original Acr/Aca candidate list
			Use isDict=True
			Pass a dict as the first argument where the dicts keys are the original position and the values are the list of proteins in locus.
	Args:
		acrs - list of lists containing acr/aca regions. Each sublist contains protein objects
	Returns:
		Nothing
'''
def print_acrs(acrs, header, isDict=False):
	print(header)
	print('\033[92m')	# color output
	if not isDict:
		for position, locus in enumerate(acrs):
			for protein in locus:
				print("[{4}] {0}-pos({5})\t{1}\t{2}\t{3}".format(protein.wp, protein.nc, protein.start, protein.end, position, protein.position))
			print("")
	else:	# prints preserving position form candidate proteins
		for position, locus in acrs.items():
			for protein in locus:
				print("[{4}] {0}-pos({5})\t{1}\t{2}\t{3}".format(protein.wp, protein.nc, protein.start, protein.end, position, protein.position))
			print("")

	print('\033[0m')	# makes output normal color
	print('\n\n')



'''
	Functionality:
		Removes files created while trying to find Acr/Aca proteins
		Files being removed are:
			faa file generated for the candidate Acr/Aca proteins
			HTH HMM output file produced by hmmscan
			parsed HTH HMM file produced by hmmparser.sh
	Args:
		GCF - The prefix given to all output files
		OUTPUT_DIR - the directory where output files are stored
	Returns:
		Nothing
'''
def cleanup_acr_id_files(GCF, OUTPUT_DIR):
	CANDIDATES_FAA_FILE = OUTPUT_DIR + GCF + '_candidate_acr.faa'
	CANDIDATES_HTH_FILE = OUTPUT_DIR + 'candidate_hth.txt'  # temp file to hold all hth hits from hmmscan
	PARSED_HTH_FILE = OUTPUT_DIR + 'parsed_hth.txt' # temp file that contains hth hits that have been condensed using hmmparser


	execute(['rm', CANDIDATES_FAA_FILE, CANDIDATES_HTH_FILE, PARSED_HTH_FILE])

'''
	Functionality:
		Extract the acr hit sequence from the diamond output file and FAA_FILE
	Args:
		FAA_FILE - The FAA file user uploaded or input
		DIAMOND_ACRHOMOLOG_FILE - the output file of diamond (ground_truth acr against FAA_FILE)
		INTERMEDIATES - path to directory used to store intermediate files
		GCF - The prefix given to all output files
		isProdigalUsed - flag the .gff and .faa files generated by prodigal
	Returns:
		acr_hit_record - the homolog sequence dict (protein id -> SeqRecord data (defined in Biopython https://biopython-cn.readthedocs.io/zh_CN/latest/cn/chr03.html ) )
		ACR_INTERMEDIATES_FASTA_FILE - the homolog sequence data (KNOWN-ACR_DATABASE against FAA_FILE)
'''
def acr_homolog(FAA_FILE, DIAMOND_ACRHOMOLOG_FILE, INTERMEDIATES, GCF, isProdigalUsed):
	acr_hit_record = dict()
	record_dict = SeqIO.to_dict(SeqIO.parse(FAA_FILE, 'fasta'))
	ACR_INTERMEDIATES_FASTA_FILE = INTERMEDIATES + GCF + '_acr_homolog_result.fasta'

	with open(DIAMOND_ACRHOMOLOG_FILE, 'r', 512) as handle:
		for line in handle:
			cols = line.rstrip().split('\t')
			acr, wp, pident, slen, evalue = cols[0], cols[1], cols[2], cols[3], float(cols[10])
			if isProdigalUsed:
				protein_info_list = record_dict[wp].description.split('#')
				protein_start = protein_info_list[1].strip()
				protein_end = protein_info_list[2].strip()
				protein = 'Protein({0}-{1})'.format(protein_start, protein_end)
				nc_id = wp[0:wp.rfind('_')]

				if nc_id + protein in acr_hit_record:
					if evalue < acr_hit_record[nc_id + protein]['evalue']:
						acr_hit_record[nc_id + protein]['record'].id = '|'.join([nc_id, protein, slen, pident])
						acr_hit_record[nc_id + protein]['record'].description = acr
						acr_hit_record[nc_id + protein]['evalue'] = evalue
				else:
					acr_hit_record[nc_id + protein] = {'record': record_dict[wp], 'evalue': evalue}
					acr_hit_record[nc_id + protein]['record'].id = '|'.join([nc_id, protein, slen, pident])
					acr_hit_record[nc_id + protein]['record'].description = acr
				
			else:
				if wp in acr_hit_record:
					if evalue < acr_hit_record[wp]['evalue']:
						acr_hit_record[wp]['record'].id = '|'.join([wp, slen, pident])
						acr_hit_record[wp]['record'].description = acr
						acr_hit_record[wp]['evalue'] = evalue
				else:
					acr_hit_record[wp] = {'record': record_dict[wp], 'evalue': evalue}
					acr_hit_record[wp]['record'].id = '|'.join([wp, slen, pident])
					acr_hit_record[wp]['record'].description = acr
	
	with open(ACR_INTERMEDIATES_FASTA_FILE, 'w') as out_handle:
		for wp_id in acr_hit_record:
			SeqIO.write(acr_hit_record[wp_id]['record'], out_handle, 'fasta')

	return acr_hit_record, ACR_INTERMEDIATES_FASTA_FILE



'''
	Functionality:
		Finalize the Acr homolog results and Write the results to the output file
	Args:
		acr_hit_record - the homolog sequence dict (protein id -> SeqRecord data (defined in Biopython https://biopython-cn.readthedocs.io/zh_CN/latest/cn/chr03.html ) )
		finalAcrs - list of list with only the candidate Acr/Aca's that had an index number found in uniqueHits.
		OUTPUT_DIR - the directory where output files are stored
		GCF - The prefix given to all output files
		isProdigalUsed - flag the .gff and .faa files generated by prodigal
	Returns:
		finalHomologFile - the result file of acr homolog search methods (contains acr homologs located in Acr-Aca locus)
'''
def finalizeHomolog(acr_hit_record, finalAcrs, OUTPUT_DIR, GCF, isProdigalUsed):
	finalHomologFile = OUTPUT_DIR + GCF + '_final_acr_homolog.fasta'
	out_handle = open(finalHomologFile, 'w')

	for locus in finalAcrs:
		locus_start = math.inf; locus_end = 0
		isAcrHitinLocus = False
		Acr_hit_homolog = list()
		Acr_Aca_cluster = list()
		nc_id = ''

		for protein in locus:
			locus_start = protein.start if protein.start < locus_start else locus_start
			locus_end = protein.end if protein.end > locus_end else locus_end

			Acr_Aca_cluster.append(str(protein.wp))
			nc_id = protein.nc
			if isProdigalUsed:
				if nc_id + protein.wp in acr_hit_record.keys():
					Acr_hit_homolog.append(nc_id + protein.wp)
					isAcrHitinLocus = True
			else:
				if protein.wp in acr_hit_record.keys():
					Acr_hit_homolog.append(protein.wp)
					isAcrHitinLocus = True

		if isAcrHitinLocus:
			for wp in Acr_hit_homolog:
				acr_hit_record[wp]['record'].description += ' {0}|{1}|{2}|{3}'.format(nc_id, '-'.join(Acr_Aca_cluster), locus_start, locus_end)

	for wp_id in acr_hit_record:
		SeqIO.write(acr_hit_record[wp_id]['record'], out_handle, 'fasta')

	out_handle.close()

	return finalHomologFile

