#!/usr/bin/python3
'''
	*****************************************************************************************************
	Purpose:
		Contains methods to find crispr spacers from a result.json file created by CRISPRCasFinder.
		Obtains the sequence of the crispr spacers and makes a fasta file containing all the spacers found in result.json

	Author:
		Javi Gomez - https://github.com/rtomyj wrote the version 1.0.
		Haidong Yi - https://github.com/haidyi revised the codes, remove some bugs, wrote the version 2.0. 
	Project:
		This work is advised by Dr. Yanbin Yin at UNL - yyin@unl.edu
	*****************************************************************************************************
'''
from sys import path as sys_path
sys_path.append('dependencies/PyGornism/')
from formated_output import faa_a_sequence



'''
	Purpose:
		Parses result.json file (output from CRISPRCasFinder) from given directory.
		Finds CRISPR spacers and writes them to a new file in fasta format.
	Arguments:
		CRISPR_CAS_OUTPUT -	Where result.json is located
		SELECT_SPACER_FASTA -	path to fasta file contianing spacers
		EVIDENCE_LEVEL -	Evidence level of spacers desired
	Returns:
		None
'''
def parse_result_json(CRISPR_CAS_OUTPUT: str, SELECT_SPACER_FASTA: str, EVIDENCE_LEVEL: int):
	from json import load

	faa = ''	# temp variable holding fasta string
	arraysWithWantedEvidence, totalCrisprs = 0, 0	# stats
	hasCRISPRCasSystem = False	# CRISPR Cas = has Array and Cas system
	try:
		with open(CRISPR_CAS_OUTPUT + '/result.json', 'r') as inHandle, open(SELECT_SPACER_FASTA, 'w', 512) as outHandle:
			JSON = load(inHandle)	# loads result json file
			for sequence in JSON['Sequences']:	# traverses ncid's
				casSystemTypes = []	# list to hold all Cas system names and their position
				hasCas, hasCRISPRArrays = False, False
				numCas, numCrisprs = len(sequence['Cas']), len(sequence['Crisprs'])

				print('\nAnalyzing sequence ' + sequence['Id'])
				print('Cas found - {0}, Crisprs found - {1}\n'.format(numCas, numCrisprs))


				'''
					See's if there is a CRISPR Cas system anywhere in organism
				'''
				if numCas > 0:
					hasCas = True
				if numCrisprs > 0:
					hasCRISPRArrays = True

				if hasCas and hasCRISPRArrays:
					hasCRISPRCasSystem = True

				totalCrisprs += len(sequence['Crisprs'])

				if not hasCas and hasCRISPRArrays:	# only gathers crisprs with cas systems
					continue


				'''
					Finds all Cas systems in current sequence (ncid)
				'''
				for casSystem in sequence['Cas']:
					startEnd, systemType = str(casSystem['Start']) + '-' + str(casSystem['End']), casSystem['Type']
					if systemType == "":
						systemType = 'Inconclusive Cas System'

					casSystemTypes.append([startEnd, systemType])	# appends Cas system type and its start and end to a running list


				casSystemHeaderStr = ''	# header info that contains all Cas systems found in this sequence
				'''
					Cycles through all Cas systems found.
					Appends all systems found to header.
				'''
				for cas in casSystemTypes:
					startEnd, casSystem = cas
					if casSystemHeaderStr != "":
						casSystemHeaderStr += '+'
					casSystemHeaderStr = casSystemHeaderStr + '{0}({1})'.format(casSystem, startEnd)	# appends new system to header


				'''
					Cycles through all CRISPR arrays.
				'''
				for crispr in sequence['Crisprs']:
					crisprEvidence = int(crispr['Evidence_Level'])

					if crisprEvidence < EVIDENCE_LEVEL:
						continue
					arraysWithWantedEvidence += 1	# running total of all arrays that have an evidence level equal to or greater than the wanted evidence level

					arrayStart, arrayEnd = crispr['Start'], crispr['End']

					_id = crispr['Name']
					for region in crispr['Regions']:
						if region['Type'] == 'Spacer':
							'''
								Obtains information about all the spacers.
							'''
							sequence = region['Sequence']
							start = region['Start']
							end = region['End']
							headerElements = [str(_id), 'Array_Start={0}|Array_End={1}|Evidence_Level={2}|Spacer_Start={3}|Spacer_end={4}'.format(arrayStart, arrayEnd, crisprEvidence, start, end), casSystemHeaderStr]
							faa += faa_a_sequence(sequence, headerElements)
			outHandle.write(faa)
	except IOError as e:
		'''
			Catches file not found error and error that can occur when parsing json file.
			result.json file might have a double quote between a string. This will cause an exception.
			THIS IS BECAUSE CRISPRCASFINDER DOESN'T ESCAPE THE QUOTES IN THE STRING. ESCAPE THE STRING IN THE FILE AND RERUN IF THE EXCEPTION HAPPENS.
		'''
		print('IOError ' + str(e))
	except ValueError as e:
		print(CRISPR_CAS_OUTPUT + '/result.json', 'has err:', e, 'CHECK TO SEE IF STRING HAS A DOUBLE QUOTE, IF SO ESCAPE IT USING UNIX COMMAND OR OTHER RESOURCES')

	print('Total CRISPRCas spacers with evidence level {0} or greater: {1}/{2}'.format(EVIDENCE_LEVEL, arraysWithWantedEvidence, totalCrisprs))


	'''
		If no CRISPR Cas systems were found then we cannot continue the finding of Acr/Aca proteins.
		Terminate program.
	'''
	if not hasCRISPRCasSystem:
		print('No CRISPRCas systems found. Terminating...')
		exit(0)



'''
	Purpose:
		Creates spacer fasta file.
	Arguments:
		CRISPR_CAS_OUTPUT =	directory containing all output files created by CRISPRCasFinder
		INTERMEDIATES - directory to store intermediate files
		EVIDENCE_LEVEL -	Evidence level of spacers desired
		SELECT_SPACER_FASTA - file to store spacers in fasta format
	Returns:
		SELECT_SPACER_FASTA - str, path to fasta file containing CRIPSR spacers with desired evidence level.
'''
def fastafy_select_spacers(CRISPR_CAS_OUTPUT, INTERMEDIATES, EVIDENCE_LEVEL, SELECT_SPACER_FASTA='spacers_with_desired_evidence.fna'):
	from os import path as os_path

	SELECT_SPACER_FASTA = INTERMEDIATES + SELECT_SPACER_FASTA

	print('Using evidence level of {0} to parse CRISPRCasFinder results found here -> {1}'.format(
		EVIDENCE_LEVEL, os_path.abspath(CRISPR_CAS_OUTPUT)))

	parse_result_json(CRISPR_CAS_OUTPUT, SELECT_SPACER_FASTA, EVIDENCE_LEVEL)

	if os_path.getsize(SELECT_SPACER_FASTA) == 0:
		print('No spacers found with evidence level {0}'.format(EVIDENCE_LEVEL))
		print('Proceeding with Acr/Aca identification.\n\n')
		return None

	return SELECT_SPACER_FASTA

