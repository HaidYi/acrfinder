#!/usr/bin/python3
'''
	*****************************************************************************************************
	Purpose:	Parses Aca/Acr result file with the aide of gi/pai database information.
				Limits the number of Aca/Acr results to only include the ones found near (lax) or strictly within (strict) either a known gi or a pai.
				Supports various options so user can specifiy the type of output (lax/strict) and which database to use.

	Author:
		Javi Gomez - https://github.com/rtomyj
	*****************************************************************************************************
'''
from collections import defaultdict

'''
	Purpose:	Judging whether two intervals overlap

	Arguments:	Interval_1: the first interval (start, end)
				Interval_2: the second interval (start, end)
'''
def isItvOverlap(interval_1, interval_2):
	return max(interval_1[0], interval_2[0]) < min(interval_1[1], interval_2[1])

'''
	Purpose:	Parses GI DB file.
				Creates a dict using unique ncid's to map to start and end positions of regions found in GI DB.

	Arguments:	GI_NC_ID_maps_START_END = defaultdict with an empty list as its default value. Stores ncids found in GI DB and all the GI's found within that NCID
				GI_DB_FILE = file name of the GI DB file

	Returns:	None
'''
def parse_gi_file(GI_NC_ID_maps_START_END, GI_DB_FILE):
	with open(GI_DB_FILE, 'r', 512) as handle:
		'''
			Cycles through GI DB file
		'''
		for line in handle:
			if line.startswith('Accession_number'):
				continue

			cols = line.split('\t')
			ncid, start, end = cols[0], cols[1], cols[2]	# extracts an assumed NC ID and its start and end position
			tupStartEnd = (int(start), int(end))	# converts start and end from string to int and makes a tuple out of it
			GI_NC_ID_maps_START_END[ncid].append(tupStartEnd)	# appends start/end tuple to ncid dict with using the tuples ID as key



'''
	Purpose:	Parses PAI DB file
				Creates a dict using unique ncid's to map to start and end positions of regions found in PAI DB.

	Arguments:	PAI_NC_ID_maps_START_END = defaultdict with an empty list as its default value. Stores ncids found in PAI DB and all the PAI's found within that NCID
				PAI_DB_FILE = file name of the PAI DB file

	Returns:	None
'''
def parse_pai_file(PAI_NC_ID_maps_START_END, PAI_DB_FILE):
	from regex import Regex	# contains often used regex's
	with open(PAI_DB_FILE) as handle:
		'''
			Cycles through PAI DB file.
		'''
		for line in handle:
			ncid = ""

			match = Regex.NC_REGEX.search(line)	# tries to find NC ID in line
			if match != None:	# if NC regex matches something in the line then extraction of start and end can continue
				ncid = match.group(0)

				match = Regex.START_END_REGEX2.search(line) 	# uses a regex to find the start and end position of PAI region
				start, end = int(match.group(1)), int(match.group(2))	# obtains start and end from regex match and converts from string to int
				startEnd = (start, end)	# creates tuple of start and end
				PAI_NC_ID_maps_START_END[ncid].append(startEnd)	# appends tuple to dict



'''
	Purpose:	Reads in Acr/Aca results containted in IN_FILE.
				Picks out only the loci that are near or within (lax/strict) GI or PAI DB.

	Arguments:	IN_FILE = Acr/Aca results file user wants parsed
				OUT_FILE = file name to use for output (might not be exactly as they entered)
				NC_ID_maps_START_END = dictionary contianing ncids found either in GI or PAI
				uniqueNcids = unique ncids found either in GI or PAI
				PRINT_STRICT = whether user wants to strict results
				PRINT_LAX = whether user wants to lax results
				_type = string identifying which DB is being used eg) 'PAI', 'GI'

	Returns:	None
'''
def limit_acr(candidateAcrs, NC_ID_maps_START_END, uniqueNcids, PRINT_STRICT, PRINT_LAX, queue, _type):
	'''
		Using Interval Overlapping Algorithm to decide the Acr hit against 'PAI' and 'GI' Database.
	'''
	for position, locus in enumerate(candidateAcrs):
		locusNcid = locus[0].nc  # nc id of the locus
		locus_start_end_list = list()
		
		'''
			We only care about loci found in a DB
			Therefore if the ncid of the loci doesn't appear in the keys of the dict that contains unique ncids then we don't need to parse that locus. It won't be in a region the DB identified.
		'''
		
		if locusNcid in uniqueNcids:
			startEndList = NC_ID_maps_START_END[locusNcid]

			isLaxHit, isStrictHit = False, False
			locus_hit_set = set()
			
			for protein in locus:
				locus_start_end_list.append( (protein.start, protein.end) )  # add (protein.start, protein.end) into list
			
			startEndList.extend( locus_start_end_list )  # extend the locus start-end list
			
			# Sort the list using the start of different intervals
			startEndList = sorted(startEndList, key = lambda x: x[0])

			# put protein (in locus) hitted in the db into locus_hit_set list
			for index, _ in enumerate(startEndList[:-1]):
				if isItvOverlap(startEndList[index], startEndList[index + 1]):
					isLaxHit = True
					if startEndList[index] in locus_start_end_list:
						locus_hit_set.add(startEndList[index])
					if startEndList[index + 1] in locus_start_end_list:
						locus_hit_set.add(startEndList[index + 1])

			# isStrictHit = true if all proteins in the loci are hitted in the db
			if len(locus) == len(locus_hit_set):
				isStrictHit = True
			
			if PRINT_STRICT: # Strict Mode
				if isStrictHit:
					lociHits = queue.get()
					lociHits.add(position)
					queue.put(lociHits)
			elif PRINT_LAX: # Lax Mode
				if isLaxHit:
					lociHits = queue.get()
					lociHits.add(position)
					queue.put(lociHits)




'''
	Purpose:	Builds dict of ncids that contain GI's using the GI DB file.
				Uses dict to limit the acr results of input file.

	Arguments:	IN_FILE = Acr/Aca results file user wants parsed
				OUT_FILE = file name to use for output (might not be exactly as they entered)
				PRINT_STRICT = whether user wants to strict results
				PRINT_LAX = whether user wants to lax results
				queue - contains synchronized data (such as a list of loci indices that should be kept)

	Returns:	None
'''
def use_gi_db_on_acr(candidateAcrs, PRINT_STRICT, PRINT_LAX, GI_DB_FILE, queue):
	GI_NC_ID_maps_START_END = defaultdict(lambda: list())
	parse_gi_file(GI_NC_ID_maps_START_END, GI_DB_FILE)	# parses GI DB file to build GI dict
	mobilomeNcids = list(GI_NC_ID_maps_START_END.keys())
	print('Total number of ncids found in gi db: ' + str(len(mobilomeNcids)))	# prints number of unique ncids found in GI file

	limit_acr(candidateAcrs, GI_NC_ID_maps_START_END, mobilomeNcids, PRINT_STRICT, PRINT_LAX, queue, 'gi')	# parses the input file and creates new files (lax/strict) with only the results that lie within/strictly in a known GI



'''
	Purpose:	Builds dict of ncids that contain GI's using the GI DB file.
				Uses dict to limit the acr results of input file.

	Arguments:	IN_FILE = Acr/Aca results file user wants parsed
				OUT_FILE = file name to use for output (might not be exactly as they entered)
				PRINT_STRICT = whether user wants to strict results
				PRINT_LAX = whether user wants to lax results
				queue - contains synchronized data (such as a list of loci indices that should be kept)

	Returns:	None
'''
def use_pai_db_on_acr(candidateAcrs, PRINT_STRICT, PRINT_LAX, PAI_DB_FILE, queue):
	PAI_NC_ID_maps_START_END = defaultdict(lambda: list())
	parse_pai_file(PAI_NC_ID_maps_START_END, PAI_DB_FILE)	# parses PAI DB file to build PAI dict
	paiNcids = list(PAI_NC_ID_maps_START_END.keys())
	print('Total number of ncids found in PAI db: ' + str(len(paiNcids)))	# prints number of unique ncids found in PAI file

	limit_acr(candidateAcrs, PAI_NC_ID_maps_START_END, paiNcids, PRINT_STRICT, PRINT_LAX, queue, 'pai')	# parses the input file and creates new files (lax/strict) with only the results that lie within/strictly in a known PAI

