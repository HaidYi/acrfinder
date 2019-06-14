'''
	*****************************************************************************************************
	Purpose:
		Contains methods used to parse and handle files produced by acr_aca_finder.py

	Author:
		Javi Gomez - https://github.com/rtomyj wrote the version 1.0.
		Haidong Yi - https://github.com/haidyi revised the codes, remove some bugs, wrote the version 2.0. 
	Project:
		This work is advised by Dr. Yanbin Yin at UNL - yyin@unl.edu
	*****************************************************************************************************
'''

'''
	# GENERATOR #
	Background:
		A locus is a region that contains Acr/Aca proteins and satisfies all 4 filters used to find Acr/Aca loci.
		Primarily, proteins in a locus are close to each other.
		A loci is more than one locus.
	Summary:
		Parses file line by line to obtain a locus. Once an empty line is reached a locus has ended.
		Yields next locus of handle to caller.
	Params:
		Handle - file handle to parse
	Yields:
		proteins of a locus where each protein is a line of the main file.
'''
def getLocus(handle):
	locus = ""
	for line in handle:
		if line.strip() == "" and  locus != "":	# blank line means there are no more proteins in that locus
			yield locus.strip()	# yeild locus to caller
			locus = ""	# reset
		else:
			locus += line	# append new protein to locus



'''
	Summary:
		Obtains the starting BP and the ending BP of a given locus.
	Params:
		locusProteins - string that contains all proteins in a locus where each line contains one protein.
	Returns:
		(starting BP of loci), (ending BP of loci)
'''
def getLocusStartAndEnd(locusProteins:str):
	start = locusProteins[0].split('\t')[3]
	end = locusProteins[len(locusProteins) - 1].split('\t')[4]
	return int(start), int(end)
