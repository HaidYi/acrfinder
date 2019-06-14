#/usr/bin/python3
'''
	Class holds often used regular expressions.
	All patterns are stored as class variables not instance variables which means Regex.(REGEX_NAME) will work and an object doesn't need to be created.
'''

import re



class Regex():
	'''
		Regex's to use in gbff file
	'''
	TRANSLATION_REGEX = re.compile(r'translation="([A-Z\s]+)"')	# finds the sequence of gene within CDS. group[1] returns just the sequence
	PROTEIN_ID_REGEX = re.compile(r'protein_id="([A-Z0-9_.]+)"')	# finds the protein id (WP id) within CDS. .group[1] returns just the protien id
	START_END_REGEX1 = re.compile(r'([0-9]+)\.\.([0-9]+)')	# finds the start and end position of an element if it has the patter #..#. group[1] returns start, group[2] returns end


	'''
		Regex's to find GCF number, NC ID and WP ID. Captures just the ID (minus .[versionnumber]) in group 1 of match (see re.search()). 
	'''
	NC_REGEX = re.compile(r'([AN][CGMRTWZ]_[A-Z0-9]+)\.[0-9]')    # finds nc id regardless of first and second letter. group[1] returns just the nc id without the version/accesnision number (no .#)
	GCF_REGEX = re.compile(r'(GC[FA]_[0-9]+)\.[0-9]')	#finds the gcf id for an organism. group[1] returns just the gcf id without .#
	WP_REGEX = re.compile(r'[ANWY]P_[0-9]+\.[0-9]')


	'''
		Other regex's
	'''
	START_END_REGEX2 = re.compile('([0-9]+)-([0-9]+)')	# start and end position are seperated by a dash


def string_with_limited_width(subject, width=80):
	return re.sub('(.{' + str(width) + '})', '\\1\n', subject, 0, re.DOTALL)
