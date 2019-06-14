#!/usr/bin/python3
'''
	Prints formatted output given organism information.
'''
from organism import *


'''
	Returns an faa string.

	Params:
		proteinList -  LIST. Contains protein objects
		delimiter - STRING. How to seperate the header elements of each protein
'''
def get_faa(proteinList, delimiter = '|'):
	faa = ""
	for protein in proteinList:
		headerElements = list()
# self, wp, nc, sequence, position, start, end
		headerElements.extend((protein.nc, protein.wp, str(protein.start), str(protein.end), str(protein.position)))
		header = ">{0}".format(delimiter.join(headerElements))
		niceSequence = re.sub("(.{80})", "\\1\n", protein.sequence, 0, re.DOTALL)	# Maximum sequence length is 80
		faa += '\n'.join((header, niceSequence))
		if not faa.endswith('\n'):
			faa += '\n'

	return faa





def faa_a_sequence(sequence, headerElements, delimiter = '|'):
	header = delimiter.join(headerElements)
	if not header.startswith('>'):
		header = '>' + header

	faa = header + '\n'
	niceSequence = re.sub("(.{80})", "\\1\n", sequence, 0, re.DOTALL)   # Maximum sequence length is 80
	faa += niceSequence
	faa += '\n'

	return faa




'''
	Writes faa file to file

	Params:
		proteinList - LIST. Contains protein objects
		outFile - STRING. Where to save faa sequence
		delimiter - STRING. How to seperate the header elements of each protein
'''
def write_faa(proteinList, outFile, delimiter = '|'):
	d = delimiter
	with open(outFile, 'w') as handle:
		handle.write(get_faa(proteinList, delimiter = d))
