#!/usr/bin/env python
'''
Testing the installation of p3d
'''
import sys, p3d
from p3d import *

TESTFILES = ['../pdbs/1BPH.pdb.gz']

def module_info():
	print('p3d module information:')
	print('  * Version {0}'.format(p3d.__version__.__version__))
	print('  * Installation path {0}\n'.format(p3d.__file__))
	return

def protein_info():
	for pdb in TESTFILES:
		protein = p3d.protein.Protein(pdb)
		for line in protein.info():
			print(line)
	return

def main():
	module_info()
	protein_info()
	return
	
if (__name__ == '__main__'):
	main()
	