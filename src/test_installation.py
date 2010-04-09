#!/usr/bin/env python
'''
Testing the installation of p3d
'''
import sys, p3d.protein


TESTFILES = ['../pdbs/1BPH.pdb.gz']


if (__name__ == '__main__'):
	for pdb in TESTFILES:
		protein = p3d.protein.Protein(pdb)
		for line in protein.info():
			print(line)