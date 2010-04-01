#!/usr/bin/env python
'''Extracts part of a PDB file using p3d query methods by Ch. Fufezan 2009

usage: PDBextract.py <pdb file> <"querystring">
e.g.   PDBextract.py pdb.pdb "residue 1 and backbone"

see online documentation for the syntax of the query string
'''
import sys, p3d.protein

if (__name__ == '__main__'):
	if (len(sys.argv) != 3):
		print (__doc__)
		sys.exit(1)
	atomSet = p3d.protein.Protein(sys.argv[1]).query(sys.argv[2])
	for atom in atomSet:
		print(atom.output()) 