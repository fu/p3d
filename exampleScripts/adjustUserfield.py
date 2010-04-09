#!/usr/bin/env python
'''set user field of certain residues to 1.00 by Ch. Fufezan 2010

usage: adjustUserfield.py <pdb file> "selection-query-string" '''

import sys 
from p3d import *

if (__name__ == '__main__'):
	if (len(sys.argv) != 3):
		print (__doc__)
		sys.exit(1)
	pdb = protein.Protein(sys.argv[1])
	for atom in pdb.query(sys.argv[2]):
		atom.user = 1.00
		print(atom.info(lvl='max'))