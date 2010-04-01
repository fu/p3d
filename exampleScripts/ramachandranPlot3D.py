#!/usr/bin/env python
'''Ramachandran Plot in 3D by Ch. Fufezan 2009

usage: ramachandranPlot3D.py <bin size for angles> <pdb file>'''

import sys, os, time
from p3d import protein as protein
from collections import defaultdict as ddict

if (__name__ == '__main__'):
	if (len(sys.argv) < 3):
		print (__doc__)
		sys.exit(1)
	binSize = float(sys.argv[1])
	if 180 % binSize != 0: sys.exit(1)
	ramachandranPlot3D = ddict(int)
	pdbs = sys.argv[2:]
	StartTime = time.time()
	for k,entry in enumerate(pdbs):
		print('Analysing',entry,end='...')
		try:
			pdb = protein.Protein(entry,DunbrackNaming=True,BSPTree=False)
			''' 
			query Set alphas to have one atom per residue from both chains
			'''
			for i,alpha in enumerate(pdb.query('alpha and model 1')):
				PhiPsi = alpha.calcPhiPsi(allowAlternativeConfs=False)
				ramachandranPlot3D[((PhiPsi[0][0]//binSize)*binSize,(PhiPsi[0][1]//binSize)*binSize)] += 1
				'''
				calcPhiPsi returns a list but allowAlternativeConfs= is set to False
				the list has only one element and thus we can directly access 
				phi and psi with PhiPsi[0][0] and PhiPsi[0][1] respectively ...
				'''
			print('done - ',k)
		except:
			print('FAILED - ',k)
	for coords,value in ramachandranPlot3D.items():
		print(coords[0],'\t',coords[1],'\t',value)