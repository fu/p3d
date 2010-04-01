#!/usr/bin/env python
'''set beta field to RMSD of residues in different NMR models by Ch. Fufezan 2009

usage: NMR_rmsd.py <pdb file>'''

import sys, os, time
from p3d import protein as protein
import p3d.vector
from collections import defaultdict as ddict

if (__name__ == '__main__'):
	if (len(sys.argv) != 2):
		print (__doc__)
		sys.exit(1)
	pdb = protein.Protein(sys.argv[1],BSPTree=False)
	''' 
	query Set alphas to have one atom per residue from both chains
	'''
	for i,alpha in enumerate(pdb.query('alpha and model 1')):
		allAtomsOfResidue = alpha.allAtomsOfSameResidue()
		#print(alpha.output())
		#print('Residue',alpha.resid,'counted Atoms',len(allAtomsOfResidue))
		geoCentres = []
		meanGeoCentre = p3d.vector.Vector()
		for model,hashed in pdb.hash['model'].items():
			zeroVector = p3d.vector.Vector()
			for atom in (allAtomsOfResidue & hashed):
				zeroVector += atom
			#print('model',model,'counted Atoms',len(allAtomsOfResidue & hashed),'centre of geo',zeroVector.info())
			geoCentres.append(zeroVector)
			meanGeoCentre += zeroVector
		meanGeoCentre = meanGeoCentre/len(pdb.hash['model'].keys())
		#print('mean geo centre',meanGeoCentre.info())
		rmsd = 0
		for geoCentre in geoCentres:
			#print(rmsd)
			d = meanGeoCentre.distanceTo(geoCentre)
			if d > rmsd:
				rmsd = d
		alpha.beta = rmsd
		#print(alpha.output())
	pdb.writeToFile(pdb.fullname+'_rmsdOnBeta.pdb')
