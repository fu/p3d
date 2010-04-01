#!/usr/bin/env python
'''contact surface by Ch. Fufezan 2009

usage: contactSurface.py <pdb file> <chain1> <chain2> <maxDistance to take into account>'''

import sys, os
from p3d import protein as protein
from collections import defaultdict as ddict

def estimateMinContacts(Protein=None,ResidueAtom=None,targetChain=None,DistanceMatrix=None,maxDistance=None):
	residueAtoms = Protein.query('protein and resid {0} and chain {1}'.format(ResidueAtom.resid,ResidueAtom.chain))
	surounding = Protein.query('protein and within {0} of'.format(maxDistance),residueAtoms,'and chain {0}'.format(targetChain))
	for residAtom in residueAtoms:
		for atom in surounding:
			d = atom.distanceTo(residAtom)
			ordered_tuple = (residAtom.resid,atom.resid) if targetChain == 'B' else (atom.resid,residAtom.resid)
			if DistanceMatrix[ordered_tuple] > d or DistanceMatrix[ordered_tuple] == 0:
				DistanceMatrix[ordered_tuple] = d
	return DistanceMatrix

if (__name__ == '__main__'):
	if (len(sys.argv) != 5):
		print (__doc__)
		sys.exit(1)
	chain_A , chain_B, maxDistance = sys.argv[2:5]
	pdb = protein.Protein(sys.argv[1],chains=[chain_A,chain_B])
	DistanceMatrix = ddict(float)
	''' query Set alphas to have one atom per residue from both chains '''
	chain_A_alphas = pdb.query('chain {0} and alpha'.format(chain_A))
	chain_B_alphas = pdb.query('chain {0} and alpha'.format(chain_B))
	for alpha in chain_A_alphas:
		DistanceMatrix = estimateMinContacts(Protein=pdb,ResidueAtom=alpha,targetChain=chain_B,DistanceMatrix=DistanceMatrix,maxDistance=maxDistance)
	for alpha in chain_B_alphas:
		DistanceMatrix = estimateMinContacts(Protein=pdb,ResidueAtom=alpha,targetChain=chain_A,DistanceMatrix=DistanceMatrix,maxDistance=maxDistance)
	for coords,value in DistanceMatrix.items():
		print(coords[0],'\t',coords[1],'\t',value)