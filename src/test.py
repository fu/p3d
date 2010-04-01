#!/usr/bin/env python
# encoding: utf-8
# Created by Ch. Fufezan on 2009-02-24.
# Copyright (c) 2009 Ch.Fufezan. All rights reserved.
#

"""
test.py

small test script for p3d
Change the Test_Set dictionary in this script to point to some pdb files :)
"""


"""
This file is part of p3d.

    p3d is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    p3d is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys

import os
import p3d.protein
import p3d.vector

Test_Set = {'1MZ4':'../pdbs/alternativeConfResid74_1MZ4.pdb'}
#, '7GAT':'../pdbs/7GAT.pdb','1EN2':'../pdbs/1EN2.pdb','2AXT':'../pdbs/2AXT.pdb.gz'}

TESTS_4_PROTEIN = False
TESTS_4_ATOMS = False
TESTS_4_PARSER = True
TESTS_4_PHIPSI = False
TESTS_4_QUERY = True

def main():
	for name,FileLocation in Test_Set.items():
		print('\n\n---Example Code---')
		a = p3d.protein.Protein(FileLocation)
		'''
		'''
		if TESTS_4_PROTEIN:
			print("\n----*---- Protein ----*----\n")
			info = a.info()
			for line in info:
				print(line)
			print(a.stats['model'])
			#a.writeToFile(a.fullname)
		'''
		'''
		if TESTS_4_ATOMS:
			print("\n----*---- Atoms ----*----\n")
			queries = [
			'resname LEU and resid 21',
			'resname LEU and resid 21 and oxygen'
			]
			for query in queries:
				n = a.query(query)
				print(query,':',len(n),'{0}'.format('atom' if len(n) ==1 else 'atoms'))
				surrounding = a.collectSphereAtoms(centre=n,radius=5)
				print('Surrounding:',len(surrounding))
				print('LookUp Routine ...')
				t = a.lookUpAtom(query)
				print(type(t))
			print('-------------')
			print(t.info(lvl='max'))
			residueAtoms = t.allAtomsOfSameResidue()
			for atom in residueAtoms:
				print(atom.info(lvl='max'))
		'''
		'''
		if TESTS_4_PARSER:
			print("\n----*---- Parser ----*----\n")
			queries = [\
			'reSname ALA',\
			'resname ALA and resid > 40',\
			'chain A aNd resid 13..25 and not bKb',\
			'resname GLU and not bkb and not oXYgen and resid 10 tO 90 and nOt atype CB',\
			'not protein and oxygen and not resname HOH',\
			'within 3 of (oxygen and resid 1)',\
			'fiRSt resiDUE of chain A',\
			'last residue of chain A',\
			]
			for i,query in enumerate(queries):
				print('Test #', i, '>>>', query)
				atoms = a.query(query)
				if len(atoms) == 0:
					print( "No atoms found" )
				else:
					for atom in atoms:
						print(atom.info(lvl='max'))
		if TESTS_4_PHIPSI:
			print("\n----*---- Phi Psi ----*----\n")
			for atom in a.query("backbone and oxygen and (resid 1..10 or resid > 130)"):
				print(atom.info(),atom.calcPhiPsi())
		'''
		'''
		if TESTS_4_QUERY:
			LowerNumber = 40
			distance = 3
			k = p3d.vector.Vector(40, 10 , -10)
			CustomSet = [ p3d.vector.Vector(40, 10 , -10) , p3d.vector.Vector(48, 15, -16) , p3d.vector.Vector(51, 25,-10) ]
			queries = [\
			'resname ALA',\
			'resname ALA and resid > 40',\
			'chain A and resid 13..25 and not bkb',\
			'resname GLU and not bkb and not oxygen and resid 10 to 90 and not atype CB',\
			'not protein and oxygen and not resname HOH',\
			'resname ALA, GLU, TYR',
			'residue id 40,50,60',
			'protein and not resname ALA, GLU, TYR',
			'protein',
			'first residue of chain A and not oxygen',
			'last residue of chain A',
			
			]
			for i,query in enumerate(queries):
				print("<source lang='python'>")
				print(">>> atoms = pdb.query('",query,"')")
				atoms = a.query(query)
				print(">>> print(len(atoms),' atoms found')")
				print(len(atoms)," atoms found")
				print(">>> print(list(atoms)[0].info(lvl='max'))")
				print(list(atoms)[0].info(lvl='max'))
				'''
				if len(atoms) != 0:
					print(">>> print('Sample atom:\n',list(atoms)[0].info(lvl='max'))")
					print('Sample atom:\n',list(atoms)[0].info(lvl='max'))
				'''
				print("</source>\n")
			## ---- Variable example ----
			print("<source lang='python'>")
			print(">>> LowerNumber = 40")
			print(">>> atoms = pdb.query('resid <',str(LowerNumber))")
			atoms = a.query('resid < ',str(LowerNumber))
			print(">>> print(len(atoms),' atoms found')")
			print(len(atoms)," atoms found")
			print(">>> print(list(atoms)[0].info(lvl='max'))")
			print(list(atoms)[0].info(lvl='max'))
			print("</source>\n")
			## ---- function examples first of chain A & within distance example
			print("<source lang='python'>")
			print(">>> distance = 3")
			print(">>> atoms = pdb.query('first residue of chain A and within',str(distance),'of oxygen')")
			atoms = a.query('first residue of chain A and within ',str(distance),' of oxygen')
			print(">>> print(len(atoms),' atoms found')")
			print(len(atoms)," atoms found")
			print(">>> print(list(atoms)[0].info(lvl='max'))")
			print(list(atoms)[0].info(lvl='max'))
			print("</source>\n")
			#
			for i in range(1,10):
				atoms = a.query('not protein and within ',str(i),' of oxygen')
				print(len(atoms)," {0} found".format('atom' if len(atoms) == 1 else 'atoms'))
			for i in range(1,10):
				atoms = a.query('protein and within ',str(i),' of',k)
				print(len(atoms)," atoms found")
			for i in range(1,10):
				atoms = a.query('protein and within',str(i),' of',CustomSet)
				print(len(atoms)," atoms found")
			#
			atoms = a.query('bkb and chain and resid {0}'.format(a.firstResidueOfChain('A',idOnly=True)))
	return	

if __name__ == '__main__':
	main()

