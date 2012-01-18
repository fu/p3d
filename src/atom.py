#!/usr/bin/env python
# encoding: utf-8
# Copyright (c) 2009 Ch. Fufezan.
# module is distributed under the terms of the GNU General Public License V2
# See LICENSE for more details.


'''p3d - a protein structure module for python 
'''


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

from collections import defaultdict as ddict
from copy import deepcopy as dcp
import random as r
from math import pi

import p3d.vector
import p3d.library
import p3d.geo

class InputError(Exception): pass
class NotAminoAcidError(InputError): pass
class Not3LetterAA(InputError): pass
class Not1LetterAA(InputError): pass
class NotOnlyOneAtom(InputError): pass
class MissingAtomForOperation(InputError): pass

class Atom(p3d.vector.Vector):
    '''
    Atoms are defined by their line in pdb file
    '''
    def __init__(self,line=None,protein=None,PositionInAtomslist=None,model=1,matchObject=None):
        '''
         ---6-|--5-|-4-|2|-3|2|-4-|2|---10----|-----8-|----8--|--6--|---6-|
         HETATM
         ATOM    559  CA BASP A  74      48.780  13.254  -1.818  0.50 16.34           C  
         ----.----|----.----|----.----|----.----|----.----|----.----|----.----|----.----|----.----
             5    10   15   20   25   30   35   40   45   50   55   60   65   70
         ATOM   9982  CZ  PHE D  27      14.293  79.865  39.022  1.00 85.37           C
         HETATM
        '''     
        """
        # > Old slicing style
        self.x       = float(line[30:38])
        self.y       = float(line[38:46])
        self.z       = float(line[46:54])
        self.idx     = int(line[6:11])
        self.atype   = str(line[11:16])
        self.altConf = str(line[16:17].upper()) if str(line[16:17].upper()) != ' ' else '_' # ' ', A, B and so on
        self.aa      = str(line[17:20])
        self.chain   = str(line[21:22])
        self.resid   = int(line[22:26])
        self.altConf2= str(line[26:27].upper())
        self.desc    = ''
        self.beta    = float(line[60:66])
        self.user    = float(line[54:60])
        self.type    = line[:6]
        """
        #self.mObject = matchObject
        self.type    = matchObject.group('type')
        self.idx     = int(matchObject.group('index'))
        self.atype   = matchObject.group('atype')
        self.altConf = matchObject.group('altconf').upper() if matchObject.group('altconf').upper() != ' ' else '_' # ' ', A, B and so on
        self.aa      = matchObject.group('aa') # .replace(' ','_')
        self.chain   = matchObject.group('chain')
        self.resid   = int(matchObject.group('resid'))
        self.altConf2= matchObject.group('altconf2').upper()
        self.x       = float(matchObject.group('x'))
        self.y       = float(matchObject.group('y'))
        self.z       = float(matchObject.group('z'))
        self.beta    = float(matchObject.group('beta')) if matchObject.group('beta').strip() else 1.0
        self.user    = float(matchObject.group('user')) if matchObject.group('user').strip() else 0.0
        self.elementType = matchObject.group('elementType')
        self.charge  = matchObject.group('charge')
        """
        ltype       = line[:6]
        idx         = int(line[6:11])
        atype       = str(line[11:16])
        altConf     = str(line[16:17].upper()) if str(line[16:17].upper()) != ' ' else '_'# ' ', A, B and so on
        aa          = str(line[17:20])
        chain       = str(line[21:22])
        resid        = int(line[22:26])
        altConf2    = str(line[26:27].upper())
        x            = float(line[30:38])
        y            = float(line[38:46])
        z            = float(line[46:54])
        beta        = float(line[60:66])
        user         = float(line[54:60])

        checkas = [\
        (ltype  ,self.type)     , 
        (idx        ,self.idx)   ,
        (atype  ,self.atype)     ,
        (altConf ,self.altConf) ,
        (aa         ,self.aa ),
        (chain  ,self.chain) ,   
        (resid  ,self.resid) ,   
        (altConf2,self.altConf2),
        (x      ,self.x)         ,
        (y      ,self.y )    ,
        (z      ,self.z )    ,
        (beta   ,self.beta) ,    
        (user   ,self.user) ,    

        ]
        for (a,b) in checkas:
            if a != b:
                print '>old>'+str(a)+'<>'+str(b)+'<new<'
                exit(1)
        """

        self.desc    = ''
        self.protein = protein
        self.model   = model
        # --*-- PositionInAtomslist is used for history tag during vector operation --*--
        self.pos_in_list = PositionInAtomslist
        #print(line+'\n')
        #print(self.output(format='TestPhase'))
        return

    def calcPhiPsi(self,prescision=0,allowAlternativeConfs=True):
        '''
        Returns tuple (phi,psi) in degrees for given residue.
        Input can be any amino acid atom. 
        '''
        a = self.protein
        #if self.aa not in p3d.library.AA3:
        #   raise NotAminoAcidError()
        Np1 = list(a.hash['resid'][self.resid+1] & a.hash['atype']['N']  & a.hash['chain'][self.chain] & a.hash['model'][self.model])
        C   = list(a.hash['resid'][self.resid]   & a.hash['atype']['C']  & a.hash['chain'][self.chain] & a.hash['model'][self.model])
        CA  = list(a.hash['resid'][self.resid]   & a.hash['atype']['CA'] & a.hash['chain'][self.chain] & a.hash['model'][self.model])
        N   = list(a.hash['resid'][self.resid]   & a.hash['atype']['N']  & a.hash['chain'][self.chain] & a.hash['model'][self.model])
        CB4 = list(a.hash['resid'][self.resid-1] & a.hash['atype']['C']  & a.hash['chain'][self.chain] & a.hash['model'][self.model])
        query_sets_phi = [{'C':a,'CA':b,'N':c,'CB4':d} for a in C for b in CA for c in N for d in CB4]
        query_sets_psi = [{'Np1':a,'C':b,'CA':c,'N':d} for a in Np1 for b in C for c in CA for d in N]
        psi = -777
        phi = -777
        for omg_set in [query_sets_phi,query_sets_psi]:
            if len(omg_set) == 0 :
                omg_set.append({'dummy':''})
        if allowAlternativeConfs:
            allpossibleCombos =  [(phia, psia) for phia in query_sets_phi for psia in query_sets_psi]
        else:
            allpossibleCombos = [(query_sets_phi[0],query_sets_psi[0])]
        phi_psi = []
        #print('found',len(allpossibleCombos),'combos ...')
        for phi_atoms,psi_atoms in allpossibleCombos:
            psi = -777 if 'dummy' in psi_atoms.keys() else round(p3d.geo.dihedral(psi_atoms['Np1'],psi_atoms['C'],psi_atoms['CA'],psi_atoms['N'])*180/pi,prescision)
            phi = -777 if 'dummy' in phi_atoms.keys() else round(p3d.geo.dihedral(phi_atoms['C'],phi_atoms['CA'],phi_atoms['N'],phi_atoms['CB4'])*180/pi,prescision)
            if (phi,psi) not in phi_psi:
                phi_psi.append((phi,psi))
        return phi_psi

    def isHelical(self,mode='both',allowAlternativeConfs=False):
        '''
        Returns True of False if residue is helical by calling self.calcPhiPsi()
        Helical boundaries are:
            phi = range(-82,42)
            psi = range(-60,21)

        Optinal mode='both|phi|psi' to have both or only one dihedral checked 
        '''
        #if self.aa not in p3d.library.AA3:
        #   raise NotAminoAcidError()
        ''' Helical Boundaries '''
        helical_phi = range(-82,-42)
        helical_psi = range(-60,-21)
        phi,psi =  self.calcPhiPsi(allowAlternativeConfs=allowAlternativeConfs)[0]
        if mode=='phi':
            if int(phi) in helical_phi:
                return True
            else:
                return False
        elif mode=='psi':
            if int(psi) in helical_psi:
                return True
            else:
                return False
        else:
            if int(phi) in helical_phi and int(psi) in helical_psi:
                return True
            else:
                return False

    def minDistancetoResidue(atom,ResidueAtom):
        ''' 
        Determines minimum distance from a given atom to a given residue (input is any atom of that residue)
        '''
        residue_atoms = ResidueAtom.protein.hash['resid'][ResidueAtom.resid] & ResidueAtom.protein.hash['chain'][ResidueAtom.chain] & ResidueAtom.protein.hash['aa'][ResidueAtom.aa]
        min_distance = None
        closest_atom = None
        for res_atom in residue_atoms:
            d = atom.distanceTo(res_atom)
            if min_distance == None or d < min_distance:
                min_distance = d
                closest_atom = res_atom
        return res_atom,min_distance    

    def minDistancetoChain(atom,ResidueAtom):
        ''' 
        Determines minimum distance from a given atom to a given chain (input is any atom of that residue)
        '''
        chain_atoms = ResidueAtom.protein.hash['chain'][ResidueAtom.chain]
        min_distance = None
        closest_atom = None
        for res_atom in chain_atoms:
            d = atom.eval_distance(res_atom,min_distance)
            if d:
                if min_distance == None or d < min_distance:
                    min_distance = d
                    closest_atom = res_atom
        return res_atom,min_distance    

    def allAtomsOfSameResidue(self):
        return self.protein.hash['resid'][self.resid] & self.protein.hash['chain'][self.chain] & self.protein.hash['model'][self.model] & (self.protein.hash['aa-resname'][self.aa] | self.protein.hash['non-aa-resname'][self.aa])

    def allAtomsOfSameChain(self):
        return self.protein.hash['chain'][self.chain]


if __name__ == '__main__':
    print('yes')
