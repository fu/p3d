#!/usr/bin/env python
# encoding: utf-8
# Copyright (c) 2009 Ch. Fufezan.
# module is distributed under the terms of the GNU General Public License

'''
p3d

a pdb module for python by Christian Fufezan & Michael Specht
see http://p3d.fufezan.net for documentation and the latest versions
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


import os, time, re, sys
from collections import defaultdict as ddict

from p3d.__version__ import __version__ as version
import p3d.vector
import p3d.library
import p3d.tree
import p3d.atom
import p3d.parser


hash_list = ['chain','atype','resid','resname','non-aa-resname',\
             'aa-resname','model','bkb','oxygen','altconf',
             'nitrogen','non-protein','alpha','protein','altConf']

class Protein:
    '''
    p3d.Protein object

    usage p3d.protein.Protein(pdbfile,mode='3D',chains=None,MaxAtomsPerLeaf=24)

    mode:   3D creates Atom objects with x y z coordinates and initializes a BSP Tree
            2D creates only sequences (not implemented yet)

    chains: reads only those peptide chains that are in list 

    MaxAtomsPerLeaf:    How many atoms per leaf on the BSP Tree. Smaller numbers increase
                        inital tree building time but decrease then overall query times.
    '''
    def __init__(self,pdbfile,mode='3D',chains=None,residues=None,models=None,MaxAtomsPerLeaf=96,DunbrackNaming=False,BSPTree=True):
        if not os.path.exists(pdbfile): 
            raise InputError('File does not exist!')
            sys.exit(1)
        else:
            startTime = time.time()
            self.filename = str(pdbfile)
            self.hash = {}
            self.stats = ddict(int)
            self.geoCentre = p3d.vector.Vector()
            self.massCentre = p3d.vector.Vector()
            self.atoms = [] # storage for all Atom objects
            self.chainTermini = ddict(list)
            self.init_hashes()
            self.resolution = 0.01
            self.leftOvers = []
            self.atomInfoIndex = None
            self.header = []
            self.headers_infos = ddict(list)
            self.conect = []
            # --*-- remove path info from pdb filename --*--
            self.fullname = str(pdbfile) if len(pdbfile.split('/')) == 0 else str(pdbfile.split('/')[-1])
            # --*-- head of family flag used for high throuput to distinguish redundant and non-redundant set --*--
            self.HOF = False
            if len(self.fullname.split('.')[0]) <= 4:
                self.id = self.fullname[:4]
                self.dunbrackChain = None
            elif len(self.fullname.split('.')[0]) == 5 and DunbrackNaming == True:
                # special identifiers in file name, i.e. Dunbrack chain ids
                self.id = self.fullname[:4]
                self.dunbrackChain = self.fullname[4:5].upper()
                if chains == None: chains = list(self.dunbrackChain)
            elif len(self.fullname.split('.')[0]) == 6 and DunbrackNaming == True:
                # special identifiers in filename, i.e. Dunbrack chain id and head of family hook
                if self.fullname[:1] == '_':
                    self.HOF = True
                    self.id = self.fullname[1:5]
                    self.dunbrackChain = self.fullname[5:6]
            else:
                self.id = self.fullname.split('.')[0]
                self.dunbrackChain = None
            if mode == '3D':
                self.BSPTree = p3d.tree.Tree(protein=self)
                self.read_in_3Dstructure(pdbfile,chains=chains,residues=residues,models=models)
                if BSPTree:
                    self.BSPTree.build(MaxAtomsPerLeaf=MaxAtomsPerLeaf)
            elif mode == '2D':
                self.read_in_2Dstructure(pdbfile,chains=chains)
            else:
                print("mode error")
                sys.exit(1)
            self.init_parser()
            self.time = round(time.time()-startTime,3)
        return

    def init_hashes(self):
        ''' Hashtables '''
        for hash in hash_list:
            self.hash[hash] = ddict(set)    
        for AA in p3d.library.Library['AA']:
            self.hash['aa-resname'][AA] = set()
        return

    def open_pdb(self,pdbfile):
        if pdbfile.endswith('.gz'):
            import gzip,codecs
            try:
                liste = codecs.getreader("utf-8")(gzip.open(pdbfile)).readlines()
            except IOError:
                raise IOError('File not found or corrupted')
        else:
            try:
                liste = open(pdbfile).readlines()
            except IOError:
                raise IOError('File not found or corrupted')
                
        return liste

    def read_in_3Dstructure(self,pdbfile,chains=None,residues=None,models=None):
        '''
        We will overwrite any existing structure
        '''
        pdbLinePattern = """
        (?P<type>[A,H][T,E][O,T][M,A].{2})  # ATOM or HETATOM and not HETSYN
        (?P<index>[0-9, ]{5})                   # index number 5 digits
        (?P<atype>.{5})                         # Atom type e.g. CA 
        (?P<altconf>.)                          # altConf on res lvl
        (?P<aa>[A-Z,0-9, ]{3})                  # aa identifier
        ([ ]{1})                                # -Space-
        (?P<chain>.{1})                         # chain identifier
        (?P<resid>[0-9, -]{4})                  # resid
        (?P<altconf2>.)                         # altconf on chain lvl
        (?P<x>[0-9, ,\-,\.]{11})                # x coordinates
        (?P<y>[0-9, ,\-,\.]{8})                 # y coordinates
        (?P<z>[0-9, ,\-,\.]{8})                 # z coordinates
        (?P<user>[0-9, ,\-,\.]{6})              # user field
        (?P<beta>[0-9, ,\-,\.]{6})              # beta field
        ([ ]{10})                               # -Space-
        (?P<elementType>[A-Z,0-9, ]{2})         # Element Type
        (?P<charge>[A-Z,0-9,\-,\+, ]*)                # Charge
        """
        '''
        NEW and not identified :
        1-----2----3----45--6
        ----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
        ATOM   1390 HO3'  DG F  15       2.055  78.006 -19.523  1.00  0.00           H
        ATOM      1  O5'  DG C  16      -5.319  62.940 -23.530  1.00 81.62           O
        HETATM46920  C25 PQ9  5054      31.412 -13.685  -9.050  0.50 52.52           C
        '''
        """ 
        COLUMNS        DATA  TYPE    FIELD        DEFINITION
        -------------------------------------------------------------------------------------
         1 -  6        Record name   "ATOM  "
         7 - 11        Integer       serial       Atom  serial number.
        13 - 16        Atom          name         Atom name.
        17             Character     altLoc       Alternate location indicator.
        18 - 20        Residue name  resName      Residue name.
        22             Character     chainID      Chain identifier.
        23 - 26        Integer       resSeq       Residue sequence number.
        27             AChar         iCode        Code for insertion of residues.
        31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
        39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
        47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
        55 - 60        Real(6.2)     occupancy    Occupancy.
        61 - 66        Real(6.2)     tempFactor   Temperature  factor.
        77 - 78        LString(2)    element      Element symbol, right-justified.
        79 - 80        LString(2)    charge       Charge  on the atom.
        ... source : http://www.wwpdb.org/documentation/format32/sect9.html
        """
        mask = re.compile(pdbLinePattern, re.VERBOSE)

        liste = self.open_pdb(pdbfile)
        self.stats['model'] = 1
        current_model = 1
        i = 0
        for line in liste:
            ''' --*-- STRUCTURE --*-- '''
            match = mask.match(line)
            
            if "HETSYN" in line:
                # I had ValueError from 'self.idx = ...' in Atom.__init__.
                # Offending lines contained: "HETSYN     NAG NAG"
                # It works if HETSYN lines are ignored.
                continue
            
            #if match == None and line.startswith('ATOM'):
            #   print line
            if match != None:
                ''' read into structure & hash will be done on the fly '''
                #
                # ze Atom :
                CurrentAtom = p3d.atom.Atom(line=None,protein=self,PositionInAtomslist=i,model=current_model,matchObject=match)
                #
                # collect Termini resids 
                if CurrentAtom.chain not in self.chainTermini.keys():
                    self.chainTermini[CurrentAtom.chain] = [CurrentAtom.resid]
                    if i != 0:
                        self.chainTermini[self.atoms[-1].chain].append(self.atoms[-1].resid)
                # check if we should concentrate on certain chains
                if chains != None:
                    if type(chains) in [tuple,list]:
                        if CurrentAtom.chain not in chains:
                            continue
                    else:
                        raise TypeError()
                if residues != None:
                    if type(residues) in [tuple,list]:
                        if CurrentAtom.resid not in residues:
                            continue
                    else:
                        raise TypeError()
                if models != None:
                    if type(models) in [tuple,list]:
                        if CurrentAtom.model not in models:
                            continue
                    else:
                        raise TypeError()
                # --- add to atom-collection list ---
                self.atoms.append(CurrentAtom)
                # --- add to atype hash ---
                self.hash['atype'][CurrentAtom.atype.strip()].add(CurrentAtom)
                # --- add to chain hash ---
                self.hash['chain'][CurrentAtom.chain].add(CurrentAtom)
                # --- add to resid hash ---
                self.hash['resid'][CurrentAtom.resid].add(CurrentAtom)
                # --- add to altConf hash ---
                self.hash['altconf'][CurrentAtom.altConf].add(CurrentAtom)
                # --- add to all resname hash ---
                self.hash['resname'][CurrentAtom.aa.strip()].add(CurrentAtom)
                # --- add to altConf hash ---
                self.hash['altConf'][CurrentAtom.altConf].add(CurrentAtom)
                # --- add to aminoacid hash ---
                if CurrentAtom.aa in p3d.library.AA3:
                    self.hash['aa-resname'][CurrentAtom.aa].add(CurrentAtom)
                    self.hash['protein'][''].add(CurrentAtom)
                    if CurrentAtom.atype.strip() in set(['CA','C','N','O']):
                        self.hash['bkb'][''].add(CurrentAtom)
                        if CurrentAtom.atype.strip() == 'CA':
                            self.hash['alpha'][''].add(CurrentAtom)
                # --- add to non-protein hashes
                else:
                    self.hash['non-aa-resname'][CurrentAtom.aa].add(CurrentAtom)
                    self.hash['non-protein'][''].add(CurrentAtom)
                if 'O' in CurrentAtom.atype:
                    self.hash['oxygen'][''].add(CurrentAtom)
                elif 'N' in CurrentAtom.atype:
                    self.hash['nitrogen'][''].add(CurrentAtom)
                else:
                    pass
                # --- add to model hash ---
                self.hash['model'][current_model].add(CurrentAtom)
                # --- add to BSPtree sets, given index position in protein.atoms[] 
                self.BSPTree.AddAtom(CurrentAtom,i)
                i += 1
            elif (line[:21] == 'REMARK   2 RESOLUTION'):
                chopped = line.split()
                try:
                    resolution = float(chopped[3])
                except:
                    resolution = 0.01
                if resolution == 'NOT':
                    resolution = 0.00
                self.resolution = resolution
            elif (line[:5] == 'MODEL'):
                self.stats['model'] += 1
                current_model = int(line.split()[1])
            else:
                if self.atomInfoIndex == None and line.startswith('MASTER'):
                    self.atomInfoIndex = len(self.leftOvers)
                    self.header = self.leftOvers[:self.atomInfoIndex]
                    temp_list = []
                    for h in self.header:
                        header_split = h.split()
                        temp_list.append((header_split[0].lower(), h[len(header_split[0]):].strip()))
                    for k, v in temp_list:
                        self.headers_infos[k].append(v)
                if line.startswith('CONECT'):
                    tmp = line.strip()[6:]
                    tmp = [int(tmp[x:x+5]) for x in range(0, len(tmp), 5)]
                    self.conect.append((tmp))
                self.leftOvers.append('{0: <76}\n'.format(line.strip()))
                if line.startswith('ATOM'):
                    print('The following atom was missed by the reg-ex.')
                    print(line.strip())
                    print('This branch is for debugging purpose ..., ')

                    k = """
                    (?P<type>[A,H][T,E].{4})            # ATOM or HETATOM
                    (?P<index>[0-9, ]{5})               # index number 5 digits
                    (?P<atype>.{5})         # Atom type e.g. CA 
                    (?P<altconf>.)                      # altConf on res lvl
                    (?P<aa>[A-Z]{3})                    # aa identifier
                    ([ ]{1})                            # -Space-
                    """
                    '''
                    NEW and not identified :
                    1-----2----3----45--6
                    ----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
                    ATOM   1390 HO3'  DG F  15       2.055  78.006 -19.523  1.00  0.00           H
                    ATOM      1  O5'  DG C  16      -5.319  62.940 -23.530  1.00 81.62           O
                    '''
                    """
                    (?P<type>[A,H][T,E].{4})            # ATOM or HETATOM
                    (?P<index>[0-9, ]{5})               # index number 5 digits
                    (?P<atype>.{5})         # Atom type e.g. CA 
                    (?P<altconf>.)                      # altConf on res lvl
                    (?P<aa>[A-Z]{3})                    # aa identifier
                    ([ ]{1})                            # -Space-
                    (?P<chain>.{1})                     # chain identifier
                    (?P<resid>[0-9, -]{4})              # resid
                    (?P<altconf2>.)                     # altconf on chain lvl
                    (?P<x>[0-9, ,\-,\.]{11})            # x coordinates
                    (?P<y>[0-9, ,\-,\.]{8})             # y coordinates
                    (?P<z>[0-9, ,\-,\.]{8})             # z coordinates
                    (?P<user>[0-9, ,\-,\.]{6})          # user field
                    (?P<beta>[0-9, ,\-,\.]{6})          # beta field
                    """
                    mask = re.compile(k, re.VERBOSE)
                    match = mask.match(line)
                    #if match == None and line.startswith('ATOM'):
                    #   print line
                    if match != None:
                        print('adapted pattern matches line now ...')
                    else:
                        print('adapted pattern still misses line...')
                    sys.exit(1)
        self.chainTermini[self.atoms[-1].chain].append(self.atoms[-1].resid)
        return

    def init_parser(self):
        info = {}

        info['aliases'] = {
            'backbone': 'bkb',
            'residue id': 'resid',
            'residue name': 'resname',
            'atom type': 'atype',
        }

        info['functions'] = {}
        info['functions']['within {radius: float} of {centre: set}'] = self.collectSphereAtoms
        info['functions']['first residue of chain {chain: value of chain}'] = self.firstResidueOfChain
        info['functions']['last residue of chain {chain: value of chain}'] = self.lastResidueOfChain

        info['caseSensitive'] = set(['chain'])

        self.parser = p3d.parser.Parser(self.hash, info)
        return

    def query(self, *args):
        return list(self.parser.parse(*args))

    def lookUpAtom(self, *args):
        a = list(self.parser.parse(*args))
        #if len(a) > 1:
        #   print('lookUpAtom found more than one atom, will return first in list!')
        if len(a) == 0:
            return False
        elif len(a) == 0:
            return False
        return a[0]


    def info(self):
        '''
        Gives information about the protein, the hash tables and the BSPTree
           - 
        '''
        output = []
        output.append('PDBFile infos')
        output.append('  '+self.fullname+' ( '+self.id+' / '+str(self.dunbrackChain)+' )')
        output.append('    '+str(self.time)+' s to parse pdb file and create BSP Tree')
        output.append( '\n--*-- Hashtable occupancies --*--')
        for key,table in self.hash.items():
            output.append('  '+key+':'+'( '+str(len(table.keys()))+' items)')
            #output.append("".join(table.keys()))
            lineoutput = ''
            for table_key,item in sorted(self.hash[key].items()):
                if len(lineoutput) < 80:
                    lineoutput += '{0: >10} : {1: <6}'.format(table_key,len(item))
                else:
                    output.append(lineoutput)
                    lineoutput = '{0: >10} : {1: <6}'.format(table_key,len(item))
            output.append(lineoutput+'\n')
        output.append('\n--*-- BSP Tree --*--')
        output += self.BSPTree.info()
        return output

    def min_Distance_to_Residue(self, centre=None,ResidueAtom=None):
        '''
        Determines minimum distance from a given atom to a given residue (input is any atom of that residue)
        '''
        if AtomSet == None:
            raise AtLeastOneVectorIsNeeded()
        residue_atoms = self.hash['resid'][ResidueAtom.resid] & self.hash['chain'][ResidueAtom.chain] & self.hash['aa'][ResidueAtom.aa]
        return self.min_Distance_to_Set(centre,AtomSet=ResidueAtom)

    def min_Distance_to_Set(self, centre=None,AtomSet=None):
        '''
        Determines minimum distance from a given atom to a given set of atoms
        '''
        if AtomSet == None:
            raise AtLeastOneVectorIsNeeded()
        min_distance = 777
        for Setatom in AtomSet:
            d = Setatom.distanceTo(centre)
            #print(d,Setatom.info(),centre.info(),'<<<<')
            if d < min_distance:
                min_distance = d
        return min_distance

    def distanceBetweenTwoSets(self, setA, setB):
        min_distance = 777
        for atom in setA:
            d = self.min_Distance_to_Set(AtomSet=setB)
            if d < min_distance:
                min_distance = d
        return min_distance

    def collectSphereAtoms(self,centre=None,radius=0):
        if centre != None:
            collection = set([])
        else:
            raise AtLeastOneVectorIsNeeded()
        for atom in centre:
            collection |= set(self.BSPTree.query(atom,radius=radius))
        return collection

    def writeToFile(self,filename,includeOrgHeader=False):
        HEADER = [\
            'HEADER p3d-modified'+self.filename,\
            'REMARK p3d  V.'+str(version),\
            'REMARK ---*--- please cite Fufezan & Specht 2009 BMC Bioinformatics 10:258 ---*---',\
            'REMARK'\
        ]
        dump = open(filename, 'w')
        # wrting header
        for headline in HEADER:
            dump.write(headline+'\n')
        if includeOrgHeader == True:
            for entry in self.header:
                dump.write(entry)
        # --*-- NMR Structure ? --*--
        if self.stats['model'] > 1:
            current_model = None
            for atom in self.atoms:
                if atom.model != current_model:
                    if current_model != None:
                        dump.write('ENDMDL')
                    current_model = atom.model
                    dump.write('MODEL        '+str(current_model))
                dump.write(atom.output()+'\n')
        else:
            for atom in self.atoms:
                dump.write(atom.output()+'\n')

        dump.write('END')
        dump.close()
        return

    def generatePDBoutput(self):
        HEADER = [\
            'HEADER p3d-modified'+self.filename,\
            'REMARK p3d  V.'+str(version),\
            'REMARK ---*--- please cite Fufezan & Specht 2009 BMC Bioinformatics 10:258 ---*---',\
            'REMARK'\
        ]
        HEADER += self.leftOvers[:self.atomInfoIndex]
        for atom in self.atoms:
            HEADER.append('{0: <76}\n'.format(atom.output()))
        return HEADER+self.leftOvers[self.atomInfoIndex:]

    def firstResidueOfChain(self,chain,idOnly=False):
        if chain not in self.chainTermini.keys():
            print('>>>>',chain)
            raise UnknownChain()
        else:
            resid = self.chainTermini[chain][0]
            if idOnly: 
                return resid
            else:
                return self.hash['resid'][resid] & self.hash['chain'][chain]

    def lastResidueOfChain(self,chain,idOnly=False):
        if chain not in self.chainTermini.keys():
            raise UnknownChain()
        else:
            resid = self.chainTermini[chain][1]
            if idOnly:
                return resid
            else:
                return self.hash['resid'][resid] & self.hash['chain'][chain]


class TypeError(Exception): pass
class InputError(Exception): pass
class AtLeastOneVectorIsNeeded(InputError): pass
class QueryResultedInMoreThanOneAtom(InputError): pass
class UnknownChain(InputError): pass

if __name__ == "__main__":
    print("yes")
