#!/usr/bin/env python
# encoding: utf-8
# Copyright (c) 2009 Ch. Fufezan.
# module is distributed under the terms of the GNU General Public License

'''p3d.tree - Binary space partitioning tree for fast queries in 3D 
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


import sys, os, time
from p3d import vector as vector
import p3d
#import bisect
import operator
from copy import deepcopy as dcp

class Tree(dict):
    '''
    p3d Binary Space Partitioning Tree
    atoms of pdb structure are divided on tree leafs.
    This accelerats the search for neighbours.

    Information about the Tree is stored in 
    protein.BSPTree.informations
    '''
    def __init__(self,protein=None,dimensions=None,root_info=None,lookup=None):
        if protein != None:
            ## init new Tree 
            self.dimensions = {'x':[],'y':[],'z':[]}
            self.information = {}
            self.protein = protein
            self.__version__ = '2.41a'
            self.lookup = []
            self.atoms = []
            """ How many atoms on leaf ?"""
        else:
            ## init new branching of old tree
            self.dimensions = dimensions
            self.information = root_info
            self.lookup = lookup
        self.Cutted_Coords = -1
        self.Cutted_Dimension = ''
        self.planned_cut = -1
        root_info = {} if root_info == None else root_info
        self.Leaf = False
        self.lvl = 0
        self.kidz = [0 , 0]
        # smaller kid in 0, bigger in 1 ...
        return

    def info(self):
        output = []
        output.append('Tree Version '+str(self.__version__))
        for lable,info in self.information.items():
            if type(info) in [dict]:
                output.append(lable+':')
                for sub_lable,value in info.items():
                    output[-1] += "{0}: [{1:3.2f},{2:3.2f}] ".format(sub_lable,value[0],value[1])
            else:
                output.append(lable+':'+str(info))
        return output

    def AddAtom(self,atom,index):
        """
        self.dimensions['x'].append((atom.x,index))
        self.dimensions['y'].append((atom.y,index))
        self.dimensions['z'].append((atom.z,index))
        """
        self.lookup.append((index,atom.x,atom.y,atom.z))
        # ><><>< can avoid to pass object <><><>
        return

    def addDummyAtom(self,x,y,z,index):
        self.lookup.append((index,x,y,z))
        return

    def getDimensions(self):
        for k,dimension in enumerate(['x','y','z']):
            self.dimensions[dimension] = [max([e[k+1] for e in self.lookup]), min([e[k+1] for e in self.lookup])]
        return # soft cut - sort only dimension required

    def SortLists(self):
        '''Will sort dimension tuples by coordinate, descending'''
        for dimension in ['x','y','z']:
            self.dimensions[dimension] = sorted(self.dimensions[dimension], key=operator.itemgetter(0))
            self.information['maxDimensions'][dimension] = (self.dimensions[dimension][0][0],self.dimensions[dimension][-1][0])
        return

    def preFlight(self):
        ''' Routin to prepare Tree and variables before building '''
        # Geometric Centre of protein box
        #self.protein.geoCentre = vector.Vector(self.information['maxDimensions']['x'][1]-self.information['maxDimensions']['x'][0],self.information['maxDimensions']['y'][1]-self.information['maxDimensions']['y'][0],self.information['maxDimensions']['z'][1]-self.information['maxDimensions']['z'][0])
        #print self.protein.geoCentre.info()
        #print self.information['maxDimensions']['x'][0]-self.information['maxDimensions']['x'][1]
        return

    def build(self,MaxAtomsPerLeaf=96):
        StartingT = time.time()
        self.information['Tree level'] = self.lvl
        self.information['Leaf atom count'] = 0
        self.information['cuts'] = 0
        self.information['MaxAtomsPerLeaf'] =  MaxAtomsPerLeaf
        #self.information['maxDimensions'] = {'x':[],'y':[],'z':[]}
        #maxDimensions in []: 0=max 1=min
        self.getDimensions()
        self.preFlight()
        self.generate()
        self.information['Building Time'] = time.time()-StartingT
        return

    def generate(self):
        self.aim()
        if self.lvl >= 2000:
            print('Tree reached 2000th branching level, have to stop - we went nuts')
            sys.exit(1)
        if self.Leaf == False:
            self.splitup()
            self.kidz[0]= Tree(lookup=self['smaller'],root_info=self.information,dimensions=self['smallerDimensions'])
            self.kidz[1]= Tree(lookup=self['bigger'],root_info=self.information,dimensions=self['biggerDimensions'])
            self.cleanUp()
            self.kidz[0].generate()
            self.kidz[1].generate()
        else:
            self.information['Leaf atom count'] += len(self.lookup)
            if self.lvl > self.information['Tree level']:
                self.information['Tree level'] = self.lvl
        return

    def aim(self):
        ''' Deciding cut bu dimension '''
        dmax = (None,0.0,None)
        for k,dimension in enumerate(['x','y','z']):
            expansion = abs(self.dimensions[dimension][0] - self.dimensions[dimension][1])
            if dmax[1] < expansion:
                dmax = (dimension,expansion,k+1) # k+1 index position in loopuk
        if len(self.lookup) >= self.information['MaxAtomsPerLeaf']:
            # Schwartzian transform
            #self.lookup[:] = [ (x[dmax[2]] ,x) for x in self.lookup ]
            #self.lookup.sort()
            #self.lookup[:] = [ k[1] for k in self.lookup ]
            self.lookup = sorted(self.lookup, key=operator.itemgetter(dmax[2]))
            #self.lookup.sort(key=operator.itemgetter(dmax[2]),reverse=True)

            planned_cut = int((len(self.lookup)+1)/2)
            while self.lookup[planned_cut][dmax[2]] == self.lookup[planned_cut-1][dmax[2]]:
                planned_cut -= 1
            self.Cutted_Dimension = dmax[0]
            self.planned_cut = planned_cut
            self.Cutted_Coords = float((self.lookup[planned_cut][dmax[2]]+self.lookup[planned_cut-1][dmax[2]])/2.00000)
        else:
            self.Leaf = True
            #sum(map(operator.itemgetter(0),a))
        return

    def splitup(self):
        if self.Leaf == False:
            self.information['cuts'] += 1
            self['smaller']  = []
            self['bigger']   = []
            self['smaller']= self.lookup[:self.planned_cut]
            self['bigger'] = self.lookup[self.planned_cut:]
            kid_dimensions = filter(lambda x: x != self.Cutted_Dimension, self.dimensions.keys())
            self['smallerDimensions'] = {self.Cutted_Dimension:[self.Cutted_Coords,self.dimensions[self.Cutted_Dimension][1]]}
            self['biggerDimensions'] = {self.Cutted_Dimension:[self.dimensions[self.Cutted_Dimension][0],self.Cutted_Coords]}
            for k in kid_dimensions:
                self['smallerDimensions'][k] = self.dimensions[k]
                self['biggerDimensions'][k] = self.dimensions[k]
        return

    def cleanUp(self):
        '''Cleaning node of elements so less memory is used'''
        #print 'Cutted lvl',self.lvl,'in',self.Cutted_Dimension,'@',self.Cutted_Coords

        #print '    #elements:',len(self.dimensions[self.Cutted_Dimension]),self.dimensions[self.Cutted_Dimension][0],self.dimensions[self.Cutted_Dimension][-1]
        self.kidz[0].lvl = self.lvl + 1
        self.kidz[1].lvl = self.lvl + 1
        del self['bigger']
        del self['smaller']
        if self.lvl != 0:
            del self.lookup
        del self.dimensions
        return

    def query(self,Vector_a=vector.Vector(),Vector_b=vector.Vector(),radius=0,returnIndices=False):
        '''
        Tree.query(Vector a, vector b, radius=(in A))
        ------------
        	optional: 	Small-VECTOR
        				radius in A
        '''
        indcs = []
        Vector_radius = vector.Vector(radius,radius,radius)
        if radius != 0:
            # we have a sphere and evaluate distance later
            Vector_l = Vector_a - Vector_radius
            Vector_r = Vector_a + Vector_radius
            bbox = {'x':[Vector_l.x,Vector_r.x],'y':[Vector_l.y,Vector_r.y],'z':[Vector_l.z,Vector_r.z]}
        elif abs(Vector_a) != 0 and abs(Vector_b) != 0:
            # we have got a box
            bbox = {'x':[Vector_b.x,Vector_a.x],'y':[Vector_b.y,Vector_a.y],'z':[Vector_b.z,Vector_a.z]}
        else:
            print(self.query.__doc__)
            sys.exit(1)
        Nestedindcs = self.walk(bbox,indcs)
        flattend = flattenNested(Nestedindcs)
        atoms = []
        if returnIndices == False:
            for i in list(flattend):
                if radius != 0:
                    #sphere collection
                    distance = Vector_a.evalDistance(self.protein.atoms[i],radius)
                    if 0 < distance <= radius:
                        atoms.append(self.protein.atoms[i])
                else:
                    atoms.append(self.protein.atoms[i])
            return atoms
        else:
            return flattend

    def walk(self,bbox,indcs):
        #print bbox,indcs
        # {'y': [96.802999999999997, 106.803], 'x': [46.704000000000001, 56.704000000000001], 'z': [26.271000000000001, 36.271000000000001]} []

        if self.Leaf == False:
            #print '\nNew Walk lvl:',self.lvl,'Cutted around',self.Cutted_Dimension,'at:',self.Cutted_Coords
            #print '<NOTE> -   Bigcoo:','at:',self.Cutted_Dimension,bbox[self.Cutted_Dimension][1]
            #print '<NOTE> - Smallcoo:','at:',self.Cutted_Dimension,bbox[self.Cutted_Dimension][0]
            cutted = float(self.Cutted_Coords)
            #print type(self.Cutted_Coords)
            if bbox[self.Cutted_Dimension][1] < cutted:
                # BBOX completly in Smaller half
                #print self.lvl,bbox[self.Cutted_Dimension][1],'<',cutted,' aka BBOX in smaller half, axis:',self.Cutted_Dimension
                indcs = self.kidz[0].walk(bbox,indcs)
                #newindcs = self.kidz[0].walk(bbox,indcs)
            elif bbox[self.Cutted_Dimension][0] > cutted:
                # BBOX Completely in Bigger half
                #print self.kidz[1].Cutted_Dimension
                #print self.lvl,bbox[self.Cutted_Dimension][0],'>',cutted,' aka BBOX in bigger half, axis:',self.Cutted_Dimension
                indcs = self.kidz[1].walk(bbox,indcs)
            else:
                #print 'init 2 walks'
                #print self.kidz[0].Cutted_Dimension
                #print self.kidz[1].Cutted_Dimension
                #print self.lvl,bbox[self.Cutted_Dimension][0],'<>',cutted,' aka BBOX in both halfs, axis:',self.Cutted_Dimension,' - diving in both branches'
                indcs = self.kidz[0].walk(bbox,indcs)
                indcs = self.kidz[1].walk(bbox,indcs)

        else:
            #print '\nFound leaf @ level',self.lvl
            #indcs.append(list(map(operator.itemgetter(1),self.dimensions['x'])))
            indcs += [index for (index,x,y,z) in self.lookup]
            #print 'element count @',indcs
            #print 'x-items:',sorted(map(operator.itemgetter(1),self.dimensions['x']))
            #print 'y-items:',sorted(map(operator.itemgetter(1),self.dimensions['y']))
            #print 'z-items:',sorted(map(operator.itemgetter(1),self.dimensions['z']))
            #print bbox
            #sys.exit(1)
            #return indcs
        #exit(1)
        return indcs

    def generateSurface(self,ListOfVectors,MinDistance,MaxDistance,MinDistanceOfSurfaceVectors):
        assert isinstance (ListOfVectors[0],p3d.vector.Vector), "Cannot generate Surface, require a list of p3d vector(s) or atoms"
        GridDimensions = {}

        SurfaceVectors = []
        SurfaceOuter = set()
        SurfaceInner = set()

        ''' Collecting all atoms within Max and Min distance '''
        for j,atom in enumerate(ListOfVectors):
            positions = self.query(atom,radius=MaxDistance,returnIndices=True)
            for pos in positions:
                index,x,y,z = self.lookup[pos]
                assert index==pos, 'Position missmatch'
                d = atom.evalDistanceToCoordinates(x,y,z,MaxDistance)
                if d:
                    if d >= MinDistance:
                        SurfaceOuter.add((x,y,z))
                    else:
                        SurfaceInner.add((x,y,z))
        ''' Generating Surface with distance between Gripoints at list MinDistanceOfSurfaceVectors '''
        notInteresting = set()
        print('REMARK outer {0}, inner {1}, outer-inner {2}'.format(len(SurfaceOuter),len(SurfaceInner),len(SurfaceOuter-SurfaceInner)))
        for k,(x,y,z) in enumerate(SurfaceOuter-SurfaceInner):
            t = p3d.vector.Vector(x,y,z)
            if (x,y,z) not in notInteresting: 	
                positions = self.query(t,radius=MinDistanceOfSurfaceVectors,returnIndices=True)
                SurfaceVectors.append(t)
                for pos in positions:
                    index,x,y,z = self.lookup[pos]
                    if t.evalDistanceToCoordinates(x,y,z,MinDistanceOfSurfaceVectors):
                        notInteresting.add((x,y,z))
        return SurfaceVectors	

""" END OF CLASS """

def flattenNested(liste):
    for element in liste:
        if type(element) in [ tuple, list]:
            for nested in flattenNested(element):
                yield nested
        else:
            yield element

def bp(message=''):
    if message != '':
        print('<ERROR>',message)
        print('exiting ...')
        sys.exit(1)
    sys.exit(0)
