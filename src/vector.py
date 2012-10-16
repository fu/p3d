#!/usr/bin/env python3
# encoding: utf-8
# Copyright (c) 2009 Ch. Fufezan.
# module is distributed under the terms of the GNU General Public License

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


import math, random
from copy import deepcopy as dcp

class KeyError(Exception): pass

class Vector:
    def __init__(self,x=0,y=0,z=0):
        '''init as vector.Vector(x,y,z)'''
        __slots__ = ["x","y","z"]
        self.x = x
        self.y = y
        self.z = z
        self.xtra = []
        self.idx   = -1
        self.atype = 'uNk'
        self.altConf = ' '
        self.aa = 'uNk'
        self.chain = ' '
        self.resid = '-1'
        self.altConf2 = ' '
        self.uberID = ' '
        self.desc = ''
        self.type = 'p3dVec'
        self.user = 77.7
        self.beta = 77.7
        self.protein = None
        self.pos_in_list = -1
        self.__bases__ = 'Vector'
        self.elementType = '  '
        self.charge = None
        return

    def __add__(self,zeOther):
        '''
        Returns new Vector = a + b
        Vector.desc contains history of operation with index position in protein.atoms
        '''
        tmp = Vector(self.x+zeOther.x,self.y+zeOther.y,self.z+zeOther.z)
        tmp.desc = 'idx:'+str(self.pos_in_list)+'+idx:'+str(zeOther.pos_in_list)
        return tmp

    def __sub__(self,zeOther):
        '''
        Returns new Vector = a - b, resulting Vector points towards self, i.e. a
        Vector.desc contains history of operation with index position in protein.atoms
        '''
        tmp = Vector(self.x-zeOther.x,self.y-zeOther.y,self.z-zeOther.z)
        tmp.protein = self.protein
        tmp.desc = 'idx:'+str(self.pos_in_list)+'-idx:'+str(zeOther.pos_in_list)
        return tmp

    def __mul__(self,value):
        '''
        Returns new Vector = self * scalar
        Vector.desc contains history of operation with index position in protein.atoms
        '''
        protein = self.protein
        self.protein = None
        v = dcp(self)
        v.protein = protein
        v.x = self.x*float(value)
        v.y = self.y*float(value)
        v.z = self.z*float(value)
        v.desc = 'idx:'+str(v.desc)+'*'+str(value)
        return v

    def __rmul__(self,value):
        '''
        Returns new Vector = scalar * self
        Vector.desc contains history of operation with index position in protein.atoms
        '''
        protein = self.protein
        self.protein = None
        v = dcp(self)
        v.protein = protein
        v.x = self.x*float(value)
        v.y = self.y*float(value)
        v.z = self.z*float(value)
        v.desc = 'idx:'+str(v.desc)+'*'+str(value)
        return v

    def __div__(self,value):
        '''
        Returns new Vector = self / scalar
        Vector.desc contains history of operation with index position in protein.atoms
        '''
        protein = self.protein
        self.protein = None
        v = dcp(self)
        v.protein = protein
        v.x = self.x/float(value)
        v.y = self.y/float(value)
        v.z = self.z/float(value)
        v.desc = 'idx:'+str(v.desc)+'/'+str(value)
        return v

    def __rdiv__(self,value):
        '''
        Returns new Vector = self / scalar
        Vector.desc contains history of operation with index position in protein.atoms
        '''
        protein = self.protein
        self.protein = None
        v = dcp(self)
        v.protein = protein
        v.x = self.x/float(value)
        v.y = self.y/float(value)
        v.z = self.z/float(value)
        v.desc = 'idx:'+str(v.desc)+'/'+str(value)
        return v

    def __truediv__(self,value):
        '''
        Returns new Vector = self / scalar
        Vector.desc contains history of operation with index position in protein.atoms
        '''
        protein = self.protein
        self.protein = None
        v = dcp(self)
        v.protein = protein
        v.x = self.x/float(value)
        v.y = self.y/float(value)
        v.z = self.z/float(value)
        v.desc = 'idx:'+str(v.desc)+'/'+str(value)
        return v

    def __abs__(self):
        '''
        Returns length of self
        '''
        return float(math.sqrt(self.x**2 + self.y**2 + self.z**2))

    def __neg__(self):
        '''
        Returns new Vector = self * -1
        '''
        protein = self.protein
        self.protein = None
        v = dcp(self)
        v.protein = protein
        return v*(-1)

    def __len__(self):
        '''
        has to be integer in python, thus doesn't make much sense ...'
        return float(math.sqrt(self.x**2 + self.y**2 + self.z**2))

        for length use self.length()
        '''
        return

    def dot(self,zother):
        '''
        returns dot product of the two vectors
        '''
        return self.x*zother.x + self.y*zother.y + self.z*zother.z

    def normalize(self):
        '''
        Return new normalized vector, i.e. length of 1
        '''
        protein = self.protein
        self.protein = None
        v = dcp(self)
        v.protein = protein
        return v/abs(v)

    def translateBy(self,zother):
        '''
        This differs from substraction as it preserves all properties of self
        '''
        self.x = self.x+zother.x
        self.y = self.y+zother.y
        self.z = self.z+zother.z
        return self

    def rotate(self,p1,p2,phi,angleunit='degree'):
        """
        Rotates Vector by angle phi around axis defined by P1 -> P2
        Angle in degrees
        taken from : http://www.mines.edu/~gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.html
        next implementation could use internal funktions and http://en.wikipedia.org/wiki/Rodrigues'_rotation_formula

        !!! NOT FULLY TESTED !!!
        """
        if angleunit == 'degree':
            phiInrad = math.radians(phi)
        else:
            phiInrad = phi
        cosphi = math.cos(phiInrad)
        sinphi = math.sin(phiInrad)
        rotate_axis = p2-p1
        x = self.x
        y = self.y
        z = self.z
        a = p1.x
        b = p1.y
        c = p1.z
        u = rotate_axis.x
        v = rotate_axis.y
        w = rotate_axis.z
        u2 = u**2
        v2 = v**2
        w2 = w**2
        total2 = u2+v2+w2
        sumsqrt = math.sqrt(total2)
        ##>> Some abreviations
        rotated_x = (a*(v2+w2) + u*(-b*v -c*w +u*x +v*y +w*z) + ((x-a)*(v2+w2) +u*(b*v +c*w -v*y -w*z)) * cosphi + sumsqrt*( b*w -c*v -w*y +v*z)*sinphi) / total2
        rotated_y = (b*(u2+w2) + v*(-a*u -c*w +u*x +v*y +w*z) + ((y-b)*(u2+w2) +v*(a*u +c*w -u*x -w*z)) * cosphi + sumsqrt*(-a*w +c*u +w*x -u*z)*sinphi) / total2
        rotated_z = (c*(u2+v2) + w*(-a*u -b*v +u*x +v*y +w*z) + ((z-c)*(u2+v2) +w*(a*u +b*v -u*x -v*y)) * cosphi + sumsqrt*( a*v -b*u -v*x +u*y)*sinphi) / total2
        self.x = round(rotated_x,3)
        self.y = round(rotated_y,3)
        self.z = round(rotated_z,3)
        return self

    def rotate2(self,p1,p2,phi,angleunit='degree'):
        """
        Rotates Vector by angle phi around axis defined by P1 -> P2
        Angle in degrees
        taken from http://en.wikipedia.org/wiki/Rodrigues'_rotation_formula

        !!! NOT FULLY TESTED !!!
        """
        if angleunit == 'degree':
            phiInrad = math.radians(phi)
        else:
            phiInrad = phi
        cosphi = math.cos(phiInrad)
        sinphi = math.sin(phiInrad)
        rotateAxis = (p2-p1).normalize()
        self.translateBy(-p1)
        crossed = rotateAxis.cross(self)
        dotted = rotateAxis.dot(self)
        k = (self*cosphi + crossed*sinphi + dotted*(1-cosphi)*rotateAxis)
        self.x = round(k.x,5)
        self.y = round(k.y,5)
        self.z = round(k.z,5)
        self.translateBy(p1)
        return

    def rotate_aroundX(self,phi,angleunit='degree'):
        if angleunit == 'degree':
            phiInrad = math.radians(phi)
        else:
            phiInrad = phi
        cosphi = math.cos(phiInrad)
        sinphi = math.sin(phiInrad)
        new_y = round((self.y *  cosphi + self.z * sinphi),3)
        new_z = round((self.y * -sinphi + self.z * cosphi),3)
        self.y = new_y
        self.z = new_z
        return self

    def rotate_aroundY(self,phi,angleunit='degree'):
        if angleunit == 'degree':
            phiInrad = math.radians(phi)
        else:
            phiInrad = phi
        cosphi = math.cos(phiInrad)
        sinphi = math.sin(phiInrad)
        new_x = round((self.x * cosphi  - self.z * sinphi),3)
        new_z = round((self.x * sinphi  + self.z * cosphi),3)
        self.x = new_x
        self.z = new_z
        return self

    def rotate_aroundZ(self,phi,angleunit='degree'):
        if angleunit == 'degree':
            phiInrad = math.radians(phi)
        else:
            phiInrad = phi
        cosphi = math.cos(phiInrad)
        sinphi = math.sin(phiInrad)
        new_x = round((self.x *  cosphi + self.y * sinphi),3)
        new_y = round((self.x * -sinphi + self.y * cosphi),3)
        self.x = new_x
        self.y = new_y
        return self

    def cross(self,zother):
        '''
        Returns new Vector from cross product between two vectors
        Vector.desc contains history of operation with index position in protein.atoms
        '''
        v = Vector(self.y*zother.z - self.z*zother.y,self.z*zother.x-self.x*zother.z,self.x*zother.y-self.y*zother.x)
        v.desc = 'idx:'+str(self.pos_in_list)+' cross idx:'+str(zother.pos_in_list)
        return v

    def length(self):
        '''
        Return length of vector
        '''
        return float(math.sqrt(self.x**2 + self.y**2 + self.z**2))

    def distanceTo(self, zeOther):
        '''
        Returns distance between two Vectors
        '''
        return math.sqrt(((self.x-zeOther.x)**2) + ((self.y-zeOther.y)**2) + ((self.z-zeOther.z)**2))

    def angleBetween(self,zother):
        '''
        Return angle between two Vectors in RADs
        '''
        # Math domain error ?
        return math.acos(float(self.dot(zother)/(abs(self)*abs(zother))))

    def vmdOutput(self,radius=0.2):
        '''
        Returns TK Console command line that can be used to visualise vector self
        '''
        return 'graphics 0 sphere {0} {1} {2} radius {3};'.format(self.x,self.y,self.z,float(radius))


    def evalDistance(self,other,distance):
        '''
        Evaluates if two vectors are within distance and is faster than computing distance at once.
        Funktion also returns false if vector are not within distance or returns the computed distance
        if vectors are within evaluated distance.
        '''
        if distance == None:
            return self.distanceTo(other)
        else:
            dsquared = float(distance)**2
            xs_diff = (self.x-other.x)**2
            if xs_diff < dsquared:
                xsys_diff = xs_diff+(self.y-other.y)**2
                if xsys_diff < dsquared:
                    xsyszs_diff = xsys_diff + (self.z-other.z)**2
                    if xsyszs_diff < dsquared:
                        return math.sqrt(xsyszs_diff)
            return False

    def evalDistanceToCoordinates(self,x,y,z,distance):
        '''
        Evaluates distance from a Vector to a pair of coordinates
        It is faster than computing distance at once.
        Funktion also returns false if vector are not within distance
        or returns the computed distance
        if evaluated distance is within distance.
        '''
        if distance == None:
            return self.distanceTo(other)
        else:
            dsquared = float(distance)**2
            xs_diff = (self.x-x)**2
            if xs_diff < dsquared:
                xsys_diff = xs_diff+(self.y-y)**2
                if xsys_diff < dsquared:
                    xsyszs_diff = xsys_diff + (self.z-z)**2
                    if xsyszs_diff < dsquared:
                        return math.sqrt(xsyszs_diff)
            return False



    def output(self,format='pdb'):
        '''
        HETATM
         ATOM    559  CA BASP A  74      48.780  13.254  -1.818  0.50 16.34           C
         ----.----|----.----|----.----|----.----|----.----|----.----|----.----|----.----|----.----
             5    10   15   20   25   30   35   40   45   50   55   60   65   70
        '''
        if format=='pdb':
            altconf = self.altConf if self.altConf != '_' else ' '
            # return "{type: <6}{idx: >5}{atype: <5}{alt1:1}{resname: >3} {chain:1}{resid:>4}{alt2:1}   {x:> 8.3f}{y:> 8.3f}{z:> 8.3f}{user:> 6.2f}{beta:> 6.2f}          {et:2}".format(\
            # type=self.type, idx=self.idx,atype=self.atype,alt1=altconf,alt2=self.altConf2,resname=self.aa,\
            # chain=self.chain,resid=self.resid,x=self.x,y=self.y,z=self.z,user=self.user,beta=self.beta,et=self.elementType)
            return "{type: <6}{idx: >5}{atype: <5}{alt1:1}{resname: >3} {chain:1}{resid:>4}{alt2:1}   {x:>8.3f}{y:>8.3f}{z:>8.3f}{user: >6.2f}{beta: >6.2f}{et: >12}".format(\
            type=self.type, idx=self.idx,atype=self.atype,alt1=altconf,alt2=self.altConf2,resname=self.aa,\
            chain=self.chain,resid=self.resid,x=self.x,y=self.y,z=self.z,user=self.user,beta=self.beta,et=self.elementType)
        else:
            raise KeyError()
        return

    def upgradeToAtom(self,idx,atype,resname,chain,resid):
        self.type = 'ATOM'
        self.idx = idx
        self.atype = atype
        self.aa = resname
        self.chain = chain
        self.resid = resid
        return

    def info(self,lvl=None):
        '''
        Returns info string about vector/atom
        Default a string is return with atom type amino acid chain and residue id.

        key=lvl
            Options:
                max:         atom type, Amino acid, Chain, Residue id and x,y,z coodinates
                coordinates: x,y,z coodinates as list

        '''
        if lvl=='max':
            return '{atype: >10s}{alt1:1}{aa:3} {chain:1} {resid: >4}{alt2:1} [{x: >8.3f},{y: >8.3f},{z: >8.3f}]'.format(\
            x=self.x,y=self.y,z=self.z,atype=self.atype,aa=self.aa,chain=self.chain,resid=self.resid,alt1=self.altConf,\
            alt2=self.altConf2)
        elif lvl=='min':
            return '{atype: >10s} {aa:3} {chain:1} {resid: >4}'.format(\
            atype=self.atype,aa=self.aa,chain=self.chain,resid=self.resid)
        else:
            return '[{x: >8.3f},{y: >8.3f},{z: >8.3f}]'.format(x=self.x,y=self.y,z=self.z)

    def upgradeToAtom(self,idx,atype,resname,chain,resid):
        self.type = 'ATOM'
        self.idx = idx
        self.atype = atype
        self.aa = resname
        self.chain = chain
        self.resid = resid
        return




def bp(message=''):
    if message != '':
        print('<ERROR>',message)
        print('exiting ...')
        sys.exit(1)
    sys.exit(0)

def dihedral(a,b,c,d):
    p1n = PLANE(a,b,c).n
    p2n = PLANE(b,c,d).n
    angle = p1n.AngleBetween(p2n)
    angle_sign = (c-b).dot(p1n.cross(p2n))
    if angle_sign < 0:
        return -1*angle
    else:
        return angle

def strTr( text, dic ):
    pat = "(%s)" % "|".join( map(re.escape, dic.keys()) )
    return re.sub( pat, lambda m:dic[m.group()], text)

