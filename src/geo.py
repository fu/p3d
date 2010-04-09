#!/usr/bin/env python
# encoding: utf-8
# Copyright (c) 2009 Ch. Fufezan.
# module is distributed under the terms of the GNU General Public License


"""geo part of p3d
see. p3d.fufezan.net for more detail
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
import time
from copy import deepcopy as dcp
import p3d.vector

class TransformationMatrix:
	def __init__(self,isv,itv):
		'''
		Returns a transformationmatrix object
		input has to be a set of 3 source vectors and 3 destination vectors  
		'''
		assert len(isv) == 3, "Need three source vectors"
		assert len(itv) == 3, "Need three target vectors"
		self.a = 0
		self.b = 0
		self.c = 0
		self.d = 0
		self.e = 0
		self.f = 0
		self.g = 0
		self.h = 0
		self.i = 0
		sb1 = (isv[1]-isv[0])
		sb2 = (isv[2]-isv[0])
		su3 = (sb1.cross(sb2)).normalize()
		su2 = sb2.normalize()
		su1 = (su2.cross(su3)).normalize()
		
		tb1 = (itv[1]-itv[0])
		tb2 = (itv[2]-itv[0])
		tu3 = (tb1.cross(tb2)).normalize()
		tu2 = tb2.normalize()
		tu1 = (tu2.cross(tu3)).normalize()
		
		self.source = isv[0]
		self.target = itv[0]
		sourcevectors = [ su1, su2, su3 ]
		targetvectors = [ tu1, tu2, tu3 ]
		
		s1x = sourcevectors[0].x
		s2x = sourcevectors[1].x
		s3x = sourcevectors[2].x
		s1y = sourcevectors[0].y
		s2y = sourcevectors[1].y
		s3y = sourcevectors[2].y
		s1z = sourcevectors[0].z
		s2z = sourcevectors[1].z
		s3z = sourcevectors[2].z
		t1x = targetvectors[0].x
		t2x = targetvectors[1].x
		t3x = targetvectors[2].x
		t1y = targetvectors[0].y
		t2y = targetvectors[1].y
		t3y = targetvectors[2].y
		t1z = targetvectors[0].z
		t2z = targetvectors[1].z
		t3z = targetvectors[2].z
		self.a = -(s2z*s3y*t1x - s2y*s3z*t1x - s1z*s3y*t2x + s1y*s3z*t2x + s1z*s2y*t3x - s1y*s2z*t3x)/(-s1z*s2y*s3x + s1y*s2z*s3x + s1z*s2x*s3y - s1x*s2z*s3y - s1y*s2x*s3z + s1x*s2y*s3z)
		self.b = -(s2z*s3x*t1x - s2x*s3z*t1x - s1z*s3x*t2x + s1x*s3z*t2x + s1z*s2x*t3x - s1x*s2z*t3x)/(s1z*s2y*s3x - s1y*s2z*s3x - s1z*s2x*s3y + s1x*s2z*s3y + s1y*s2x*s3z - s1x*s2y*s3z)
		self.c = -(-s2y*s3x*t1x + s2x*s3y*t1x + s1y*s3x*t2x - s1x*s3y*t2x - s1y*s2x*t3x + s1x*s2y*t3x)/(s1z*s2y*s3x - s1y*s2z*s3x - s1z*s2x*s3y + s1x*s2z*s3y + s1y*s2x*s3z - s1x*s2y*s3z)
		self.d = -(s2z*s3y*t1y - s2y*s3z*t1y - s1z*s3y*t2y + s1y*s3z*t2y + s1z*s2y*t3y - s1y*s2z*t3y)/(-s1z*s2y*s3x + s1y*s2z*s3x + s1z*s2x*s3y - s1x*s2z*s3y - s1y*s2x*s3z + s1x*s2y*s3z)
		self.e = -(s2z*s3x*t1y - s2x*s3z*t1y - s1z*s3x*t2y + s1x*s3z*t2y + s1z*s2x*t3y - s1x*s2z*t3y)/(s1z*s2y*s3x - s1y*s2z*s3x - s1z*s2x*s3y + s1x*s2z*s3y + s1y*s2x*s3z - s1x*s2y*s3z)
		self.f = -(-s2y*s3x*t1y + s2x*s3y*t1y + s1y*s3x*t2y - s1x*s3y*t2y - s1y*s2x*t3y + s1x*s2y*t3y)/(s1z*s2y*s3x - s1y*s2z*s3x - s1z*s2x*s3y + s1x*s2z*s3y + s1y*s2x*s3z - s1x*s2y*s3z)
		self.g = -(s2z*s3y*t1z - s2y*s3z*t1z - s1z*s3y*t2z + s1y*s3z*t2z + s1z*s2y*t3z - s1y*s2z*t3z)/(-s1z*s2y*s3x + s1y*s2z*s3x + s1z*s2x*s3y - s1x*s2z*s3y - s1y*s2x*s3z + s1x*s2y*s3z)
		self.h = -(s2z*s3x*t1z - s2x*s3z*t1z - s1z*s3x*t2z + s1x*s3z*t2z + s1z*s2x*t3z - s1x*s2z*t3z)/(s1z*s2y*s3x - s1y*s2z*s3x - s1z*s2x*s3y + s1x*s2z*s3y + s1y*s2x*s3z - s1x*s2y*s3z) 
		self.i = -(-s2y*s3x*t1z + s2x*s3y*t1z + s1y*s3x*t2z - s1x*s3y*t2z - s1y*s2x*t3z + s1x*s2y*t3z)/(s1z*s2y*s3x - s1y*s2z*s3x - s1z*s2x*s3y + s1x*s2z*s3y + s1y*s2x*s3z - s1x*s2y*s3z)
		return
	
	def __mul__(self,other):
		'''
		Returns new Vector with old properties, but transformed by transformation matrix
		'''
		protein = other.protein
		other.protein = None # must add to avoid deepcopy recursive ultra madness		
		transformed_Vector = dcp(other)
		transformed_Vector.protein = protein
		transformed_Vector = transformed_Vector.translateBy(-self.source)
		new_x = (self.a*transformed_Vector.x + self.b*transformed_Vector.y + self.c*transformed_Vector.z) 
		new_y = (self.d*transformed_Vector.x + self.e*transformed_Vector.y + self.f*transformed_Vector.z) 
		new_z = (self.g*transformed_Vector.x + self.h*transformed_Vector.y + self.i*transformed_Vector.z)
		transformed_Vector.x = round(new_x,7)
		transformed_Vector.y = round(new_y,7)
		transformed_Vector.z = round(new_z,7)
		transformed_Vector = transformed_Vector.translateBy(self.target)
		transformed_Vector.desc += 'Transformed'
		return transformed_Vector

class Plane:
	def __init__(self,a,b,c):
		self.a = a
		self.b = b
		self.c = c
		tmp = (self.a-self.b).cross(self.a-self.c)
		#tmp = self.a.sub(self.b).cross(self.a.sub(self.c))
		self.n = tmp/(tmp.length())
		#self.n.uberID = ",".join([self.a.uberID,self.b.uberID,self.c.uberID])
	# Distance to origine
		self.d = self.a.dot(self.n)
	# Plane Basis Vector = Plane Normal times distance
		self.D = self.n*self.d
		return
	
	def projectionOfVector(self,vector):
		'''self is plane and vector vector tp project onto membrane '''
		projected = p3d.vector.Vector((self.n.y**2+self.n.z**2)*vector.x + (-self.n.x*self.n.y)*vector.y + (-self.n.x*self.n.z) * vector.z,
		                                        -self.n.y*self.n.x*vector.x + (self.n.x**2+self.n.z**2)*vector.y +(-self.n.y*self.n.z)*vector.z,
		                                        -self.n.z*self.n.x*vector.x+ (-self.n.z*self.n.y)*vector.y + (self.n.x**2+self.n.y**2) * vector.z)
		return projected+self.D
		'''
		Old style here self vector and zeother plane.n
		return Vector((zeother.y**2+zeother.z**2)*self.x + (-zeother.x*zeother.y)*self.y + (-zeother.x*zeother.z) * self.z,
		                                        -zeother.y*zeother.x*self.x + (zeother.x**2+zeother.z**2)*self.y +(-zeother.y*zeother.z)*self.z,
		                                        -zeother.z*zeother.x*self.x+ (-zeother.z*zeother.y)*self.y + (zeother.x**2+zeother.y**2) * self.z)
		'''
	
	def minDistance2Vector(self,vector):
		projected = self.projectionOfVector(vector)
		toZeroPlane = vector-projected
		return abs(toZeroPlane-self.D)
	
	def info(self):
		print('Vectors used to define plane:')
		for v in [self.a,self.b,self.c]:
			print(v.info(lvl='coordinates'))
		print('Plane normal vector:',self.n.info(lvl='coordinates'))
		print('Plane basis vector:',self.D.info(lvl='coordinates'))
		print('Distance to origine:',self.d)
		return
	

def bp(message=''):
	if message != '':
		print('<ERROR>',message)
		print('exiting ...')
		sys.exit(1)
	sys.exit(0)

def dihedral(a,b,c,d):
	#startT = time.time()
	p1n = Plane(a,b,c).n
	p2n = Plane(b,c,d).n
	angle = p1n.angleBetween(p2n)
	angle_sign = (c-b).dot(p1n.cross(p2n))
	#print('<--*-->',b.info(),time.time()-startT,'s')
	if angle_sign < 0:
		return -1*angle
	else:
		return angle

def strTr( text, dic ):
	pat = "(%s)" % "|".join( map(re.escape, dic.keys()) )
	return re.sub( pat, lambda m:dic[m.group()], text)

def test():
	a = p3d.vector.Vector(1,0,0)
	b = p3d.vector.Vector(0,1,0)
	c = p3d.vector.Vector(1,1,0)
	plane = p3d.geo.Plane(a,b,c)
	plane.info() # <<< should be like .info() in the other modules ....
	print(type(plane))
	k = p3d.vector.Vector(2,2,2)
	print('geo.plane: Projecting',k.info(lvl='coordinates'),'onto plane results into',plane.projectionOfVector(k).info(lvl='coordinates'))
	exit(0)
	
if __name__ == '__main__':
	print('yes!')