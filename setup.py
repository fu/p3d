#!/usr/bin/env python

from distutils.core import setup
setup(name='p3d',
      version='0.4.3',
      packages = ['p3d'],
      package_dir = {'p3d': 'src'},
      description='protein structure module',
      long_description='p3d - python module for structural bioinformatics',
      author='Christian Fufezan & Michael Specht',
      author_email='p3d@fufezan.net',
      url='http://p3d.fufezan.net',
      license='GNU General Public License (GPL)',
      platforms='any that supports python 3+',
      classifiers=[
	       'Development Status :: 4 - Beta',
	       'Environment :: Console',
	       'Intended Audience :: Education',
	       'Intended Audience :: Science/Research',
	       'Intended Audience :: Developers',
	       'License :: OSI Approved :: GNU General Public License (GPL)',
	       'Operating System :: MacOS :: MacOS X',
	       'Operating System :: Microsoft :: Windows',
	       'Operating System :: POSIX',
           'Operating System :: POSIX :: SunOS/Solaris',
           'Operating System :: Unix'
	       'Programming Language :: Python :: 3',
	       'Topic :: Scientific/Engineering :: Bio-Informatics',
	       'Topic :: Scientific/Engineering :: Chemistry',
	       'Topic :: Scientific/Engineering :: Medical Science Apps.',
	       'Topic :: Software Development :: Libraries :: Python Modules'
	          ],
      )
