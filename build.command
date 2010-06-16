#!/bin/bash

## create all packages
rm MANIFEST
python setup.py sdist --formats=bztar,gztar,zip
#python setup.py bdist --formats=bztar,gztar,bztar,zip
cd dist
tar xvfj *.bz2
