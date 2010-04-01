#!/bin/bash

## create all packages

python setup.py sdist --formats=bztar,gztar,zip
#python setup.py bdist --formats=bztar,gztar,bztar,zip
