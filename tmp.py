from setuptools import setup, find_packages
import setuptools
import setuptools.command.develop
import setuptools.command.install
import subprocess
import pkg_resources
import sys, os
from setuptools import Extension

def compileKepcart():
    loc = pkg_resources.Environment()['mpcutilities'][0].location
    cwd = os.path.join(loc, 'mpcutilities')
    print("loc ", loc)
    print("cwd ", cwd)
    _ = subprocess.check_call(
                          'bash README_COMPILE.txt',
                          cwd=os.path.join(loc, 'mpcutilities'),
                          shell=True)
    print("_",_)

compileKepcart()
