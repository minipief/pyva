# -*- coding: utf-8 -*-
"""
Created on Sun Mar  6 12:01:53 2022

"""

from setuptools import setup,find_packages

setup(
    name = "pyva",
    version = "1.0.0",
    author = 'Dr. Alexander Peiffer',
    author_email = "author@alexanderpeiffer.de",
    py_modules = ["models","useful"],
    packages = find_packages(),
    install_requires = ['pint>=0.17','numpy>=1.20.2','matplotlib>=3.3.4','scipy>=1.6.2','pandas>=1.2.5']
    )

