# -*- coding: utf-8 -*-
"""
Created on Sun Mar  6 12:01:53 2022

"""

from setuptools import setup,find_packages

setup(
    name = "pyva",
    version = "1.0",
    author = 'Dr. Alexander Peiffer',
    author_email = "author@alexanderpeiffer.de",
    py_modules = ["models","useful"],
    packages = find_packages(),
    install_requires = ['pint>=0.17','numpy>=1.202','matplotlib>=3.3.4','scipy','pandas']
    )

