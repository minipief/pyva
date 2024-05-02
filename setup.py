# -*- coding: utf-8 -*-
"""
Created on Sun Mar  6 12:01:53 2022

"""

from setuptools import setup,find_packages

setup(
    name = "pyva-toolbox",
    version = "1.2.3",
    author = 'Dr. Alexander Peiffer',
    author_email = "author@alexanderpeiffer.de",
    py_modules = ["models","useful"],
    packages = find_packages()
    )

