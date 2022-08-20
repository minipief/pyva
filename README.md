# pyva
Python toolbox for vibroacoustics

This toolbox is the code extension for the book 
A. Peiffer: Vibroacoustic Simulation: In Introduction to Statistical Energy Analysis and Hybrid Methods, John Wiley, 2022

[Wiley book details](https://www.wiley.com/en-us/Vibroacoustic+Simulation%3A+An+Introduction+to+Statistical+Energy+Analysis+and+Hybrid+Methods-p-9781119849841)

[Order book](http://www.wiley.com/buy/9781119849841)

This toolbox allows the reader to follow the examples in the book and people that are interested in vibroacoustic simulation with a focus on 
statistical energy analysis (SEA) and hybrid methods (hybrid FEM/SE), to perform SEA simulation without heavy licence costs.

# Documentation

Direct link to the documentation can be found on [pyva.eu](https://pyva.eu)

# Contributing to pyva

The current version is the baseline for further extensions. It was mainly driven by providing examples and test cases for the bool
on [Vibroacoustic Simulation](http://www.wiley.com/buy/9781119849841). Contributions are highly welcome.
The following extensions shall be included mid-term

A major step is a GUI development! Here, suggestions would be very helpful. I am currently thinking about a combination of pyqt 
and vtk, eventually integrating cadquery. However, I am still scanning the options. 

## Main Module

- Engineering units contraints, e.g. setting pressure or velocity level of one SEA subsystem

## Coupling Package

- point junctions
- inclusion of beams in line junctions

## Properties Package

### Geometrical Property Module

- Further Beam Cross Sections

### Material Module

- Anisotropic materials
- Porous Absorber Models (Biot - JCA)

### Structural Property Module

- Linear Laminates
- Sandwich

## Systems Package

- Singly curved shells
- Doubly curved shells
- Infinite layers with porous Biot material
- Elastic Solid as infinite layer

## Load package

- Diffuse wave field as load
- Cross spectral density excitation
 
# Author page and contact information

If you are not sure if pyva suits your needs in the current status feel free to contact me via.
[Author Page - Alexander Peiffer](https://docpeiffer.com)


