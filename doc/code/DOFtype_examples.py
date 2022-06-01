import pyva.data.dof as dof
import copy

force = dof.DOFtype(typestr='force')

area = dof.DOFtype(typestr='area')

pressure = force/area
p2 = copy.deepcopy(pressure)

p2._exp = 2

print(p2)

p2_per_frequency = dof.DOFtype(typestr = 'pressure',exponent = 2, xdata_exponent = 1, xtypestr = 'frequency')

freq      = dof.DOFtype(typestr = 'frequency')
pressure_ = dof.DOFtype(typestr = 'pressure')

p2_ = p2_per_frequency*freq

no_unit = pressure_/pressure

print(no_unit)

# switch to dynamic matrix






