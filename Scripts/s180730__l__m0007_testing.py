"""
180725
Testing m0007

"""

import M as x
import Models.m0007 as m0007

x.datadirec = '../../ModelData'

######### Single simulation <- good

model = m0007.Model(m0007.p0)
x.alg_singlesim(model, 6, 0, 0, compression=0, funcs=[x.mse, x.asi_a, x.asi_p])
x.sliderplot(6, 0, 0)