"""
180725
Testing m0000

"""

import M as x
import Models.m0010 as m0010

x.datadirec = '../../ModelData'


########## Single simulaton <- good

p = m0010.p0
model = m0010.Model(p)
x.alg_singlesim(model, 2, 0, 0, compression=0, funcs=x.all_analysis)
x.sliderplot(2, 0, 0)

res = x.loaddata(2, 0, 0)
x.plt.plot(res.atot)
x.plt.plot(res.ptot)
x.plt.show()


# model = m.Model(m.p0)
# model.params.pA = 0
# alg_singlesim(model, 2, 0, 1, compression=0)
# sliderplot(2, 0, 1)





