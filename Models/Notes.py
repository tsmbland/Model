"""
Model structure:

- class Params

- class Model
    - __init__(self, params, starts)
    - get_all(self)
    - run(self, res)

- class Res
    - __init__(self, params)
    - update(self, t, res)
    - save(self, direc, compression)


To run:
    p = x.Params(...)
    m = x.Model(params)
    r = x.Res(params, starts)
    res = x.Model.run(r)
    res.save(direc, compression)

"""
