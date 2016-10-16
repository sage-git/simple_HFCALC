#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
env = Environment(ENV=os.environ)
IsDebug = ARGUMENTS.get('debug', 0)
Use_my_ERF = ARGUMENTS.get('my_erf', 1)
env.Append(F90FLAGS=["-cpp"])
if int(IsDebug) == 1:
    env.Append(F90FLAGS=["-g", "-fbacktrace"])
else:
    env.Append(F90FLAGS=["-O2"])

if int(Use_my_ERF) == 1:
    env.Append(F90FLAGS=["-DUSE_MY_ERF"])

env.Program(target='szabo_test', source=Glob('*.f90'))
