# -*- coding: utf-8 -*-
#py.test tests/AnalyticCalc1_test.py
#C:\d from old\RadiaBeam\RadSoft\radtrack\tests>py.test AnalyticCalc1_test.py
u"""PyTest for :mod:`radiasoft.AnalyticCalc`

:copyright: Copyright (c) 2015 Bivio Software, Inc.  All Rights Reserved.
:license: http://www.apache.org/licenses/LICENSE-2.0.html
"""
from __future__ import absolute_import, division, print_function, unicode_literals
from io import open

import os
import numpy as np
import pytest

from pykern.pkdebug import pkdc, pkdp
from pykern import pkunit
from pykern import pkyaml

from radtrack.srw import AnalyticCalc

_EPSILON = 1e-300

def _run(fname):
    y=ReadYaml(fname)
    for fn in y:
        f = getattr(AnalyticCalc, fn)
        for case in y[fn]:
            kwargs = case['Input']
            actual = f(**kwargs)
            print (['***', fname])
            print (actual)
            print (case['Expect'])
            _assert_array(case['Expect'], actual)

def _assert(expect, actual, expected_error=0.01):
    if _EPSILON > abs(expect):
        assert _EPSILON > abs(actual)
        return
    elif _EPSILON > abs(actual):
        raise AssertionError(
            'expect {} != {} actual'.format(expect, actual))
    assert expected_error > abs(expect/actual - 1)

def _assert_array(expect, actual):
    if np.shape(expect):
        for e, a in zip(expect, actual):
            _assert(e, a)
    else:
        _assert(expect, actual)

def ReadYaml(filename):
    d = pkunit.data_dir()
    y = pkyaml.load_file(d.join(filename))
    return y

#ReadYaml('nsls2.yml')
def test_1():
    _run('IDWaveLengthPhotonEnergy.yml')
    _run('CriticalEnergyWiggler.yml')
    _run('RadiatedPowerPlanarWiggler.yml')
    _run('CentralPowerDensityPlanarWiggler.yml')
    _run('UndulatorSourceSizeDivergence.yml')
    _run('SpectralFlux.yml')
    _run('SpectralCenBrightness.yml')
    _run('UndulatorAngleCoordinateOscillation.yml')
#    assert 0

