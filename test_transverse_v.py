#file located in Gaia/
#filename = test_transverse_v.py

from complete_code_velocity_distribution_2D.py import transverse_velocity
import pytest

#assert checks true/false
def test_transverse_v_1():
    assert transverse_velocity(1) == 'Transverse velocity bigger than 0'