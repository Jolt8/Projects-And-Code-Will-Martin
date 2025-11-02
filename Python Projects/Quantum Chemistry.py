from CoolProp.CoolProp import PropsSI
#from CoolProp.CoolProp import IProps, get_Fluid_index
#from CoolProp.CoolProp import PhaseSI
import CoolProp.CoolProp as cp
from CoolProp.HumidAirProp import HAPropsSI # type: ignore

import numpy as np

import math 

import scipy
from scipy import optimize

from pyfluids import Fluid, FluidsList, Input 

from pyscf import gto, scf, cc

dm = gto.M(atom="H 0 0 0; H 0 0 1")