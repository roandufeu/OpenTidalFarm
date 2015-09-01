"""
Initialises the various functionals. All functionals should overload the
Prototype, as this allows them to be combined (added, subtracted and scaled).
Functionals are grouped by type; Power / Cost / Environment etc.
"""

from sediment_functionals import SedimentFunctional
from ecology_functionals import EcologyFunctional
from fishing_functionals import FishingFunctional
from power_functionals import PowerFunctional
from cost_functionals import CostFunctional
from prototype_functional import PrototypeFunctional
from time_integrator import TimeIntegrator
