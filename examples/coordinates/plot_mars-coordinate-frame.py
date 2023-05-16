r"""
============================================
Create a new coordinate frame class for Mars
============================================

This example describes how to subclass and define a custom geodetic
coordinate frame, as discussed in :ref:`astropy:astropy-coordinates-design` and
:ref:`astropy-coordinates-create-geodetic`.

To do this, first we need to define a subclass of a
`~astropy.coordinates.BaseGeodeticRepresentation`, then a subclass of
`~astropy.coordinates.BodyBaseCoordinateFrame` using the previous defined
representation.

*By: Chiara Marmo*

*License: BSD*


"""

##############################################################################
# Make `print` work the same in all versions of Python, set up numpy,
# matplotlib, and use a nicer set of plot parameters:

import matplotlib.pyplot as plt
import numpy as np

from astropy.visualization import astropy_mpl_style

plt.style.use(astropy_mpl_style)


##############################################################################
# Import the packages necessary for coordinates

import astropy.units as u
from astropy.coordinates.representation import CartesianRepresentation
from astropy.coordinates.representation.geodetic import BaseGeodeticRepresentation
from astropy.coordinates.builtin_frames.bodyframes import BodyBaseCoordinateFrame

##############################################################################
# The first step is to create a new class, which we'll call
# ``MarsBestFit`` and make it a subclass of
# `~astropy.coordinates.BaseGeodeticRepresentation`.
# It represent a best fit of the Mars spheroid to the martian geoid (areoid):

class MarsBestFitAeroid(BaseGeodeticRepresentation):
    """
    A Spheroidal representation of Mars that minimized deviations with respect to the
    areoid following
        Ardalan A. A, R. Karimi, and E. W. Grafarend (2010)
        https://doi.org/10.1007/s11038-009-9342-7
    """

    _equatorial_radius = 3395.4280 * u.km
    _flattening = 0.5227617843759314 * u.percent


#############################################################################
# The new planetary body-fixed reference system will use the previous defined
# representation with a longitude spanning from 0 to 360 degrees.

class MarsCoordinateFrame(BodyBaseCoordinateFrame):
    """
    A spheroidal reference system with longitude spanning from 0 to 360 degree
    est positive.

    """

    default_representation = MarsBestFitAeroid
    wrap_angle = 360.0 * u.deg

mars = MarsCoordinateFrame(lon=np.linspace(0, 2*np.pi, 128)*u.radian,
                           lat=np.linspace(-0.5 * np.pi, 0.5 * np.pi, 128)*u.radian)

#############################################################################
# As a comparison we define a new spherical frame representation,

class MarsSphere(BaseGeodeticRepresentation):
    """
    A Spherical representation of Mars
    """

    _equatorial_radius = 3395.4280 * u.km
    _flattening = 0.0 * u.percent

#############################################################################
# and the corresponding reference frame.

class MarsSphereCoordinateFrame(BodyBaseCoordinateFrame):
    """
    A spheric reference system with longitude spanning from 0 to 360 degree
    est positive.

    """

    default_representation = MarsSphere
    wrap_angle = 360.0 * u.deg

#############################################################################
# Now we plot the differences between each component of the cartesian
# representation

mars = MarsCoordinateFrame(lon=np.linspace(0, 2*np.pi, 128)*u.radian,
                           lat=np.linspace(-0.5 * np.pi, 0.5 * np.pi, 128)*u.radian)

mars_sphere = MarsSphereCoordinateFrame(lon=np.linspace(0, 2*np.pi, 128)*u.radian,
                           lat=np.linspace(-0.5 * np.pi, 0.5 * np.pi, 128)*u.radian)

xyz = mars.represent_as(CartesianRepresentation)
xyz_sphere = mars_sphere.represent_as(CartesianRepresentation)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

ax.scatter((xyz._x - xyz_sphere._x) << u.km,
           (xyz._y - xyz_sphere._y) << u.km,
           (xyz._z - xyz_sphere._z) << u.km)
ax.tick_params(labelsize=8)
ax.set(xlabel='x Difference [Km]',
       ylabel='y Difference [Km]',
       zlabel='z difference [Km]')

ax.set_title("Mars spheroid comparison")

plt.show()
