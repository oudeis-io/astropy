r"""
============================================
Create a new coordinate frame class for Mars
============================================

This example describes how to subclass and define a custom coordinate frame for a
planetary body which can be described by a geodetic or bodycentric representation,
as discussed in :ref:`astropy:astropy-coordinates-design` and
:ref:`astropy-coordinates-create-geodetic`.

To do this, first we need to define a subclass of a
`~astropy.coordinates.BaseGeodeticRepresentation` and
`~astropy.coordinates.BaseBodycentricRepresentation`, then a subclass of
`~astropy.coordinates.BaseBodyCoordinateFrame` using the previous defined
representations.

*By: Chiara Marmo*

*License: BSD*


"""

##############################################################################
# Set up numpy, matplotlib, and use a nicer set of plot parameters:

import matplotlib.pyplot as plt
import numpy as np

from astropy.visualization import astropy_mpl_style

plt.style.use(astropy_mpl_style)


##############################################################################
# Import the packages necessary for coordinates

import astropy.units as u
from astropy.coordinates.builtin_frames.bodyframes import BaseBodyCoordinateFrame
from astropy.coordinates.representation import CartesianRepresentation
from astropy.coordinates.representation.geodetic import (
    BaseBodycentricRepresentation,
    BaseGeodeticRepresentation,
)

##############################################################################
# The first step is to create a new class, and make it a subclass of
# `~astropy.coordinates.BaseGeodeticRepresentation`.
# Geodetic latitudes are used and longitudes span from 0 to 360 degrees east positive
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

#####################################################################################
# Now let's define a new geodetic representation obtained from MarsBestFitAeroid but
# described by planetocentric latitudes.

class MarsBestFitOcentricAeroid(BaseBodycentricRepresentation):
    """
    A Spheroidal planetocentric representation of Mars that minimized deviations with
    respect to the areoid following
        Ardalan A. A, R. Karimi, and E. W. Grafarend (2010)
        https://doi.org/10.1007/s11038-009-9342-7
    """

    _equatorial_radius = 3395.4280 * u.km
    _flattening = 0.5227617843759314 * u.percent

#############################################################################
# As a comparison we define a new spherical frame representation, we could
# have based it on `~astropy.coordinates.BaseBodycentricRepresentation` too.

class MarsSphere(BaseGeodeticRepresentation):
    """
    A Spherical representation of Mars
    """

    _equatorial_radius = 3395.4280 * u.km
    _flattening = 0.0 * u.percent

#############################################################################
# The new planetary body-fixed reference system will be described using the
# previous defined representations.

class MarsCoordinateFrame(BaseBodyCoordinateFrame):
    """
    A reference system for Mars.
    """

    object_name = "Mars"


#############################################################################
# Now we plot the differences between each component of the cartesian
# representation with respect to the spherical model, assuming the point on the
# surface of the body (``height = 0``)

mars_sphere = MarsCoordinateFrame(lon=np.linspace(0, 2*np.pi, 128)*u.radian,
                           lat=np.linspace(-0.5 * np.pi, 0.5 * np.pi, 128)*u.radian,
                           representation_type=MarsSphere)
mars = MarsCoordinateFrame(lon=np.linspace(0, 2*np.pi, 128)*u.radian,
                           lat=np.linspace(-0.5 * np.pi, 0.5 * np.pi, 128)*u.radian,
                           representation_type=MarsBestFitAeroid)
mars_ocentric = MarsCoordinateFrame(lon=np.linspace(0, 2*np.pi, 128)*u.radian,
                           lat=np.linspace(-0.5 * np.pi, 0.5 * np.pi, 128)*u.radian,
                           representation_type=MarsBestFitOcentricAeroid)

xyz_sphere = mars_sphere.represent_as(CartesianRepresentation)
xyz = mars.represent_as(CartesianRepresentation)
xyz_ocentric = mars_ocentric.represent_as(CartesianRepresentation)

fig, ax = plt.subplots(2, subplot_kw={"projection": "3d"})

ax[0].scatter((xyz_sphere._x - xyz._x) << u.km,
              (xyz_sphere._y - xyz._y) << u.km,
              (xyz_sphere._z - xyz._z) << u.km)
ax[0].tick_params(labelsize=8)
ax[0].set(xlabel='x [Km]',
          ylabel='y [Km]',
          zlabel='z [Km]')

ax[0].set_title("Mars sphere - odetic spheroid difference")

ax[1].scatter((xyz_sphere._x - xyz_ocentric._x) << u.km,
              (xyz_sphere._y - xyz_ocentric._y) << u.km,
              (xyz_sphere._z - xyz_ocentric._z) << u.km)
ax[1].tick_params(labelsize=8)
ax[1].set(xlabel='x [Km]',
          ylabel='y [Km]',
          zlabel='z [Km]')

ax[1].set_title("Mars sphere - ocentric spheroid difference")

plt.show()
