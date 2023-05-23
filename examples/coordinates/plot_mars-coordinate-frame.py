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
from astropy.coordinates.builtin_frames.bodyframes import BodyBaseCoordinateFrame
from astropy.coordinates.representation import CartesianRepresentation
from astropy.coordinates.representation.geodetic import BaseGeodeticRepresentation

##############################################################################
# The first step is to create a new class, which we'll call
# ``MarsBestFit`` and make it a subclass of
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
    _wrap_angle = 360


#############################################################################
# The new planetary body-fixed reference system will use the previous defined
# representation.

class MarsCoordinateFrame(BodyBaseCoordinateFrame):
    """
    A spheroidal reference system with longitude spanning from 0 to 360 degree
    est positive.

    """

    default_representation = MarsBestFitAeroid

mars = MarsCoordinateFrame(lon=np.linspace(0, 2*np.pi, 128)*u.radian,
                           lat=np.linspace(-0.5 * np.pi, 0.5 * np.pi, 128)*u.radian)

#####################################################################################
# Now let's define a new geodetic representation obtained from MarsBestFitAeroid but
# described by planetocentric latitudes.

class MarsBestFitOcentricAeroid(BaseGeodeticRepresentation):
    """
    A Spheroidal planetocentric representation of Mars that minimized deviations with
    respect to the areoid following
        Ardalan A. A, R. Karimi, and E. W. Grafarend (2010)
        https://doi.org/10.1007/s11038-009-9342-7
    """

    _equatorial_radius = 3395.4280 * u.km
    _flattening = 0.5227617843759314 * u.percent
    _wrap_angle = 360
    _ographic = False

#############################################################################
# and a related planetary body-fixed reference system.

class MarsCoordinateOcentricFrame(BodyBaseCoordinateFrame):
    """
    A spheroidal reference system with longitude spanning from 0 to 360 degree
    est positive.

    """

    default_representation = MarsBestFitOcentricAeroid

#############################################################################
# As a comparison we define a new spherical frame representation,

class MarsSphere(BaseGeodeticRepresentation):
    """
    A Spherical representation of Mars
    """

    _equatorial_radius = 3395.4280 * u.km
    _flattening = 0.0 * u.percent
    _wrap_angle = 360

#############################################################################
# and the corresponding reference frame.

class MarsSphereCoordinateFrame(BodyBaseCoordinateFrame):
    """
    A spherical reference system with longitude spanning from 0 to 360 degree
    east positive.

    """

    default_representation = MarsSphere

#############################################################################
# Now we plot the differences between each component of the cartesian
# representation with respect to the spherical model, assuming the point on the
# surface of the body (``height = 0``)

mars_sphere = MarsSphereCoordinateFrame(lon=np.linspace(0, 2*np.pi, 128)*u.radian,
                           lat=np.linspace(-0.5 * np.pi, 0.5 * np.pi, 128)*u.radian)
mars = MarsCoordinateFrame(lon=np.linspace(0, 2*np.pi, 128)*u.radian,
                           lat=np.linspace(-0.5 * np.pi, 0.5 * np.pi, 128)*u.radian)
mars_ocentric = MarsCoordinateOcentricFrame(lon=np.linspace(0, 2*np.pi, 128)*u.radian,
                           lat=np.linspace(-0.5 * np.pi, 0.5 * np.pi, 128)*u.radian)

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
