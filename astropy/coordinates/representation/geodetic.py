# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import erfa

from astropy import units as u
from astropy.coordinates.angles import Latitude, Longitude
from astropy.utils.decorators import format_doc

from .base import BaseRepresentation
from .cartesian import CartesianRepresentation


def _compute_longitude(longitude, positive_longitude, wrap_angle):
    if positive_longitude == "east":
        orient = 1
    elif positive_longitude == "west":
        orient = -1
    return Longitude(
        orient * longitude, u.deg, wrap_angle=wrap_angle * u.deg, copy=False
    )


ELLIPSOIDS = {}
"""Available ellipsoids (defined in erfam.h, with numbers exposed in erfa)."""
# Note: they get filled by the creation of the geodetic classes.


geodetic_base_doc = """{__doc__}

    Parameters
    ----------
    lon, lat : angle-like
        The longitude and latitude of the point(s), in angular units. The
        latitude should be between -90 and 90 degrees, and the longitude will
        be wrapped to an angle between 0 and 360 degrees. These can also be
        instances of `~astropy.coordinates.Angle` and either
        `~astropy.coordinates.Longitude` not `~astropy.coordinates.Latitude`,
        depending on the parameter.

    height : `~astropy.units.Quantity` ['length']
        The height to the point(s).

    copy : bool, optional
        If `True` (default), arrays will be copied. If `False`, arrays will
        be references, though possibly broadcast to ensure matching shapes.
"""


@format_doc(geodetic_base_doc)
class BaseGeodeticRepresentation(BaseRepresentation):
    """
    Base class for geodetic representations.

    Subclasses need to set attributes ``_equatorial_radius`` and ``_flattening``
    to quantities holding correct values (with units of length and dimensionless,
    respectively), or alternatively an ``_ellipsoid`` attribute to the relevant ERFA
    index (as passed in to `erfa.eform`).
    Longitudes are east positive and span from -180 to 180 degrees by default.
    They can be made west positive setting `_positive_longitude='west'`, or spanning
    from 0 to 360 degrees setting `_wrap_angle=360`.
    Planetocentric latitudes can be obtained setting `_ographic=False`.
    """

    attr_classes = {"lon": Longitude, "lat": Latitude, "height": u.Quantity}
    _positive_longitude = "east"
    _wrap_angle = 180
    _ographic = True

    def __init_subclass__(cls, **kwargs):
        if "_ellipsoid" in cls.__dict__:
            equatorial_radius, flattening = erfa.eform(getattr(erfa, cls._ellipsoid))
            cls._equatorial_radius = equatorial_radius * u.m
            cls._flattening = flattening * u.dimensionless_unscaled
            ELLIPSOIDS[cls._ellipsoid] = cls
        elif (
            "_equatorial_radius" not in cls.__dict__
            or "_flattening" not in cls.__dict__
        ):
            raise AttributeError(
                f"{cls.__name__} requires '_ellipsoid' or '_equatorial_radius' and '_flattening'."
            )
        if cls._positive_longitude not in ["east", "west"]:
            raise ValueError(
                f"Invalid argument '{cls._positive_longitude}' for '_positive_logitude' "
                "attribute: valid argument are 'east' or 'west'."
            )
        if cls._wrap_angle not in [180, 360]:
            raise ValueError(
                f"Invalid argument '{cls._wrap_angle}' for '_wrap_angle' attribute: "
                "valid arguments are 180 or 360."
            )
        super().__init_subclass__(**kwargs)

    def __init__(self, lon, lat=None, height=None, copy=True):
        if height is None and not isinstance(lon, self.__class__):
            height = 0 << u.m

        super().__init__(lon, lat, height, copy=copy)
        if not self.height.unit.is_equivalent(u.m):
            raise u.UnitTypeError(
                f"{self.__class__.__name__} requires height with units of length."
            )

    def to_cartesian(self):
        """
        Converts geodetic coordinates to 3D rectangular (geocentric)
        cartesian coordinates.
        """
        lon = _compute_longitude(self.lon, self._positive_longitude, self._wrap_angle)

        if self._ographic:
            xyz = erfa.gd2gce(
                self._equatorial_radius,
                self._flattening,
                lon,
                self.lat,
                self.height,
            )
        else:
            x_spheroid = self._equatorial_radius * np.cos(self.lat) * np.cos(lon)
            y_spheroid = self._equatorial_radius * np.cos(self.lat) * np.sin(lon)
            z_spheroid = (self._equatorial_radius * (1 - self._flattening)) * np.sin(
                self.lat
            )
            r = (
                np.sqrt(
                    x_spheroid * x_spheroid
                    + y_spheroid * y_spheroid
                    + z_spheroid * z_spheroid
                )
                + self.height
            )
            x = r * np.cos(self.lon) * np.cos(self.lat)
            y = r * np.sin(self.lon) * np.cos(self.lat)
            z = r * np.sin(self.lat)
            xyz = np.stack([x, y, z], axis=1) << u.m
        return CartesianRepresentation(xyz, xyz_axis=-1, copy=False)

    @classmethod
    def from_cartesian(cls, cart):
        """
        Converts 3D rectangular cartesian coordinates (assumed geocentric) to
        geodetic coordinates.
        """
        # Compute geodetic/planetodetic angles
        lon, lat, height = erfa.gc2gde(
            cls._equatorial_radius, cls._flattening, cart.get_xyz(xyz_axis=-1)
        )
        lon = _compute_longitude(lon, cls._positive_longitude, cls._wrap_angle)
        if not cls._ographic:
            # Compute planetocentric angles
            xyz = cart.get_xyz()
            # Compute planetocentric latitude
            p = np.sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1])
            d = np.sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2])
            lat = np.where(
                p != 0.0,
                np.arctan(xyz[2] / p),
                np.sign(xyz[2]) * 0.5 * np.pi * u.radian,
            )
            p_spheroid = cls._equatorial_radius * np.cos(lat)
            z_spheroid = (cls._equatorial_radius * (1 - cls._flattening)) * np.sin(lat)
            r_spheroid = np.sqrt(p_spheroid * p_spheroid + z_spheroid * z_spheroid)
            height = np.where(
                p_spheroid != 0.0,
                (d - r_spheroid),
                (np.abs(xyz[2]) - np.abs(z_spheroid)),
            )

        return cls(lon, lat, height, copy=False)


@format_doc(geodetic_base_doc)
class WGS84GeodeticRepresentation(BaseGeodeticRepresentation):
    """Representation of points in WGS84 3D geodetic coordinates."""

    _ellipsoid = "WGS84"


@format_doc(geodetic_base_doc)
class WGS72GeodeticRepresentation(BaseGeodeticRepresentation):
    """Representation of points in WGS72 3D geodetic coordinates."""

    _ellipsoid = "WGS72"


@format_doc(geodetic_base_doc)
class GRS80GeodeticRepresentation(BaseGeodeticRepresentation):
    """Representation of points in GRS80 3D geodetic coordinates."""

    _ellipsoid = "GRS80"
