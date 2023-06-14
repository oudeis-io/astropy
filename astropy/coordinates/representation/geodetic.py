# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import erfa

from astropy import units as u
from astropy.coordinates.angles import Latitude, Longitude
from astropy.utils.decorators import format_doc

from .base import BaseRepresentation
from .cartesian import CartesianRepresentation


class WestLongitudeMixin:
    def to_cartesian(self):
        self.lon = self.lon.value
        return super().to_cartesian()

    @classmethod
    def from_cartesian(cls, cart):
        angular = super().from_cartesian(cart)
        return cls(
            Longitude(
                -1 * angular.lon.value,
                u.deg,
                wrap_angle=angular.lon.wrap_angle,
                copy=False,
            ),
            angular.lat.to(u.deg),
            angular.height,
            copy=False,
        )

    @property
    def lon(self):
        return self._lon

    @lon.setter
    def lon(self, value):
        self._lon = Longitude(
            -1 * value, u.deg, wrap_angle=self._lon.wrap_angle, copy=False
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
    They can be made west positive setting ``_positive_longitude='west'``, or spanning
    from 0 to 360 degrees setting ``_wrap_angle=360 * u.deg``.
    """

    attr_classes = {"lon": Longitude, "lat": Latitude, "height": u.Quantity}
    _wrap_angle = 180 * u.deg

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
        if not hasattr(
            cls._wrap_angle, "unit"
        ) or not cls._wrap_angle.unit.is_equivalent(u.deg):
            raise u.UnitTypeError("Attribute _wrap_angle requires deg units.")
        super().__init_subclass__(**kwargs)

    def __init__(self, lon, lat=None, height=None, copy=True):
        if height is None and not isinstance(lon, self.__class__):
            height = 0 << u.m

        if not lon.unit.is_equivalent(u.deg) or not lat.unit.is_equivalent(u.deg):
            raise u.UnitTypeError(
                f"{self.__class__.__name__} requires lat and lon with units of angle."
            )
        lon = Longitude(lon << u.deg, wrap_angle=self._wrap_angle)
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
        xyz = erfa.gd2gce(
            self._equatorial_radius,
            self._flattening,
            self.lon,
            self.lat,
            self.height,
        )
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
        return cls(
            Longitude(lon << u.deg, wrap_angle=cls._wrap_angle),
            Latitude(lat << u.deg),
            height,
            copy=False,
        )


@format_doc(geodetic_base_doc)
class BaseBodycentricRepresentation(BaseRepresentation):
    """Representation of points in bodycentric 3D coordinates."""

    attr_classes = {"lon": Longitude, "lat": Latitude, "height": u.Quantity}
    _wrap_angle = 180 * u.deg

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
        super().__init_subclass__(**kwargs)

    def __init__(self, lon, lat=None, height=None, copy=True):
        if height is None and not isinstance(lon, self.__class__):
            height = 0 << u.m

        super().__init__(lon, lat, height, copy=copy)
        if not self.height.unit.is_equivalent(u.m):
            raise u.UnitTypeError(
                f"{self.__class__.__name__} requires height with units of length."
            )
        if not self.lon.unit.is_equivalent(u.deg) or not self.lat.unit.is_equivalent(
            u.deg
        ):
            raise u.UnitTypeError(
                f"{self.__class__.__name__} requires lat and lon with units of angle."
            )

    def to_cartesian(self):
        """
        Converts bodycentric coordinates to 3D rectangular (geocentric)
        cartesian coordinates.
        """
        x_spheroid = self._equatorial_radius * np.cos(self.lat) * np.cos(self.lon)
        y_spheroid = self._equatorial_radius * np.cos(self.lat) * np.sin(self.lon)
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
        bodycentric coordinates.
        """
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
        lon = np.where(
            xyz[0] != 0.0,
            np.arctan(xyz[1] / xyz[0]),
            np.sign(xyz[1]) * 0.5 * np.pi * u.radian,
        )
        return cls(
            Longitude(lon << u.deg, wrap_angle=cls._wrap_angle),
            Latitude(lat << u.deg),
            height,
            copy=False,
        )


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
