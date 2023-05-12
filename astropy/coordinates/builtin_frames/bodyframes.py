# Licensed under a 3-clause BSD style license - see LICENSE.rst


import astropy.units as u
from astropy.coordinates.attributes import TimeAttribute
from astropy.coordinates.baseframe import BaseCoordinateFrame, base_doc
from astropy.coordinates.representation.geodetic import (
    BaseGeodeticRepresentation,
)

from astropy.utils.decorators import format_doc

__all__ = ["BodyBaseCoordinateFrame"]

doc_components_body = """
    lon : `~astropy.coordinates.Angle`, optional, keyword-only
        The geodetic longitude for this object (``lat`` must also be given).
    lat : `~astropy.coordinates.Angle`, optional, keyword-only
        The geodetic latitude for this object (``lon`` must also be given).
    height : `~astropy.units.Quantity` ['length'], optional, keyword-only
        The distance for this object from the surface.
    representation : `~astropy.coordinates.BaseGeodeticRepresentation`, keyword-only
        The geodetic representation used to describe the body.
"""


doc_footer_geo = """
    Other parameters
    ----------------
    obstime : `~astropy.time.Time`, optional
        The time at which the observation is taken.  Used for determining the
        position of the Earth. Defaults to J2000.
"""


@format_doc(
    base_doc, components=doc_components_body.format("specified location"), footer=""
)
class BodyBaseCoordinateFrame(BaseCoordinateFrame):
    """
    Base class for body coordinate frames.

    This class is not intended to be used directly and has no transformations defined.

    * Defines the frame attribute ``obstime`` for observation time.
    * Defines a default wrap angle of 180 degrees for longitude in geodetic representation,
      which can be overridden via the class variable ``_wrap_angle``.
    * Defines a default positive direction East for longitude in geodetic representation,
      which can be overridden via the class variable ``_positive_longitude``.
    """

    default_representation = BaseGeodeticRepresentation

    obstime = TimeAttribute(default=None)

    _wrap_angle = 180 * u.deg  # for longitude
    _positive_longitude = "east"  # for longitude

    def __init__(self, *args, **kwargs):
        self.object_name = None

        # If wrap_longitude=False is passed in, do not impose a specific wrap angle for the frame
        if not kwargs.pop("wrap_longitude", True):
            self._wrap_angle = None

        super().__init__(*args, **kwargs)
        if self.has_data:
            self._set_data_lon_wrap_angle(self.data)

    @staticmethod
    def _set_data_lon_wrap_angle(data):
        if hasattr(data, "lon"):
            data.lon.wrap_angle = 180.0 * u.deg
        return data

    def represent_as(self, base, s="base", in_frame_units=False):
        """
        Ensure the wrap angle for any geodetic
        representations.
        """
        data = super().represent_as(base, s, in_frame_units=in_frame_units)
        self._set_data_lon_wrap_angle(data)
        return data
