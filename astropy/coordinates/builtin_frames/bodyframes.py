# Licensed under a 3-clause BSD style license - see LICENSE.rst


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
    default_representation : `~astropy.coordinates.BaseGeodeticRepresentation`, keyword-only
        The geodetic representation used to describe the body.
"""


doc_footer_geo = """
    Other parameters
    ----------------
    object_name : str, optional
        The name of the body the frame is attached to.
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
    """

    default_representation = BaseGeodeticRepresentation

    obstime = TimeAttribute(default=None)
    object_name = None
