# Licensed under a 3-clause BSD style license - see LICENSE.rst


from astropy.coordinates.attributes import CoordinateAttribute, TimeAttribute
from astropy.coordinates.baseframe import BaseCoordinateFrame, base_doc
from astropy.coordinates.representation import CartesianRepresentation
from astropy.utils.decorators import format_doc

from .icrs import ICRS

__all__ = ["BaseBodyCoordinateFrame"]


doc_footer = """
    Other parameters
    ----------------
    naif_code : int, optional
        The body NAIF code as defined in
        <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/naif_ids.html#NAIF%20Object%20ID%20numbers>_
    object_name : str, optional
        The name of the body the frame is attached to.
    observer : `~astropy.coordinates.SkyCoord`, optional
        The coordinates of the observer.
    obstime : `~astropy.time.Time`, optional
        The time at which the observation is taken.  Used for determining the
        position of the planetary body. Defaults to J2000.
"""


@format_doc(base_doc, components="", footer=doc_footer)
class BaseBodyCoordinateFrame(BaseCoordinateFrame):
    """
    Base class for body coordinate frames.

    This class is not intended to be used directly and has no transformations defined.

    * Defines the frame attribute ``obstime`` for observation time.
    """

    default_representation = CartesianRepresentation

    naif_code = None
    object_name = None
    observer = CoordinateAttribute(ICRS, default=None)
    obstime = TimeAttribute(default=None)
