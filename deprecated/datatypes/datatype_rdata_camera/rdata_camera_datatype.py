# -*- coding: utf-8 -*-
"""
xcms datatypes

Author: lecorguille@sb-roscoff.fr
"""

import logging

from galaxy.datatypes.binary import RData

log = logging.getLogger(__name__)


class RdataCameraPositive(RData):
    file_ext = "rdata.camera.positive"
    blurb = "Output of CAMERA.annotate using positif acquisition mode\n\
        Input of ProbMetab or CAMERA.combinexcAnnos"

    # TODO: sniff


class RdataCameraNegative(RData):
    file_ext = "rdata.camera.negative"
    blurb = "Output of CAMERA.annotate using negatif acquisition mode\n\
        Input of ProbMetab or CAMERA.combinexcAnnos"

    # TODO: sniff


class RdataCameraQuick(RData):
    file_ext = "rdata.camera.quick"
    blurb = "Output of CAMERA.annotate using the quick mode"

    # TODO: sniff
