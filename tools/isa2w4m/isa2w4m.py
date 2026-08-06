#!/usr/bin/env python3

import warnings

from isatools.convert import isatab2w4m

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    isatab2w4m.main()
