#!/usr/bin/env python3

from isatools.convert import isatab2w4m
import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    isatab2w4m.main()
