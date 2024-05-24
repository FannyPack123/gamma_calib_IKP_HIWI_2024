# -*- coding: utf-8 -*-
"""
Created on Fri May 24 23:27:21 2024

@author: phili
"""

import io
import matplotlib.pyplot as plt
import argparse
import numpy as np
import os

from parse import quickParse
from parse import genSpectrum

quickParse(voltage="1300", src="Bi207", offset="0000")

genSpectrum(voltage="1300", src="Bi207", offset="0000")

print("hello world")