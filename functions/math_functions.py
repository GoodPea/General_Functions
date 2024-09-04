#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 10:36:47 2024

@author: ethanschneider
"""
import numpy as np

def gauss(x, A, u, sigma):
    gauss_values = A*(1/(sigma*np.sqrt(2*np.pi))*
                      np.exp(-.5*(((x-u)/sigma)**2)))
    return gauss_values

def PolyCoefficients(x, coeffs):
    """ Returns a polynomial for ``x`` values for the ``coeffs`` provided.

    The coefficients must be in ascending order (``x**0`` to ``x**o``).
    """
    o = len(coeffs)
    print(f'# This is a polynomial of order {o}.')
    y = 0
    for i in range(o):
        y += coeffs[i]*x**i
    return y