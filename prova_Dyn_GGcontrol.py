# -*- coding: utf-8 -*-
"""
Created on Sat Oct 29 13:23:43 2022

@author: stesi
"""

#% GEO dynamics With GG 

import daceypy_import_helper  # noqa: F401

from time import perf_counter
from typing import Callable

import numpy as np
from daceypy import DA, array


#adim 
ll = 1e5
tt = 1e5
R = 6378.137/ll 
mue = 398600*tt**2/ll**3
Rgeo = 42165.8/ll                              
om = np.sqrt(mue/Rgeo**3)


# values taken from 'vallado' 
c20 = -1.083e-3  
c21 = -2.414e-10  
c22 = 1.574e-6  
c30 = 2.532e-6  
c31 = 2.191e-6 
c32 = 3.089e-7  
c33 = 1.006e-7  
s21 = 1.543e-9  
s22 = -9.038e-7  
s31 = 2.687e-7 
s32 = -2.115e-7
s33 = 1.972e-7



def SPHmotionGG_CONTROL(dyn: array, t: float) -> array:   # defines x vector of states as a DA array by calling function array
    
    r     = dyn[0]
    l     = dyn[1]
    phi   = dyn[2]
    v     = dyn[3]
    csi   = dyn[4]
    n     = dyn[5]
    
    l_r   = dyn[6]
    l_l   = dyn[7]  
    l_phi = dyn[8]
    l_v   = dyn[9]
    l_csi = dyn[10]
    l_n   = dyn[11]
                        
    #PERTURBATIONS
    aPg_r = -mue*((3*R**2 * c20*(3*(phi.sin())**2 - 1))/(2*r**4) + (9*R**2 * (phi.cos())**2*(s22*(2*l).sin() + c22*(2*l).cos()))/r**4 + (60*R**3 * (phi.cos())**3 * (s33*(3*l).sin() + c33*(3*l).cos()))/r**5 + (2*R**3 * c30*phi.sin()*(5*(phi.sin())**2 - 3))/r**5 + (9*R**2*phi.cos()*phi.sin()*(c21*l.cos() + s21*l.sin()))/r**4 + (60*R**3 * (phi.cos())**2 * phi.sin()*(s32*(2*l).sin() + c32*(2*l).cos()))/r**5 + (2*R**3 * phi.cos()*(15*(phi.sin())**2 - 3)*(c31*l.cos() + s31*l.sin()))/r**5)
    aPg_l =  -mue*((3*R**2 * (phi.cos())**2 * (2*c22*(2*l).sin() - 2*s22*(2*l).cos()))/r**3 + (15*R**3 * (phi.cos())**3 * (3*c33*(3*l).sin() - 3*s33*(3*l).cos()))/r**4 + (3*R**2 * phi.cos()*phi.sin()*(c21*l.sin() - s21*l.cos()))/r**3 + (15*R**3 * (phi.cos())**2 * phi.sin()*(2*c32*(2*l).sin() - 2*s32*(2*l).cos()))/r**4 + (R**3 * phi.cos()*(15*(phi.sin())**2 - 3)*(c31*l.sin() - s31*l.cos()))/(2*r**4))
    aPg_phi = mue*((15*R**3 * (phi.cos())**3 * (s32*(2*l).sin() + c32*(2*l).cos()))/r**4 + (3*R**2 * (phi.cos())**2 * (c21*l.cos() + s21*l.sin()))/r**3 - (3*R**2 * (phi.sin())**2 * (c21*l.cos() + s21*l.sin()))/r**3 + (R**3 * c30*phi.cos()*(5*(phi.sin())**2 - 3))/(2*r**4) - (30*R**3 * phi.cos()*(phi.sin())**2 * (s32*(2*l).sin() + c32*(2*l).cos()))/r**4 - (45*R**3 * (phi.cos())**2 * phi.sin()*(s33*(3*l).sin() + c33*(3*l).cos()))/r**4 + (5*R**3 * c30*phi.cos()*(phi.sin())**2)/r**4 - (R**3 * phi.sin()*(15*(phi.sin())**2 - 3)*(c31*l.cos() + s31*l.sin()))/(2*r**4) + (15*R**3 * (phi.cos())**2 * phi.sin()*(c31*l.cos() + s31*l.sin()))/r**4 - (6*R**2 * phi.cos()*phi.sin()*(s22*(2*l).sin() + c22*(2*l).cos()))/r**3 + (3*R**2 * c20*phi.cos()*phi.sin())/r**3)
    
    
    Ddyn = array.identity(12)
    
    #DYNAMICS (CONTROLLED)
    Ddyn[0] = v
    Ddyn[1] = csi
    Ddyn[2] = n
    Ddyn[3] = -mue/r**2 + r*n**2 + r*(csi + om)**2 *(phi.cos())**2 + aPg_r - l_v 
    Ddyn[4] = 2*n*(csi + om)*phi.tan() - 2*v/r *(csi + om) + aPg_l/(r * phi.cos())**2 - l_csi/(r * phi.cos())**2 
    Ddyn[5] = -2*v/r *n - (csi + om)**2 *phi.sin()*phi.cos() + aPg_phi/r**2 - l_n/r**2 
        
    #COSTATE DYAMICS --- from matlab symbolic toolbox
    Ddyn[6] = - l_n*((2*l_n)/r**3 + (2*n*v)/r**2 - (2*mue*((15*R**3*phi.cos()**3*(s32*(2*l).sin() + c32*(2*l).cos()))/r**4 + (3*R**2*phi.cos()**2*(c21*l.cos() + s21*l.sin()))/r**3 - (3*R**2*phi.sin()**2*(c21*l.cos() + s21*l.sin()))/r**3 + (R**3*c30*phi.cos()*(5*phi.sin()**2 - 3))/(2*r**4) - (30*R**3*phi.cos()*phi.sin()**2*(s32*(2*l).sin() + c32*(2*l).cos()))/r**4 - (45*R**3*phi.cos()**2*phi.sin()*(s33*(3*l).sin() + c33*(3*l).cos()))/r**4 + (5*R**3*c30*phi.cos()*phi.sin()**2)/r**4 - (R**3*phi.sin()*(15*phi.sin()**2 - 3)*(c31*l.cos() + s31*l.sin()))/(2*r**4) + (15*R**3*phi.cos()**2*phi.sin()*(c31*l.cos() + s31*l.sin()))/r**4 - (6*R**2*phi.cos()*phi.sin()*(s22*(2*l).sin() + c22*(2*l).cos()))/r**3 + (3*R**2*c20*phi.cos()*phi.sin())/r**3))/r**3 - (mue*((60*R**3*phi.cos()**3*(s32*(2*l).sin() + c32*(2*l).cos()))/r**5 + (9*R**2*phi.cos()**2*(c21*l.cos() + s21*l.sin()))/r**4 - (9*R**2*phi.sin()**2*(c21*l.cos() + s21*l.sin()))/r**4 + (2*R**3*c30*phi.cos()*(5*phi.sin()**2 - 3))/r**5 - (120*R**3*phi.cos()*phi.sin()**2*(s32*(2*l).sin() + c32*(2*l).cos()))/r**5 - (180*R**3*phi.cos()**2*phi.sin()*(s33*(3*l).sin() + c33*(3*l).cos()))/r**5 + (20*R**3*c30*phi.cos()*phi.sin()**2)/r**5 - (2*R**3*phi.sin()*(15*phi.sin()**2 - 3)*(c31*l.cos() + s31*l.sin()))/r**5 + (60*R**3*phi.cos()**2*phi.sin()*(c31*l.cos() + s31*l.sin()))/r**5 - (18*R**2*phi.cos()*phi.sin()*(s22*(2*l).sin() + c22*(2*l).cos()))/r**4 + (9*R**2*c20*phi.cos()*phi.sin())/r**4))/r**2) - l_v*(phi.cos()**2*(csi + om)**2 + mue*((6*R**2*c20*(3*phi.sin()**2 - 1))/r**5 + (36*R**2*phi.cos()**2*(s22*(2*l).sin() + c22*(2*l).cos()))/r**5 + (300*R**3*phi.cos()**3*(s33*(3*l).sin() + c33*(3*l).cos()))/r**6 + (10*R**3*c30*phi.sin()*(5*phi.sin()**2 - 3))/r**6 + (36*R**2*phi.cos()*phi.sin()*(c21*l.cos() + s21*l.sin()))/r**5 + (300*R**3*phi.cos()**2*phi.sin()*(s32*(2*l).sin() + c32*(2*l).cos()))/r**6 + (10*R**3*phi.cos()*(15*phi.sin()**2 - 3)*(c31*l.cos() + s31*l.sin()))/r**6) + (2*mue)/r**3 + n**2) - l_csi*((2*l_csi)/(r**3*phi.cos()**2) + (2*v*(csi + om))/r**2 + (2*mue*((3*R**2*phi.cos()**2*(2*c22*(2*l).sin() - 2*s22*(2*l).cos()))/r**3 + (15*R**3*phi.cos()**3*(3*c33*(3*l).sin() - 3*s33*(3*l).cos()))/r**4 + (3*R**2*phi.cos()*phi.sin()*(c21*l.sin() - s21*l.cos()))/r**3 + (15*R**3*phi.cos()**2*phi.sin()*(2*c32*(2*l).sin() - 2*s32*(2*l).cos()))/r**4 + (R**3*phi.cos()*(15*phi.sin()**2 - 3)*(c31*l.sin() - s31*l.cos()))/(2*r**4)))/(r**3*phi.cos()**2) + (mue*((9*R**2*phi.cos()**2*(2*c22*(2*l).sin() - 2*s22*(2*l).cos()))/r**4 + (60*R**3*phi.cos()**3*(3*c33*(3*l).sin() - 3*s33*(3*l).cos()))/r**5 + (9*R**2*phi.cos()*phi.sin()*(c21*l.sin() - s21*l.cos()))/r**4 + (60*R**3*phi.cos()**2*phi.sin()*(2*c32*(2*l).sin() - 2*s32*(2*l).cos()))/r**5 + (2*R**3*phi.cos()*(15*phi.sin()**2 - 3)*(c31*l.sin() - s31*l.cos()))/r**5))/(r**2*phi.cos()**2))
    Ddyn[7] = (l_csi*mue*((3*R**2*phi.cos()**2*(4*s22*(2*l).sin() + 4*c22*(2*l).cos()))/r**3 + (15*R**3*phi.cos()**3*(9*s33*(3*l).sin() + 9*c33*(3*l).cos()))/r**4 + (3*R**2*phi.cos()*phi.sin()*(c21*l.cos() + s21*l.sin()))/r**3 + (15*R**3*phi.cos()**2*phi.sin()*(4*s32*(2*l).sin() + 4*c32*(2*l).cos()))/r**4 + (R**3*phi.cos()*(15*phi.sin()**2 - 3)*(c31*l.cos() + s31*l.sin()))/(2*r**4)))/(r**2*phi.cos()**2) - (l_n*mue*((3*R**2*phi.sin()**2*(c21*l.sin() - s21*l.cos()))/r**3 - (3*R**2*phi.cos()**2*(c21*l.sin() - s21*l.cos()))/r**3 - (15*R**3*phi.cos()**3*(2*c32*(2*l).sin() - 2*s32*(2*l).cos()))/r**4 + (30*R**3*phi.cos()*phi.sin()**2*(2*c32*(2*l).sin() - 2*s32*(2*l).cos()))/r**4 + (45*R**3*phi.cos()**2*phi.sin()*(3*c33*(3*l).sin() - 3*s33*(3*l).cos()))/r**4 + (R**3*phi.sin()*(15*phi.sin()**2 - 3)*(c31*l.sin() - s31*l.cos()))/(2*r**4) - (15*R**3*phi.cos()**2*phi.sin()*(c31*l.sin() - s31*l.cos()))/r**4 + (6*R**2*phi.cos()*phi.sin()*(2*c22*(2*l).sin() - 2*s22*(2*l).cos()))/r**3))/r**2 - l_v*mue*((9*R**2*phi.cos()**2*(2*c22*(2*l).sin() - 2*s22*(2*l).cos()))/r**4 + (60*R**3*phi.cos()**3*(3*c33*(3*l).sin() - 3*s33*(3*l).cos()))/r**5 + (9*R**2*phi.cos()*phi.sin()*(c21*l.sin() - s21*l.cos()))/r**4 + (60*R**3*phi.cos()**2*phi.sin()*(2*c32*(2*l).sin() - 2*s32*(2*l).cos()))/r**5 + (2*R**3*phi.cos()*(15*phi.sin()**2 - 3)*(c31*l.sin() - s31*l.cos()))/r**5)
    Ddyn[8] = l_n*(phi.cos()**2*(csi + om)**2 - phi.sin()**2*(csi + om)**2 + (mue*((6*R**2*phi.cos()**2*(s22*(2*l).sin() + c22*(2*l).cos()))/r**3 + (45*R**3*phi.cos()**3*(s33*(3*l).sin() + c33*(3*l).cos()))/r**4 - (3*R**2*c20*phi.cos()**2)/r**3 - (6*R**2*phi.sin()**2*(s22*(2*l).sin() + c22*(2*l).cos()))/r**3 - (30*R**3*phi.sin()**3*(s32*(2*l).sin() + c32*(2*l).cos()))/r**4 + (3*R**2*c20*phi.sin()**2)/r**3 + (5*R**3*c30*phi.sin()**3)/r**4 - (15*R**3*phi.cos()**3*(c31*l.cos() + s31*l.sin()))/r**4 + (R**3*c30*phi.sin()*(5*phi.sin()**2 - 3))/(2*r**4) + (12*R**2*phi.cos()*phi.sin()*(c21*l.cos() + s21*l.sin()))/r**3 + (105*R**3*phi.cos()**2*phi.sin()*(s32*(2*l).sin() + c32*(2*l).cos()))/r**4 - (90*R**3*phi.cos()*phi.sin()**2*(s33*(3*l).sin() + c33*(3*l).cos()))/r**4 + (R**3*phi.cos()*(15*phi.sin()**2 - 3)*(c31*l.cos() + s31*l.sin()))/(2*r**4) - (15*R**3*c30*phi.cos()**2*phi.sin())/r**4 + (45*R**3*phi.cos()*phi.sin()**2*(c31*l.cos() + s31*l.sin()))/r**4))/r**2) - l_csi*(2*n*(csi + om)*(phi.tan()**2 + 1) - (2*l_csi*phi.sin())/(r**2*phi.cos()**3) + (mue*((3*R**2*phi.sin()**2*(c21*l.sin() - s21*l.cos()))/r**3 - (3*R**2*phi.cos()**2*(c21*l.sin() - s21*l.cos()))/r**3 - (15*R**3*phi.cos()**3*(2*c32*(2*l).sin() - 2*s32*(2*l).cos()))/r**4 + (30*R**3*phi.cos()*phi.sin()**2*(2*c32*(2*l).sin() - 2*s32*(2*l).cos()))/r**4 + (45*R**3*phi.cos()**2*phi.sin()*(3*c33*(3*l).sin() - 3*s33*(3*l).cos()))/r**4 + (R**3*phi.sin()*(15*phi.sin()**2 - 3)*(c31*l.sin() - s31*l.cos()))/(2*r**4) - (15*R**3*phi.cos()**2*phi.sin()*(c31*l.sin() - s31*l.cos()))/r**4 + (6*R**2*phi.cos()*phi.sin()*(2*c22*(2*l).sin() - 2*s22*(2*l).cos()))/r**3))/(r**2*phi.cos()**2) - (2*mue*phi.sin()*((3*R**2*phi.cos()**2*(2*c22*(2*l).sin() - 2*s22*(2*l).cos()))/r**3 + (15*R**3*phi.cos()**3*(3*c33*(3*l).sin() - 3*s33*(3*l).cos()))/r**4 + (3*R**2*phi.cos()*phi.sin()*(c21*l.sin() - s21*l.cos()))/r**3 + (15*R**3*phi.cos()**2*phi.sin()*(2*c32*(2*l).sin() - 2*s32*(2*l).cos()))/r**4 + (R**3*phi.cos()*(15*phi.sin()**2 - 3)*(c31*l.sin() - s31*l.cos()))/(2*r**4)))/(r**2*phi.cos()**3)) + l_v*(mue*((60*R**3*phi.cos()**3*(s32*(2*l).sin() + c32*(2*l).cos()))/r**5 + (9*R**2*phi.cos()**2*(c21*l.cos() + s21*l.sin()))/r**4 - (9*R**2*phi.sin()**2*(c21*l.cos() + s21*l.sin()))/r**4 + (2*R**3*c30*phi.cos()*(5*phi.sin()**2 - 3))/r**5 - (120*R**3*phi.cos()*phi.sin()**2*(s32*(2*l).sin() + c32*(2*l).cos()))/r**5 - (180*R**3*phi.cos()**2*phi.sin()*(s33*(3*l).sin() + c33*(3*l).cos()))/r**5 + (20*R**3*c30*phi.cos()*phi.sin()**2)/r**5 - (2*R**3*phi.sin()*(15*phi.sin()**2 - 3)*(c31*l.cos() + s31*l.sin()))/r**5 + (60*R**3*phi.cos()**2*phi.sin()*(c31*l.cos() + s31*l.sin()))/r**5 - (18*R**2*phi.cos()*phi.sin()*(s22*(2*l).sin() + c22*(2*l).cos()))/r**4 + (9*R**2*c20*phi.cos()*phi.sin())/r**4) + 2*r*phi.cos()*phi.sin()*(csi + om)**2)
    Ddyn[9] = (2*l_n*n)/r - l_r + (2*l_csi*(csi + om))/r
    Ddyn[10] = l_csi*((2*v)/r - 2*n*phi.tan()) - l_l - l_v*r*phi.cos()**2*(2*csi + 2*om) + l_n*phi.cos()*phi.sin()*(2*csi + 2*om)
    Ddyn[11] = (2*l_n*v)/r - l_phi - 2*l_v*n*r - 2*l_csi*phi.tan()*(csi + om)
 

    return Ddyn


def RK78(Y0: array, X0: float, X1: float, f: Callable[[array, float], array]):

    Y0 = Y0.copy()

    N = len(Y0)

    H0 = 0.001
    HS = 0.1
    H1 = 100.0
    EPS = 1.e-12
    BS = 20 * EPS

    Z = array.zeros((N, 16))
    Y1 = array.zeros(N)

    VIHMAX = 0.0

    HSQR = 1.0 / 9.0
    A = np.zeros(13)
    B = np.zeros((13, 12))
    C = np.zeros(13)
    D = np.zeros(13)

    A = np.array([
        0.0, 1.0/18.0, 1.0/12.0, 1.0/8.0, 5.0/16.0, 3.0/8.0, 59.0/400.0,
        93.0/200.0, 5490023248.0/9719169821.0, 13.0/20.0,
        1201146811.0/1299019798.0, 1.0, 1.0,
    ])

    B[1, 0] = 1.0/18.0
    B[2, 0] = 1.0/48.0
    B[2, 1] = 1.0/16.0
    B[3, 0] = 1.0/32.0
    B[3, 2] = 3.0/32.0
    B[4, 0] = 5.0/16.0
    B[4, 2] = -75.0/64.0
    B[4, 3] = 75.0/64.0
    B[5, 0] = 3.0/80.0
    B[5, 3] = 3.0/16.0
    B[5, 4] = 3.0/20.0
    B[6, 0] = 29443841.0/614563906.0
    B[6, 3] = 77736538.0/692538347.0
    B[6, 4] = -28693883.0/1125000000.0
    B[6, 5] = 23124283.0/1800000000.0
    B[7, 0] = 16016141.0/946692911.0
    B[7, 3] = 61564180.0/158732637.0
    B[7, 4] = 22789713.0/633445777.0
    B[7, 5] = 545815736.0/2771057229.0
    B[7, 6] = -180193667.0/1043307555.0
    B[8, 0] = 39632708.0/573591083.0
    B[8, 3] = -433636366.0/683701615.0
    B[8, 4] = -421739975.0/2616292301.0
    B[8, 5] = 100302831.0/723423059.0
    B[8, 6] = 790204164.0/839813087.0
    B[8, 7] = 800635310.0/3783071287.0
    B[9, 0] = 246121993.0/1340847787.0
    B[9, 3] = -37695042795.0/15268766246.0
    B[9, 4] = -309121744.0/1061227803.0
    B[9, 5] = -12992083.0/490766935.0
    B[9, 6] = 6005943493.0/2108947869.0
    B[9, 7] = 393006217.0/1396673457.0
    B[9, 8] = 123872331.0/1001029789.0
    B[10, 0] = -1028468189.0/846180014.0
    B[10, 3] = 8478235783.0/508512852.0
    B[10, 4] = 1311729495.0/1432422823.0
    B[10, 5] = -10304129995.0/1701304382.0
    B[10, 6] = -48777925059.0/3047939560.0
    B[10, 7] = 15336726248.0/1032824649.0
    B[10, 8] = -45442868181.0/3398467696.0
    B[10, 9] = 3065993473.0/597172653.0
    B[11, 0] = 185892177.0/718116043.0
    B[11, 3] = -3185094517.0/667107341.0
    B[11, 4] = -477755414.0/1098053517.0
    B[11, 5] = -703635378.0/230739211.0
    B[11, 6] = 5731566787.0/1027545527.0
    B[11, 7] = 5232866602.0/850066563.0
    B[11, 8] = -4093664535.0/808688257.0
    B[11, 9] = 3962137247.0/1805957418.0
    B[11, 10] = 65686358.0/487910083.0
    B[12, 0] = 403863854.0/491063109.0
    B[12, 3] = - 5068492393.0/434740067.0
    B[12, 4] = -411421997.0/543043805.0
    B[12, 5] = 652783627.0/914296604.0
    B[12, 6] = 11173962825.0/925320556.0
    B[12, 7] = -13158990841.0/6184727034.0
    B[12, 8] = 3936647629.0/1978049680.0
    B[12, 9] = -160528059.0/685178525.0
    B[12, 10] = 248638103.0/1413531060.0

    C = np.array([
        14005451.0/335480064.0, 0.0, 0.0, 0.0, 0.0, -59238493.0/1068277825.0,
        181606767.0/758867731.0, 561292985.0/797845732.0,
        -1041891430.0/1371343529.0, 760417239.0/1151165299.0,
        118820643.0/751138087.0, -528747749.0/2220607170.0, 1.0/4.0,
    ])

    D = np.array([
        13451932.0/455176623.0, 0.0, 0.0, 0.0, 0.0, -808719846.0/976000145.0,
        1757004468.0/5645159321.0, 656045339.0/265891186.0,
        -3867574721.0/1518517206.0, 465885868.0/322736535.0,
        53011238.0/667516719.0, 2.0/45.0, 0.0,
    ])

    Z[:, 0] = Y0

    H = abs(HS)
    HH0 = abs(H0)
    HH1 = abs(H1)
    X = X0
    RFNORM = 0.0
    ERREST = 0.0

    while X != X1:

        # compute new stepsize
        if RFNORM != 0:
            H = H * min(4.0, np.exp(HSQR * np.log(EPS / RFNORM)))
        if abs(H) > abs(HH1):
            H = HH1
        elif abs(H) < abs(HH0) * 0.99:
            H = HH0
            print("--- WARNING, MINIMUM STEPSIZE REACHED IN RK")

        if (X + H - X1) * H > 0:
            H = X1 - X

        for j in range(13):

            for i in range(N):

                Y0[i] = 0.0
                # EVALUATE RHS AT 13 POINTS
                for k in range(j):
                    Y0[i] = Y0[i] + Z[i, k + 3] * B[j, k]

                Y0[i] = H * Y0[i] + Z[i, 0]

            Y1 = f(Y0, X + H * A[j])

            for i in range(N):
                Z[i, j + 3] = Y1[i]

        for i in range(N):

            Z[i, 1] = 0.0
            Z[i, 2] = 0.0
            # EXECUTE 7TH,8TH ORDER STEPS
            for j in range(13):
                Z[i, 1] = Z[i, 1] + Z[i, j + 3] * D[j]
                Z[i, 2] = Z[i, 2] + Z[i, j + 3] * C[j]

            Y1[i] = (Z[i, 2] - Z[i, 1]) * H
            Z[i, 2] = Z[i, 2] * H + Z[i, 0]

        Y1cons = Y1.cons()

        # ESTIMATE ERROR AND DECIDE ABOUT BACKSTEP
        RFNORM = np.linalg.norm(Y1cons, np.inf)  # type: ignore
        if RFNORM > BS and abs(H / H0) > 1.2:
            H = H / 3.0
            RFNORM = 0
        else:
            for i in range(N):
                Z[i, 0] = Z[i, 2]
            X = X + H
            VIHMAX = max(VIHMAX, H)
            ERREST = ERREST + RFNORM

    Y1 = Z[:, 0]

    return Y1



# MAIN(): 

exp_order = 4
DA.init(exp_order, 12)

# Set initial conditions

dyn0 = array.identity(12)       # initialize state + costate augmented vector 
dyn0[0]  += 42165.8/ll
dyn0[1]  += 60*np.pi/180
dyn0[2]  += 0
dyn0[3]  += 0
dyn0[4]  += 0
dyn0[5]  += 0

dyn0[6]  += 0
dyn0[7]  += 0
dyn0[8]  += 0
dyn0[9]  += 0
dyn0[10] += 0
dyn0[11] += 0


#Set control/free-drift time
hours = 3600                               
days = 24*hours                                
time = 2.5*days/tt       #SET 0.5 days for checking the controlled case                                     



# Compute final state given perturbed initial state 
T0 = perf_counter()       # time RquiRd for integration.. elapsed time 
with DA.cache_manager():  # optional, for efficiency
       Xf = RK78(dyn0, 0.0, time, SPHmotionGG_CONTROL)
T1 = perf_counter()

#print(str(Xf))

# DIMENSIONALIZE
# Xf[0] = Xf[0]*ll
# Xf[4] = Xf[4]/tt
# Xf[5] = Xf[5]/tt

print(f"Info: time RquiRd for integration = {T1 - T0} s")

deltaXi = np.array([0, -0.04*np.pi/180, 0, 
                    0,               0, 0, 
                    0,               0, 0, 0, 0, 0])  #NON CONTROLLED CASE

#from matlab..end of free drift conditions -- Perturbation on ref cond were POL must be evaluated 
# deltaXi = np.array([42162.5027715973/ll - dyn0[0].cons(), 1.04771524217833 - dyn0[1].cons(), -6.11361154669619e-08, 
                    # 3.02803748412535e-06, 1.12365365870834e-08*tt, 9.63339629885044e-14*tt, 
                   # -2e-12, -1e-12, 0, 2e-7, 1e-7, 0])  


# final_cond = Xf.eval(deltaXi)
print(f"Final conditions:\n{Xf.eval(deltaXi)}\n")   
