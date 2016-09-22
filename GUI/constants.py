import math
def constants(TC, pCO2):
    #calculates activity coefficients, mass action constants,
    #reaction rate constants and concentrations of all species [mol/l] comprised in the
    #CO2-H2O-CaCO3 system.
    
    TK = TC + 273.16
    if pCO2<=0:
        pCO2=0.0000000000001
    #temperature (only) dependent variables

    #reaction rate constants
    k1m = 10**(13.558 - 3617.1/TK)
    k1p = 10**(329.850 - 110.54 * math.log(TK,10) - 17265.4/TK)
    k2m = 10**(14.09 - 5308/TK)
    k2p = 10**(13.635 - 2985/TK)
    
    #mass action constants
    K2 = 10**(-107.8871 - 0.03252849 * TK + 5151.79 / TK + 38.92561 * math.log(TK,10) - 563713.9 / (TK**2))
    K5 = 1.707e-4
    K6 = 10**(-356.3094 + 21834.37 / TK - 0.060919964 * TK + 126.8339 * math.log(TK,10) - 1684915 / (TK**2))   

    K0 = K5 / K6                              #BUHMANN1985
    K1 = 10**(-356.3094 - 0.06091964 * TK + 21834.37 / TK + 126.8339 * math.log(TK,10) - 1684915 / (TK**2)) #KAUFMANN2003
    
    KC = 10**(-171.9065 - 0.077993 * TK + 2839.319 / TK + 71.595 * math.log(TK,10))
    KH = 10**(108.3865 + 0.01985076 * TK - 6919.53 / TK - 40.45154 * math.log(TK,10) + 669365 / (TK**2))
    KS = 10**-(8.15087602 + 0.0136633623 * TC - 3.5812701e-5 * TC**2)
    KW = 10**(22.801 - 4787.3 / TK - 0.010365 * TK - 7.1321 * math.log(TK,10))

    #constants for activity coefficients
    A = 0.48809 + 8.074e-4 * TC
    B = 0.3241 + 1.600e-4 * TC

    #calculation of the concentrations in equilibrium at a given pCO2

    Ca = 1e-3 #any starting concentration of calcium
    
    #calculation of the activity coefficients acoording to the boundary
    #condition IS = 3 * Ca. Equilibrium values establish after a few runs
    #(less than 10)
    
    gammaCO2 = 1
    
    for i in range(1,11):
        
        #ionic strength
        IS = 3 * Ca
        
        #activity coefficients
        gammaH = 10**(-A * math.sqrt(IS) / (1 + B * 9 * math.sqrt(IS)))       #Nach Dreybrodt1988/S15
        gammaCa = 10**(-A * 4 * math.sqrt(IS) / (1 + B * 5 * math.sqrt(IS)))  #Nach Dreybrodt1988/S15
        gammaHCO = 10**(-A * math.sqrt(IS) / (1 + B * 5.4 * math.sqrt(IS)))   #Nach Dreybrodt1988/S15
        gammaOH = 10**(-A * math.sqrt(IS) / (1 + B * 3.5 * math.sqrt(IS)))
        gammaCO3 = 10**(-A * 4 * math.sqrt(IS) / (1 + B * 5.4 * math.sqrt(IS)))
        
        #calcium
        Ca = (pCO2 * K1 * KC * KH / (4 * K2 * gammaCa * gammaHCO**2))**(1.0/3)
    
    #carbondioxide
    CO2 = KH * pCO2
    
    #hydrogen
    H = (K1**2 * K2 * KH**2 * gammaCa * pCO2**2 / (2 * KC * gammaHCO))**(1.0/3) / (gammaH)

    #bicarbonate
    HCO = (pCO2 * KH * K1) / (gammaHCO * gammaH * H)
    
    #pH value
    pH = -math.log(H,10)
    
    #OH
    OH = KW / (gammaOH * gammaH * H)

    #more reaction rate constants
    km = k1m * gammaH * gammaHCO * H / ( K1 * (1 + K0)) + k2m
    kp = k1p + k2p * OH
    
    kHCO = k2m + k1m * H

    #y parameter
    C1 = kp * KH * pCO2 / km * gammaCO2
    C2 = km
    
    ac = [gammaCa, gammaCO3, gammaH, gammaHCO, gammaOH]
    rc = [k1m, k1p, k2m, k2p, km, kp, kHCO]
    cc = [Ca, H, HCO, OH, CO2]
    mp = [C1, C2]
    ks = [K0, K1, K2, KH]
    
    return [ac, rc, cc, mp, pH, ks]