import math
def cmodel_frac(TC):
    TK = 273.15 + TC

    e18_c4 = math.exp((- 2590000./(TK**2) - 1.89)/1000.)-1              #Beck 2005 %HCO3 -> H2O                 (1000ln alpha)
    e18_b5 = math.exp((2520000./(TK**2) + 12.12)/1000.)-1               #Beck 2005 %H2O -> CO_2**{aq}            (1000ln alpha) 
    e18_a1 = (-160515./(TK**2)+1441.76/TK - 1.9585)/1000.             #Thorstenson 2004 %CO_2**{g} -> CO_2**{aq}

    e18_s1 = ((e18_c4+1)*(e18_b5+1)/(e18_a1+1))-1                      #HCO3 -> CO_2{g}

    e18_d1 = math.exp((18030./TK - 32.42)/1000.)-1                           #Kim and O'Neil 1997 %H2O -> CaCO3      (1000ln alpha)
 
    e18_s2 = ((e18_c4+1)*(e18_d1+1))-1                                 #HCO3 -> CaCO3

    e18_m = 2./6 * e18_s1 + 3./6 * e18_s2 + 1./6 * e18_c4                 #mean enrichment for Rayleighfractionation HCO3- -> CO2 + CaCO3 + H2O
    a18_m = e18_m + 1                                                  #mean fractionationfacor for Rayleighfractionation HCO3- -> CO2 + CaCO3 + H2O

    #Output: fractionation and enrichment factors
    return {'e18_hco_caco':e18_s2, 'e18_hco_h2o': e18_c4, 'a18_m':a18_m}