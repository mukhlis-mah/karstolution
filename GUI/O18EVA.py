import math, evaporation, cmodel_frac, constants
def O18EVA(tmax, TC, pCO2, pCO2cave, h, v, R18_hco_ini, R18_h2o_ini, R18v, HCOMIX, h2o_new,tt):
    
    # Sourcecode to develope the evolution of the isotopic ratio of the oxygen
    # compostion of the oxygen isotopes 16O and 18O as a function of
    # temperature TC, supersaturation (pCO2), relative humidity (h) and wind
    # velocity (v). %(08.12.2010/m)
    
    eva = evaporation.evaporation(TC, h, v)
    fracs = cmodel_frac.cmodel_frac(TC)
    TK = 273.15 + TC
    #Tau precipitation (s); t=d/a (according to Baker 98)
    alpha_p = (1.188e-011 * TC**3 - 1.29e-011 * TC**2 + 7.875e-009 * TC + 4.844e-008)
    #Tau buffering, after Dreybrodt and Scholz (2010)
    T = 125715.87302 - 16243.30688*TC + 1005.61111*TC**2 - 32.71852*TC**3 + 0.43333*TC**4
    
    #Concentrations and mol mass
    outputcave = constants.constants(TC, pCO2cave)      #Concentrations of the spezies in the solution, with respect to cave pCO2                            
    HCOCAVE = outputcave[2][2]/math.sqrt(0.8)    #HCO3- concentration, with respect to cave pCO2 (mol/l)

    h2o_ini = h2o_new                          #Mol mass of the water, with respect to the volume of a single box (mol)
    H20_ini = h2o_ini*18/1000.
    hco_ini = HCOMIX*H20_ini                   #Mol mass of the HCO3-(soil), with respect to the volume of a single box (mol)
    #hco_eq = HCOCAVE*H20_ini                   #Mol mass of the HCO3-(cave), with respect to the volume of a single box (mol)

    #Fractionation facors for oxygen isotope
    eps_m = fracs['a18_m'] - 1                                                 
    avl = ((-7356./TK + 15.38)/1000. + 1)
    abl = 1/(fracs['e18_hco_h2o'] + 1)
    a = 1/1.008*1.003
    f = 1/6

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    #Calculation of the 18R

    r_hco18 = [R18_hco_ini]
    r_h2o18 = [R18_h2o_ini]
    
    if tmax > math.floor(h2o_ini/eva):
        tmax = int(math.floor(h2o_ini/eva))
        print 'DRIPINTERVALL IS TOO LONG, THE WATERLAYER EVAPORATES COMPLETLY FOR THE GIVEN d, run # %i' %(tt)
    
    dt = 1
    t=range(1,tmax+1)
            
    HCO = [HCOMIX]                                  #Konzentration von HCO3-
    hco = [hco_ini]                                 #Menge an HCO3-
    H2O = [H20_ini]
    h2o = [h2o_new]

    #"Restwassermenge" und daher konstant
 
    HCO_EQ = HCOCAVE                       #Equilibriumconcentration 
    
    for i in t:

        delta = (H2O[-1]/1000.)/0.001
        
        #Neue Wassermenge
        h2o.append(h2o[-1] - eva*dt)              #Water (mol)
        d_h2o = -eva                              #Evaporationrate (mol/l)
        H2O.append(h2o[-1]*18*1e-3)               #Water (l)

        #Verdundstung
        HCO_temp = (HCO[-1] - HCO_EQ) * math.exp(-dt/(delta/alpha_p)) + HCO_EQ          #HCO3- concentration after timeintervall dt
        hco.append(HCO_temp * H2O[-2])                                         #HCO3- mass (mol)
                
        HCO.append(HCO_temp * (H2O[-2]/H2O[-1]))      #HCO3- concentration after timeintervall dt and the evaporation of water
                
        r_hco18.append(r_hco18[-1] + ((eps_m*(hco[-1]-hco[-2])/hco[-1]-1/T) * r_hco18[-1] + abl/T*r_h2o18[-1]) * dt)
        r_h2o18.append(r_h2o18[-1] + ((hco[-1]/h2o[-1]/T - f/abl/h2o[-1]*(hco[-1]-hco[-2]) * r_hco18[-1] + (d_h2o/h2o[-1]*(a*avl/(1-h)-1) - hco[-1]/h2o[-1]*abl/T) * r_h2o18[-1] - a*h/(1-h)*R18v/h2o[-1]*d_h2o) * dt))
        
    return [r_hco18[-1], r_h2o18[-1], HCO[-1], hco[-1], H2O[-1], h2o[-1]]

