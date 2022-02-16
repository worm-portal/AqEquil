import pandas as pd
import numpy as np
import math
import copy
import warnings
import rpy2.robjects as ro
from rpy2.robjects.packages import importr

def water(props, water_model, T=25, P=1):
    
    if not isinstance(props, list):
        props = [props]
    
    chnosz = importr("CHNOSZ")
    chnosz.water(water_model)
    return chnosz.water(ro.StrVector([str(p) for p in props]), T=T, P=P, P1=False)

def entropy(formula):
    chnosz = importr("CHNOSZ")
    return float(chnosz.entropy(formula))

def calc_logK(OBIGT_df, Tc, P, TP_i, water_model):
    
    OBIGT_TP, rows_added = calc_G_TP(OBIGT_df, Tc, P, water_model)
    
    dissrxn2logK_out = []
    for i in OBIGT_TP.index:
        dissrxn2logK_out.append(dissrxn2logK(OBIGT_TP, i, Tc))
    assert len(dissrxn2logK_out) == OBIGT_TP.shape[0]
    
    OBIGT_TP['dissrxn_logK_'+str(TP_i)] = dissrxn2logK_out
    
    # remove any rows added by calc_G_TP
    OBIGT_TP.drop(OBIGT_TP.tail(rows_added).index, inplace = True)
    
    return OBIGT_TP


def calc_G_TP(OBIGT, Tc, P, water_model):
    
    aq_out, H2O_Pt = hkf(property=["G"], parameters=OBIGT,
                         T=273.15+Tc, P=P, contrib=["n", "s", "o"],
                         H2O_props=["rho"], water_model=water_model)
    
    cgl_out = cgl(property=["G"], parameters=OBIGT, T=273.15+Tc, P=P)
    
    aq_col = pd.DataFrame.from_dict(aq_out, orient="index")
    cgl_col = pd.DataFrame.from_dict(cgl_out, orient="index")

    G_TP_df = pd.concat([aq_col, cgl_col], axis=1)
    G_TP_df.columns = ['aq','cgl']
    
    OBIGT["G_TP"] = G_TP_df['aq'].combine_first(G_TP_df['cgl'])
    
    rows_added = 0
    
    # add a row for water
    if "H2O" not in list(OBIGT["name"]):
        OBIGT = OBIGT.append(pd.DataFrame([[float('NaN')] * len(OBIGT.columns)], columns=OBIGT.columns), ignore_index=True)
        OBIGT.iloc[-1, OBIGT.columns.get_loc('name')] = "H2O"
        OBIGT.iloc[-1, OBIGT.columns.get_loc('tag')] = "nan"
        OBIGT.iloc[-1, OBIGT.columns.get_loc('G_TP')] = float(water("G", water_model, T=Tc+273.15, P=P)["G"])
        rows_added += 1

    # add a row for protons
    if "H+" not in list(OBIGT["name"]):
        OBIGT = OBIGT.append(pd.DataFrame([[float('NaN')] * len(OBIGT.columns)], columns=OBIGT.columns), ignore_index=True)
        OBIGT.iloc[-1, OBIGT.columns.get_loc('name')] = "H+"
        OBIGT.iloc[-1, OBIGT.columns.get_loc('tag')] = "nan"
        OBIGT.iloc[-1, OBIGT.columns.get_loc('G_TP')] = 0
        rows_added += 1
    
    return OBIGT, rows_added
    
    
def G2logK(G, Tc):
    # Gas constant R is in cal/mol K
    return G / (-math.log(10) * 1.9872 * (273.15+Tc))


def dissrxn2logK(OBIGT, i, Tc):
    
    this_dissrxn = OBIGT.iloc[i, OBIGT.columns.get_loc('dissrxn')]
    
    try:
        split_dissrxn = this_dissrxn.split(" ")
    except:
        return float('NaN')
    
    coeff = [float(n) for n in split_dissrxn[::2]]
    species = split_dissrxn[1::2]
    try:
        G = sum([float(c*OBIGT.loc[OBIGT["name"]==sp, "G_TP"]) for c,sp in zip(coeff, species)])
    except:
        G_list = []
        for ii, sp in enumerate(species):
            G_TP = OBIGT.loc[OBIGT["name"]==sp, "G_TP"]
            if len(G_TP) == 1:
                G_list.append(float(coeff[ii]*OBIGT.loc[OBIGT["name"]==sp, "G_TP"]))
            else:
                ### check valid polymorph T

                # get polymorph entries of OBIGT that match mineral
                poly_df = copy.copy(OBIGT.loc[OBIGT["name"]==sp,:])
                # ensure polymorph df is sorted according to cr, cr2, cr3... etc.
                poly_df = poly_df.sort_values("state")

                z_Ts = list(poly_df.loc[poly_df["name"]==sp, "z.T"])

                last_t = float('-inf')
                appended=False
                for iii,t in enumerate(z_Ts):

                    if Tc+273.15 > last_t and Tc+273.15 < t:
                        G_list.append(float(coeff[ii]*list(poly_df.loc[poly_df["name"]==sp, "G_TP"])[iii]))
                        appended=True
                    if not appended and z_Ts[-1] == t:
                        G_list.append(float(coeff[ii]*list(poly_df.loc[poly_df["name"]==sp, "G_TP"])[iii]))
                    last_t = t

        G = sum(G_list)

    return G2logK(G, Tc)



def convert_cm3bar(value):
    return value*4.184 * 10

# GB conversion notes:
# - T and P are supposed to be lists but can only be single value lists until pyCHNOSZ's water can be made to accept lists for T and P
def hkf(property=None, parameters=None, T=298.15, P=1,
    contrib = ["n", "s", "o"], H2O_props=["rho"], water_model="SUPCRT92"):
    # calculate G, H, S, Cp, V, kT, and/or E using
    # the revised HKF equations of state
    # H2O_props - H2O properties needed for subcrt() output
    # constants
    Tr = 298.15 # K
    Pr = 1      # bar
    Theta = 228 # K
    Psi = 2600  # bar
    
    # GB conversion notes: hkf function is not vectorized...
#     # make T and P equal length
#     if not P == "Psat":
#         if len(P) < len(T):
#             P = [P]*len(T)
#         if len(T) < len(P):
#             T = [T]*len(P)
    
    # GB conversion note: handle error messages later
#     # nonsolvation, solvation, and origination contribution
#     notcontrib <- ! contrib %in% c("n", "s", "o")
#     if(TRUE %in% notcontrib) stop(paste("contrib must be in c('n', 's', 'o); got", c2s(contrib[notcontrib])))
    
    # get water properties
    # rho - for subcrt() output and g function
    # Born functions and epsilon - for HKF calculations
    H2O_props += ["QBorn", "XBorn", "YBorn", "epsilon"]

    if water_model == "SUPCRT92":
      # using H2O92D.f from SUPCRT92: alpha, daldT, beta - for partial derivatives of omega (g function)
      H2O_props += ["alpha", "daldT", "beta"]
    
    elif water_model == "IAPWS95":
      # using IAPWS-95: NBorn, UBorn - for compressibility, expansibility
      H2O_props += ["NBorn", "UBorn"]
    
    elif water_model == "DEW":
      # using DEW model: get beta to calculate dgdP
      H2O_props += ["beta"]
    
    H2O_PrTr = water(H2O_props, water_model, T=Tr, P=Pr) # pyCHNOSZ's water function does not handle lists yet, hence splitting this into two steps
    H2O_PT = water(H2O_props, water_model, T=T, P=P) # pyCHNOSZ's water function does not handle lists yet, hence splitting this into two steps
    ZBorn = -1 / H2O_PT.loc["1", "epsilon"]
    ZBorn_PrTr = -1 / H2O_PrTr.loc["1", "epsilon"]
    
    # a class to store the result
    out_dict = {} # dictionary to store output
    
    for k in parameters.index:
        
        if parameters["state"][k] != "aq":
            out_dict[k] = {p:float('NaN') for p in property}
        else:
            sp = parameters["name"][k]

            # loop over each species
            PAR = copy.copy(parameters.loc[k, :])

            PAR["a1.a"] = copy.copy(PAR["a1.a"]*10**-1)
            PAR["a2.b"] = copy.copy(PAR["a2.b"]*10**2)
            PAR["a4.d"] = copy.copy(PAR["a4.d"]*10**4)
            PAR["c2.f"] = copy.copy(PAR["c2.f"]*10**4)
            PAR["omega.lambda"] = copy.copy(PAR["omega.lambda"]*10**5)

            # substitute Cp and V for missing EoS parameters
            # here we assume that the parameters are in the same position as in thermo()$OBIGT
            # we don't need this if we're just looking at solvation properties (Cp_s_var, V_s_var)

            # GB conversion note: this block checks various things about EOS parameters.
            # for now, just set hasEOS to True
            hasEOS = True # delete this once the following block is converted to python
    #         if "n" in contrib:
    #             # put the heat capacity in for c1 if both c1 and c2 are missing
    #             if all(is.na(PAR[, 18:19])):
    #                 PAR[, 18] = PAR["Cp"]
    #             # put the volume in for a1 if a1, a2, a3 and a4 are missing
    #             if all(is.na(PAR[, 14:17])):
    #                 PAR[, 14] = convert(PAR["V"], "calories")
    #             # test for availability of the EoS parameters
    #             hasEOS = any(!is.na(PAR[, 14:21]))
    #             # if at least one of the EoS parameters is available, zero out any NA's in the rest
    #             if hasEOS:
    #                 PAR[, 14:21][, is.na(PAR[, 14:21])] = 0

            # compute values of omega(P,T) from those of omega(Pr,Tr)
            # using g function etc. (Shock et al., 1992 and others)
            omega = PAR["omega.lambda"]  # omega_PrTr
            # its derivatives are zero unless the g function kicks in
            dwdP, dwdT, d2wdT2 = 0, 0, 0
            Z = PAR["z.T"]

            omega_PT = PAR["omega.lambda"]
            if Z != 0 and Z != "NA" and PAR["name"] != "H+":
                # compute derivatives of omega: g and f functions (Shock et al., 1992; Johnson et al., 1992)
                rhohat = H2O_PT["rho"]/1000  # just converting kg/m3 to g/cm3

                
                # temporarily filter out Python's warnings about dividing by zero, which is possible
                # with the equations in the gfunction
                # Possible complex output is acounted for in gfun().
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    g = gfun(rhohat, T-273.15, P, H2O_PT["alpha"], H2O_PT["daldT"], H2O_PT["beta"])


                # after SUPCRT92/reac92.f
                eta = 1.66027E5
                reref = (Z**2) / (omega/eta + Z/(3.082 + 0))
                re = reref + abs(Z) * g["g"]
                omega_PT = eta * (Z**2/re - Z/(3.082 + g["g"]))
                Z3 = abs(Z**3)/re**2 - Z/(3.082 + g["g"])**2
                Z4 = abs(Z**4)/re**3 - Z/(3.082 + g["g"])**3
                dwdP = (-eta * Z3 * g["dgdP"])
                dwdT = (-eta * Z3 * g["dgdT"])
                d2wdT2 = (2 * eta * Z4 * g["dgdT"]**2 - eta * Z3 * g["d2gdT2"])

            # loop over each property
            w = float('NaN')
            for i,PROP in enumerate(property) :

                # over nonsolvation, solvation, or origination contributions
                hkf_p = 0

                for icontrib in contrib :
                    # various contributions to the properties
                    if icontrib == "n":
                        # nonsolvation ghs equations
                        if PROP == "H":
                            p_c = PAR["c1.e"]*(T-Tr) - PAR["c2.f"]*(1/(T-Theta)-1/(Tr-Theta))
                            p_a = PAR["a1.a"]*(P-Pr) + PAR["a2.b"]*math.log((Psi+P)/(Psi+Pr)) + \
                              ((2*T-Theta)/(T-Theta)**2)*(PAR["a3.c"]*(P-Pr)+PAR["a4.d"]*math.log((Psi+P)/(Psi+Pr)))
                            p = p_c + p_a
                        elif PROP == "S":
                            p_c = PAR["c1.e"]*math.log(T/Tr) - \
                              (PAR["c2.f"]/Theta)*( 1/(T-Theta)-1/(Tr-Theta) + \
                              math.log( (Tr*(T-Theta))/(T*(Tr-Theta)) )/Theta )
                            p_a = (T-Theta)**(-2)*(PAR["a3.c"]*(P-Pr)+PAR["a4.d"]*math.log((Psi+P)/(Psi+Pr)))
                            p = p_c + p_a
                        elif PROP == "G":
                            p_c = -PAR["c1.e"]*(T*math.log(T/Tr)-T+Tr) - \
                              PAR["c2.f"]*( (1/(T-Theta)-1/(Tr-Theta))*((Theta-T)/Theta) - \
                              (T/Theta**2)*math.log((Tr*(T-Theta))/(T*(Tr-Theta))) )
                            p_a = PAR["a1.a"]*(P-Pr) + PAR["a2.b"]*math.log((Psi+P)/(Psi+Pr)) + \
                              (PAR["a3.c"]*(P-Pr) + PAR["a4.d"]*math.log((Psi+P)/(Psi+Pr)))/(T-Theta)
                            p = p_c + p_a
                            # at Tr,Pr, if the origination contribution is not NA, ensure the solvation contribution is 0, not NA
                            if not np.isnan(PAR["G"]):
                                if T==Tr and P==Pr:
                                    p = 0
                        # nonsolvation cp v kt e equations
                        elif PROP == "Cp":
                            p = PAR["c1.e"] + PAR["c2.f"] * ( T - Theta ) ** (-2)        
                        elif PROP == "V":
                            p = convert_cm3bar(PAR["a1.a"]) + \
                              convert_cm3bar(PAR["a2.b"]) / (Psi + P) + \
                              (convert_cm3bar(PAR["a3.c"]) + convert_cm3bar(PAR["a4.d"]) / (Psi + P)) / (T - Theta)
#                         elif PROP == "kT":
#                             p = (convert(PAR["a2.b"], "cm3bar") + \
#                               convert(PAR["a4.d"], "cm3bar") / (T - Theta)) * (Psi + P) ** (-2)
#                         elif PROP == "E":
#                             p = convert( - (PAR["a3.c"] + PAR["a4.d"] / convert((Psi + P), "calories")) * \
#                               (T - Theta) ** (-2), "cm3bar")
                        else:
                            print("BAD")

                    if icontrib == "s":
                        # solvation ghs equations
                        if PROP == "G":
                            p = -omega_PT*(ZBorn+1) + omega*(ZBorn_PrTr+1) + omega*H2O_PrTr.loc["1", "YBorn"]*(T-Tr)
                            # at Tr,Pr, if the origination contribution is not NA, ensure the solvation contribution is 0, not NA
                            if(np.isnan(PAR["G"])):
                                if T==Tr and P==Pr:
                                    p = 0
                        if PROP == "H": 
                            p = -omega_PT*(ZBorn+1) + omega_PT*T*H2O_PT.loc["1", "YBorn"] + T*(ZBorn+1)*dwdT + \
                                   omega*(ZBorn_PrTr+1) - omega*Tr*H2O_PrTr.loc["1", "YBorn"]
                        if PROP == "S":
                            p = omega_PT*H2O_PT.loc["1","YBorn"] + (ZBorn+1)*dwdT - omega*H2O_PrTr.loc["1","YBorn"]
                        # solvation cp v kt e equations
                        if PROP == "Cp":
                            p = omega_PT*T*H2O_PT.loc["1","XBorn"] + 2*T*H2O_PT.loc["1","YBorn"]*dwdT + T*(ZBorn+1)*d2wdT2
                        if PROP == "V":
                            p = -convert_cm3bar(omega_PT) * \
                                H2O_PT.loc["1", "QBorn"] + convert_cm3bar(dwdP) * (-ZBorn - 1)
                        # TODO: the partial derivatives of omega are not included here here for kt and e
                        # (to do it, see p. 820 of SOJ+92 ... but kt requires d2wdP2 which we don"t have yet)
                        if PROP == "kT":
                            p = convert_cm3bar(omega) * H2O_PT["NBorn"]
                        if PROP == "E":
                            p = -convert_cm3bar(omega) * H2O_PT["UBorn"]

                    if icontrib == "o":
                        # origination ghs equations
                        if PROP == "G":
                            p = PAR["G"] - PAR["S"] * (T-Tr)
                            # don"t inherit NA from PAR$S at Tr
                            if T == Tr:
                                p = PAR["G"]
                        elif PROP == "H":
                            p = PAR["H"]
                        elif PROP == "S":
                            p = PAR["S"]
                        # origination eos equations (Cp, V, kT, E): senseless
                        else:
                            p = 0 * T

                    # accumulate the contribution
                    hkf_p = hkf_p + p
                
                # species have to be numbered (k) instead of named because of name repeats in db (e.g., cr polymorphs)
                if i > 0:
                    out_dict[k][PROP] = hkf_p
                else:
                    out_dict[k] = {PROP:hkf_p}

    return(out_dict, H2O_PT)




def gfun(rhohat, Tc, P, alpha, daldT, beta):
    ## g and f functions for describing effective electrostatic radii of ions
    ## split from hkf() 20120123 jmd      
    ## based on equations in
    ## Shock EL, Oelkers EH, Johnson JW, Sverjensky DA, Helgeson HC, 1992
    ## Calculation of the Thermodynamic Properties of Aqueous Species at High Pressures 
    ## and Temperatures: Effective Electrostatic Radii, Dissociation Constants and 
    ## Standard Partial Molal Properties to 1000 degrees C and 5 kbar
    ## J. Chem. Soc. Faraday Trans., 88(6), 803-826  doi:10.1039/FT9928800803
    # rhohat - density of water in g/cm3
    # Tc - temperature in degrees Celsius
    # P - pressure in bars
    # start with an output list of zeros
    out0 = float(len(rhohat))
    out = dict(g=out0, dgdT=out0, d2gdT2=out0, dgdP=out0)
    
    rhohat = rhohat[0]
    alpha = alpha[0]
    daldT = daldT[0]
    beta = beta[0]
    
    # only rhohat less than 1 will give results other than zero
    if rhohat >= 1:
        return {"g":0, "dgdT":0, "d2gdT2":0, "dgdP":0}

    # eta in Eq. 1
    eta = 1.66027E5
    # Table 3
    ag1 = -2.037662
    ag2 = 5.747000E-3
    ag3 = -6.557892E-6
    bg1 = 6.107361
    bg2 = -1.074377E-2
    bg3 = 1.268348E-5
    # Eq. 25
    ag = ag1 + ag2 * Tc + ag3 * Tc ** 2
    # Eq. 26
    bg = bg1 + bg2 * Tc + bg3 * Tc ** 2
    # Eq. 24
    g = ag * (1 - rhohat) ** bg
    
    # Table 4
    af1 = 0.3666666E2
    af2 = -0.1504956E-9
    af3 = 0.5017997E-13
    
    # Eq. 33
    f = ( ((Tc - 155) / 300) ** 4.8 + af1 * ((Tc - 155) / 300) ** 16 ) * \
        ( af2 * (1000 - P) ** 3 + af3 * (1000 - P) ** 4 )
    
    # limits of the f function (region II of Fig. 6)
    if Tc > 155 and P < 1000 and Tc < 355:
        ifg = True
    else:
        ifg = False
    
    # Eq. 32
    if ifg and not isinstance(f, complex):
        g = g*ifg - f*ifg
    
    # at P > 6000 bar (in DEW calculations), g is zero 20170926
    if P > 6000:
        g = 0

    ## now we have g at P, T
    # put the results in their right place (where rhohat < 1)
    out["g"] = g
    
    ## the rest is to get its partial derivatives with pressure and temperature
    ## after Johnson et al., 1992
    # alpha - coefficient of isobaric expansivity (K^-1)
    # daldT - temperature derivative of coefficient of isobaric expansivity (K^-2)
    # beta - coefficient of isothermal compressibility (bar^-1)
    # if these are NULL or NA (for IAPWS-95 and DEW), we skip the calculation
    if np.isnan(alpha): alpha = float('NaN')
    if np.isnan(daldT): daldT = float('NaN')
    if np.isnan(beta): beta = float('NaN')
    # Eqn. 76
    d2fdT2 = (0.0608/300*((Tc-155)/300)**2.8 + af1/375*((Tc-155)/300)**14) * (af2*(1000-P)**3 + af3*(1000-P)**4)
    # Eqn. 75
    dfdT = (0.016*((Tc-155)/300)**3.8 + 16*af1/300*((Tc-155)/300)**15) * \
        (af2*(1000-P)**3 + af3*(1000-P)**4)
    # Eqn. 74
    dfdP = -(((Tc-155)/300)**4.8 + af1*((Tc-155)/300)**16) * \
        (3*af2*(1000-P)**2 + 4*af3*(1000-P)**3)
    d2bdT2 = 2 * bg3  # Eqn. 73
    d2adT2 = 2 * ag3  # Eqn. 72
    dbdT = bg2 + 2*bg3*Tc  # Eqn. 71
    dadT = ag2 + 2*ag3*Tc  # Eqn. 70
    
    if isinstance(d2fdT2, complex): d2fdT2 = float('NaN')
    if isinstance(dfdT, complex): dfdT = float('NaN')
    if isinstance(dfdP, complex): dfdP = float('NaN')
    if isinstance(d2bdT2, complex): d2bdT2 = float('NaN')
    if isinstance(d2adT2, complex): d2adT2 = float('NaN')
    if isinstance(dbdT, complex): dbdT = float('NaN')
    if isinstance(dadT, complex): dadT = float('NaN')
    
    if not np.isnan(alpha) and not np.isnan(daldT):
        # Eqn. 69
        dgadT = bg*rhohat*alpha*(1-rhohat)**(bg-1) + math.log(1-rhohat)*g/ag*dbdT  
        D = rhohat
        
        # transcribed from SUPCRT92/reac92.f
        dDdT = -D * alpha
        #dDdP = D * beta
        dDdTT = -D * (daldT - alpha**2)
        Db = (1-D)**bg
        dDbdT = -bg*(1-D)**(bg-1)*dDdT + math.log(1-D)*Db*dbdT
        dDbdTT = -(bg*(1-D)**(bg-1)*dDdTT + (1-D)**(bg-1)*dDdT*dbdT + \
            bg*dDdT*(-(bg-1)*(1-D)**(bg-2)*dDdT + math.log(1-D)*(1-D)**(bg-1)*dbdT)) + \
            math.log(1-D)*(1-D)**bg*d2bdT2 - (1-D)**bg*dbdT*dDdT/(1-D) + math.log(1-D)*dbdT*dDbdT
        d2gdT2 = ag*dDbdTT + 2*dDbdT*dadT + Db*d2adT2
        
        if ifg:
            d2gdT2 = d2gdT2 - d2fdT2
        
        dgdT = g/ag*dadT + ag*dgadT  # Eqn. 67
        
        if ifg:
            dgdT = dgdT - dfdT
        
        # phew! done with those derivatives
        out["dgdT"] = dgdT
        out["d2gdT2"] = d2gdT2

    if not np.isnan(beta) :
        dgdP = -bg*rhohat*beta*g*(1-rhohat)**-1  # Eqn. 66
        if ifg:
            dgdP = dgdP - dfdP
        out["dgdP"] = dgdP
    
    return(out)

# CHNOSZ/cgl.R
# calculate standard thermodynamic properties of non-aqueous species
# 20060729 jmd

def cgl(property = None, parameters = None, T = 298.15, P = 1):
    # calculate properties of crystalline, liquid (except H2O) and gas species
    Tr = 298.15
    Pr = 1
    # the number of T, P conditions
    ncond = 1 # max([len(T), len(P)])
    # initialize output dict
    out_dict = dict()
    # loop over each species
    
    for k in parameters.index:
        
        if parameters["state"][k] == "aq":
            out_dict[k] = {p:float('NaN') for p in property}
        else:
            
            # the parameters for *this* species
            PAR = copy.copy(parameters.loc[k])

            PAR["a2.b"] = copy.copy(PAR["a2.b"]*10**-3)
            PAR["a3.c"] = copy.copy(PAR["a3.c"]*10**5)
            PAR["c1.e"] = copy.copy(PAR["c1.e"]*10**-5)

            # GB: converted code won't handle Berman minerals to begin with
    #         if(all(is.na(PAR[9:21]))) {
    #             # use Berman equations (parameters not in thermo()$OBIGT)
    #             properties <- berman(PAR["name"], T=T, P=P, thisinfo=PAR)
    #             iprop <- match(property, colnames(properties))
    #             values <- properties[, iprop, drop=FALSE]
    #         } else {

            # in CHNOSZ, we have
            # 1 cm^3 bar --> convert(1, "calories") == 0.02390057 cal
            # but REAC92D.F in SUPCRT92 uses
            cm3bar_to_cal = 0.023901488 # cal
            # start with NA values
            values = dict()
            # a test for availability of heat capacity coefficients (a, b, c, d, e, f)
            # based on the column assignments in thermo()$OBIGT

            if any([True if not np.isnan(p) else False for p in list(PAR.iloc[13:19])]):
                # we have at least one of the heat capacity coefficients;
                # zero out any NA's in the rest (leave lambda and T of transition (columns 19-20) alone)
                PAR.iloc[13:19] = [0 if np.isnan(p) else p for p in list(PAR.iloc[13:19])]
                # calculate the heat capacity and its integrals
                Cp = PAR["a1.a"] + PAR["a2.b"]*T + PAR["a3.c"]*T**-2 + PAR["a4.d"]*T**-0.5 + PAR["c1.e"]*T**2 + PAR["c2.f"]*T**PAR["omega.lambda"]
                intCpdT = PAR["a1.a"]*(T - Tr) + PAR["a2.b"]*(T**2 - Tr**2)/2 + PAR["a3.c"]*(1/T - 1/Tr)/-1 + PAR["a4.d"]*(T**0.5 - Tr**0.5)/0.5 + PAR["c1.e"]*(T**3-Tr**3)/3
                intCpdlnT = PAR["a1.a"]*math.log(T / Tr) + PAR["a2.b"]*(T - Tr) + PAR["a3.c"]*(T**-2 - Tr**-2)/-2 + PAR["a4.d"]*(T**-0.5 - Tr**-0.5)/-0.5  + PAR["c1.e"]*(T**2 - Tr**2)/2

                # do we also have the lambda parameter (Cp term with adjustable exponent on T)?
                if not np.isnan(PAR["omega.lambda"]) and PAR["omega.lambda"] != 0:
                    # equations for lambda adapted from Helgeson et al., 1998 (doi:10.1016/S0016-7037(97)00219-6)
                    if PAR["omega.lambda"] == -1:
                        intCpdT = intCpdT + PAR["c2.f"]*log(T/Tr) 
                    else:
                        intCpdT = intCpdT - PAR["c2.f"]*( T**(PAR["omega.lambda"] + 1) - Tr**(PAR["omega.lambda"] + 1) ) / (PAR["omega.lambda"] + 1)
                    intCpdlnT = intCpdlnT + PAR["c2.f"]*(T**PAR["omega.lambda"] - Tr**PAR["omega.lambda"]) / PAR["omega.lambda"]

            else:
                # use constant heat capacity if the coefficients are not available
                Cp = PAR["Cp"]
                intCpdT = PAR["Cp"]*(T - Tr)
                intCpdlnT = PAR["Cp"]*math.log(T / Tr)
                # in case Cp is listed as NA, set the integrals to 0 at Tr
                if T == Tr:
                    intCpdT = 0
                    intCpdlnT = 0


            # volume and its integrals
            if PAR["name"] in ["quartz", "coesite"]:
                # volume calculations for quartz and coesite
                qtz = quartz_coesite(PAR, T, P)
                V = qtz["V"]
                intVdP = qtz["intVdP"]
                intdVdTdP = qtz["intdVdTdP"]

            else:
                # for other minerals, volume is constant (Helgeson et al., 1978)
                V = PAR["V"]
                # if the volume is NA, set its integrals to zero
                if np.isnan(PAR["V"]):
                    intVdP = 0
                    intdVdTdP = 0
                else:
                    intVdP = PAR["V"]*(P - Pr) * cm3bar_to_cal
                    intdVdTdP = 0

            # get the values of each of the requested thermodynamic properties
            for i,prop in enumerate(property):
                if prop == "Cp": values["Cp"] = Cp
                if prop == "V": values["V"] = V
                if prop == "E": values["E"] = float('NaN')
                if prop == "kT": values["kT"] = float('NaN')
                if prop == "G":
                    # calculate S * (T - Tr), but set it to 0 at Tr (in case S is NA)
                    Sterm = PAR["S"]*(T - Tr)
                    if T == Tr:
                        Sterm = 0

                    values["G"] = PAR["G"] - Sterm + intCpdT - T*intCpdlnT + intVdP
                if prop == "H":
                    values["H"] = PAR["H"] + intCpdT + intVdP - T*intdVdTdP
                if prop == "S": values["S"] = PAR["S"] + intCpdlnT - intdVdTdP

            out_dict[k] = values # species have to be numbered instead of named because of name repeats (e.g., cr polymorphs)

    return out_dict


### unexported function ###

# calculate GHS and V corrections for quartz and coesite 20170929
# (these are the only mineral phases for which SUPCRT92 uses an inconstant volume)
def quartz_coesite(PAR, T, P):
    # the corrections are 0 for anything other than quartz and coesite
    if not PAR["name"] in ["quartz", "coesite"]:
        return(dict(G=0, H=0, S=0, V=0))
    # Tr, Pr and TtPr (transition temperature at Pr)
    Pr = 1      # bar
    Tr = 298.15 # K
    TtPr = 848  # K
    # constants from SUP92D.f
    aa = 549.824
    ba = 0.65995
    ca = -0.4973e-4
    VPtTta = 23.348
    VPrTtb = 23.72
    Stran = 0.342
    # constants from REAC92D.f
    VPrTra = 22.688 # VPrTr(a-quartz)
    Vdiff = 2.047   # VPrTr(a-quartz) - VPrTr(coesite)
    k = 38.5       # dPdTtr(a/b-quartz)
    #k <- 38.45834    # calculated in CHNOSZ: dPdTtr(info("quartz"))
    # code adapted from REAC92D.f
    qphase = PAR["state"].replace("cr", "")
    
    if qphase == "2":
        Pstar = P
        Sstar = 0
        V = copy.copy(VPrTtb)
    else:
        Pstar = Pr + k * (T - TtPr)
        Sstar = copy.copy(Stran)
        V = VPrTra + ca*(P-Pr) + (VPtTta - VPrTra - ca*(P-Pr))*(T-Tr) / (TtPr + (P-Pr)/k - Tr)
    
    if T < TtPr:
        Pstar = Pr
        Sstar = 0

    if PAR["name"] == "coesite":
        VPrTra = VPrTra - Vdiff
        VPrTtb = VPrTtb - Vdiff
        V = V - Vdiff
    
    cm3bar_to_cal = 0.023901488
    GVterm = cm3bar_to_cal * (VPrTra * (P - Pstar) + VPrTtb * (Pstar - Pr) - \
        0.5 * ca * (2 * Pr * (P - Pstar) - (P**2 - Pstar**2)) - \
        ca * k * (T - Tr) * (P - Pstar) + \
        k * (ba + aa * ca * k) * (T - Tr) * math.log((aa + P/k) / (aa + Pstar/k)))
    SVterm = cm3bar_to_cal * (-k * (ba + aa * ca * k) * \
        math.log((aa + P/k) / (aa + Pstar/k)) + ca * k * (P - Pstar)) - Sstar
    
    # note the minus sign on "SVterm" in order that intdVdTdP has the correct sign
    return dict(intVdP=GVterm, intdVdTdP=-SVterm, V=V)


def OBIGT2eos(OBIGT, fixGHS=True, tocal=True):
    
    OBIGT_out = OBIGT.copy()
    for i in range(0, OBIGT.shape[0]):
        
        # we only convert column 20 for aqueous species (omega), not for cgl species (lambda)
        if tocal and OBIGT.iloc[i, :]["E_units"] == "J" and OBIGT.iloc[i, :]["state"] == "aq":
            OBIGT_out.iloc[i, 8:12] = OBIGT.iloc[i, 8:12]/4.184
            OBIGT_out.iloc[i, 13:20] = OBIGT.iloc[i, 13:20]/4.184
            OBIGT_out.iloc[i, OBIGT.columns.get_loc('E_units')] = "cal"
            
        elif tocal and OBIGT.iloc[i, :]["E_units"] == "J":
            OBIGT_out.iloc[i, 8:12] = OBIGT.iloc[i, 8:12]/4.184
            OBIGT_out.iloc[i, 13:19] = OBIGT.iloc[i, 13:19]/4.184
            OBIGT_out.iloc[i, OBIGT.columns.get_loc('E_units')] = "cal"
        
        # fill in one of missing G, H, S
        # for use esp. by subcrt because NA for one of G, H or S 
        # will preclude calculations at high T
        if fixGHS:
            # which entries are missing just one
            imiss = [np.isnan(v) for v in OBIGT.iloc[i, 8:11]]
            if sum(imiss) == 1:
                
                ii = np.where(imiss)[0][0]
                
                if OBIGT_out.columns[8+ii] == "G":
                    H = OBIGT_out.iloc[i, OBIGT_out.columns.get_loc('H')]
                    S = OBIGT_out.iloc[i, OBIGT_out.columns.get_loc('S')]
                    Selem = entropy(OBIGT_out.iloc[i, OBIGT_out.columns.get_loc('formula')])
                    T = 298.15
                    G = H - T*(S - Selem)
                    OBIGT_out.iloc[i, 8+ii] = G
                elif OBIGT_out.columns[8+ii] == "H":
                    G = OBIGT_out.iloc[i, OBIGT_out.columns.get_loc('G')]
                    S = OBIGT_out.iloc[i, OBIGT_out.columns.get_loc('S')]
                    Selem = entropy(OBIGT_out.iloc[i, OBIGT_out.columns.get_loc('formula')])
                    T = 298.15
                    H = G + T*(S - Selem)
                    OBIGT_out.iloc[i, 8+ii] = H
                elif OBIGT_out.columns[8+ii] == "S":
                    G = OBIGT_out.iloc[i, OBIGT_out.columns.get_loc('G')]
                    H = OBIGT_out.iloc[i, OBIGT_out.columns.get_loc('H')]
                    Selem = entropy(OBIGT_out.iloc[i, OBIGT_out.columns.get_loc('formula')])
                    T = 298.15
                    G = H - T*(S - Selem)
                    S = Selem + (G - H)/T
                    OBIGT_out.iloc[i, 8+ii] = S
    
    return OBIGT_out