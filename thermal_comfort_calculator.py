import math

def FNPS(TA):
    return math.exp(16.6536 - 4030.183 / (TA + 235))  # saturated vapor pressure in Pa

def calculate_PMV_PPD(CLO, MET, WME, TA, TR, VEL, RH=None, PA=None):
    if RH is None and PA is None:
        raise ValueError("Either relative humidity (RH) or water vapor pressure (PA) must be provided.")
    
    if PA is None:
        PA = RH * 10 * FNPS(TA)  # water vapor pressure, Pa
    
    ICL = 0.155 * CLO  # thermal insulation of clothing in m²K/W
    M = MET * 58.15  # metabolic rate in W/m²
    W = WME * 58.15  # external work in W/m²
    MW = M - W  # internal heat production in the human body
    
    if ICL < 0.078:
        FCL = 1 + 1.29 * ICL
    else:
        FCL = 1.05 + 0.645 * ICL  # clothing area factor
    
    HCF = 12.1 * math.sqrt(VEL)  # heat transfer coefficient by forced convection
    TAA = TA + 273  # air temperature in Kelvin
    TRA = TR + 273  # mean radiant temperature in Kelvin
    
    # Iterative calculation of surface temperature of clothing
    TCLA = TAA + (35.5 - TA) / (3.5 * (6.45 * ICL + 0.1))
    P1 = ICL * FCL
    P2 = P1 * 3.96
    P3 = P1 * 100
    P4 = P1 * TAA
    P5 = 308.7 - 0.028 * MW + P2 * (TRA / 100) ** 4
    
    XN = TCLA / 100
    XF = XN
    N = 0  # number of iterations
    EPS = 0.00015  # stop criteria in iteration
    
    while N <= 150:
        XF = (XF + XN) / 2
        HCN = 2.38 * abs(100 * XF - TAA) ** 0.25  # heat transfer coefficient by natural convection
        HC = HCF if HCF > HCN else HCN
        XN = (P5 + P4 * HC - P2 * XF ** 4) / (100 + P3 * HC)
        N += 1
        if abs(XN - XF) <= EPS:
            break
    
    TCL = 100 * XN - 273  # surface temperature of the clothing
    
    # Heat loss components
    HL1 = 3.05 * 0.001 * (5733 - 6.99 * MW - PA)  # heat loss diff. through skin
    HL2 = 0.42 * (MW - 58.15) if MW > 58.15 else 0  # heat loss by sweating (comfort)
    HL3 = 1.7 * 0.00001 * M * (5867 - PA)  # latent respiration heat loss
    HL4 = 0.0014 * M * (34 - TA)  # dry respiration heat loss
    HL5 = 3.96 * FCL * (XN ** 4 - (TRA / 100) ** 4)  # heat loss by radiation
    HL6 = FCL * HC * (TCL - TA)  # heat loss by convection
    
    # Predicted Mean Vote (PMV)
    TS = 0.303 * math.exp(-0.036 * M) + 0.028  # thermal sensation trans. Coeff.
    PMV = TS * (MW - HL1 - HL2 - HL3 - HL4 - HL5 - HL6)
    
    # Predicted Percentage Dissatisfied (PPD)
    PPD = 100 - 95 * math.exp(-0.03353 * PMV ** 4 - 0.2179 * PMV ** 2)
    
    return PMV, PPD

# Example usage
CLO = float(input("Clothing (clo): "))
MET = float(input("Metabolic rate (met): "))
WME = float(input("External work, normally around 0 (met): "))
TA = float(input("Air Temperature (°C): "))
TR = float(input("Mean radiant temperature (°C): "))
VEL = float(input("Relative air velocity (m/s): "))

RH_option = input("Do you have relative humidity (RH)? (Y/N): ")
if RH_option.upper() == "Y":
    RH = float(input("Relative humidity (%): "))
    PA = None
else:
    PA = float(input("Water vapor pressure (Pa): "))
    RH = None

PMV, PPD = calculate_PMV_PPD(CLO, MET, WME, TA, TR, VEL, RH, PA)
print("Predicted Mean Vote (PMV): {:.3f}".format(PMV))
print("Predicted Percentage of Dissatisfied (PPD): {:.3f}".format(PPD))
