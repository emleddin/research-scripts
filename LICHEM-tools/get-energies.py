#!/usr/env/python3
## Get energy breakdown for QM and MM, and convert units as needed.
import re

## How do you want results given to you?
## Enter as a list, even if it's just 1 value! Ex1: ['kcal'] Ex2: ['au', 'ev']
## Options: 'au', 'hartree', 'hartrees', 'kcal', 'kcal/mol', 'ev', 'all'
## Note: 'au' and 'hartree(s)' are equivalent, so 'all' only prints 'au'!
# result_format = ['Hartree', 'eV', 'kcal/mol']
result_format = ['au']

## Did you use AMOEBA for the Electrostatics keyword?
## Currently not supported, because I don't have a test system to verify
##  polarizaton energy, E += TINKERPolEnergy
## AMOEBA = True, CHARGES = False
# AMOEBA = False

## ========================= Behind the Curtain ========================= ##

## Constants as implemented in LICHEM
avogNum = 6.02214129e23
har2ev = 27.21138505
SI2eV = 1/(1.602176565e-19)

kcal2ev = 4184*SI2eV/avogNum
au2kcal = 627.51

def ReadGaussian():
    """Parses the Gaussian log file for the QM energy.
    Returns
    -------
    GaussE : float
        The QM energy from Gaussian, in eV.
    """
    # LICHEM calculates this as the optimization energy - self energy.
    # Self energy of the charges = {f} a.u.
    self_line = " Self energy of the charges"
    # SCF Done:  E({s}) = {f} A.U. after {d} cycles
    SCF_line = " SCF Done:"
    with open('LICHM_GaussEnergy_0.log') as f:
        for line in f:
            if line.startswith(self_line):
                # print(line)
                ## Parses unmarked positive and - int and float values
                selfE = re.findall(r'\-\d+\.*\d*', line)
                selfE = float(selfE[0])
            elif line.startswith(SCF_line):
                # print(line)
                SCFE = re.findall(r'\-\d+\.*\d*', line)
                SCFE = float(SCFE[0])
    f.close()
    GaussE = SCFE - selfE
    GaussE *= har2ev
    return GaussE

# def ReadTinker(AMOEBA):
def ReadTinker():
    """Parses the TINKER log file for the MM energy.
    Currently does not support AMOEBA polarization terms!
    Returns
    -------
    TinkE : float
        The MM energy from Tinker, in eV.
    """
    #  Total Potential Energy : {f} Kcal/mole
    total_line = " Total Potential Energy :"
    with open('LICHM_TINKEREnergy_0.log') as f:
        for line in f:
            if line.startswith(total_line):
                # print(line)
                TinkE = re.findall(r'\-\d+\.*\d*', line)
                TinkE = float(TinkE[0])
                # if AMOEBA == True:
                #     if line.startswith("Polarization"):
                #         Epol = re.findall(r'\-\d+\.*\d*', line)
                #         Epol = float(Epol[0])
                #     elif line.startswith("Implicit Solvation")
                #         Esolv = re.findall(r'\-\d+\.*\d*', line)
                #         Esolv = float(Esolv[0])
    f.close()
    # if AMOEBA == True:
    #     TINKERPolForces = EPol + ESolv
    #     TinkE += TINKERPolForces
    #
    TinkE *= kcal2ev
    return TinkE

def GetSum(GaussE, TinkE):
    """Sums the QM and MM energy.
    Returns
    -------
    sumE : float
        The combined QM and MM energy, in eV.
    """
    sumE = GaussE + TinkE
    return sumE

def PrintDetailedE(GaussE, TinkE, sumE, result_format):
    """Prints out the energy breakdown as requested.
    Formats them to 10 decimal places.
    """
    ## Force to lowercase for search
    result_format = [each_form.lower() for each_form in result_format]
    if 'all' in result_format:
        print("\n")
        print(f"Gaussian Energy (au): {GaussE/har2ev:.10f} +")
        print(f"Tinker Energy   (au): {TinkE/har2ev:.10f}")
        print("===========================================")
        print(f"QM/MM Energy    (au): {sumE/har2ev:.10f}\n")
        print("\n")
        print(f"Gaussian Energy (kcal/mol): {GaussE/kcal2ev:.10f} +")
        print(f"Tinker Energy   (kcal/mol): {TinkE/kcal2ev:.10f}")
        print("====================================================")
        print(f"QM/MM Energy    (kcal/mol): {sumE/kcal2ev:.10f}\n")
        print("\n")
        print(f"Gaussian Energy (eV): {GaussE:.10f} +")
        print(f"Tinker Energy   (eV): {TinkE:.10f}")
        print("============================================")
        print(f"QM/MM Energy    (eV): {sumE:.10f}\n")
    else:
        one_printed = False
        print("")
        for form in result_format:
            if form.lower() == 'au':
                    print(f"Gaussian Energy (au): {GaussE/har2ev:.10f} +")
                    print(f"Tinker Energy   (au): {TinkE/har2ev:.10f}")
                    print("==========================================")
                    print(f"QM/MM Energy    (au): {sumE/har2ev:.10f}\n")
                    one_printed = True
            elif form.lower() in ('hartree', 'hartrees', 'hartree(s)'):
                    print(f"Gaussian Energy (Hartrees): {GaussE/har2ev:.10f} +")
                    print(f"Tinker Energy   (Hartrees): {TinkE/har2ev:.10f}")
                    print("================================================")
                    print(f"QM/MM Energy    (Hartrees): {sumE/har2ev:.10f}\n")
                    one_printed = True
            elif form.lower() in ('kcal', 'kcal/mol'):
                    print(f"Gaussian Energy (kcal/mol): {GaussE/kcal2ev:.10f} +")
                    print(f"Tinker Energy   (kcal/mol): {TinkE/kcal2ev:.10f}")
                    print("====================================================")
                    print(f"QM/MM Energy    (kcal/mol): {sumE/kcal2ev:.10f}\n")
                    one_printed = True
            elif form.lower() == 'ev':
                    print(f"Gaussian Energy (eV): {GaussE:.10f} +")
                    print(f"Tinker Energy   (eV): {TinkE:.10f}")
                    print("============================================")
                    print(f"QM/MM Energy    (eV): {sumE:.10f}\n")
                    one_printed = True
            else:
                if one_printed == True:
                    print(f"I can't parse your requested result_format, '{form}'.")
                    print("I only accept a list of values, such as:")
                    print(" 'au', 'hartrees', 'ev', 'kcal', 'all'")
                else:
                    print(f"I can't parse your requested result_format, '{form}'.")
                    print("I only accept a list of values, such as:")
                    print(" 'au', 'hartrees', 'ev', 'kcal', 'all'")
                    print("Printing au instead...\n")
                    print(f"Gaussian Energy (au): {GaussE/har2ev:.10f} +")
                    print(f"Tinker Energy   (au): {TinkE/har2ev:.10f}")
                    print("============================================")
                    print(f"QM/MM Energy    (au): {sumE/har2ev:.10f}\n")

GaussE = ReadGaussian()
# TinkE = ReadTinker(AMOEBA)
TinkE = ReadTinker()
sumE = GetSum(GaussE, TinkE)
PrintDetailedE(GaussE, TinkE, sumE, result_format)
