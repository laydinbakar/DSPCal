#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#   Author        : Levent Aydinbakar

from collections import Counter
import argparse
import numpy
import math
import streamlit as st
import pandas as pd
import sys
import re
from st_btn_select import st_btn_select

pd.options.display.float_format = '{:.10f}'.format

st.write("""
# Differential Scattering Parameters Calculator (DSPCal)
DSPCal web app has been developed in the following journal article. To use this web app, we would like you to agree to properly cite the following journal article in any resulting publications or presentations. 
""")

citation='M. Buyukyildiz, L. Aydinbakar, "DSPCal app", Journal name, Numbers'
st.text_area("", value=citation, height=3)

cite=st.checkbox("I understand and agree that any use of DSPCal web app shall be subject to the terms and conditions set forth in the paper above, and I agree to properly cite the paper in any resulting publications or presentations.")



st.sidebar.header("Chemical formula or composition")
chemical = st.sidebar.selectbox("Which one you would like to compute DSP for?", ("A chemical formula", "A chemical composition"))

if chemical == "A chemical formula":
  formula = str(st.sidebar.text_input("Write the chemical formula.", "H2O"))
elif chemical == "A chemical composition":
  formula = str(st.sidebar.text_input("The elements in the composition (e.g. HO, FeO, CHNOS.) ", "HO"))
#######################################################################################################
# This is taken from: https://codereview.stackexchange.com/questions/232630/parsing-molecular-formula
def parse_molecule(molecule):
    array = [[]]
    for token in re.findall('[A-Z][a-z]?|\d+|.', molecule):
        if token.isalpha() and token.istitle():
            last = [token]
            upper = token
            array[-1].append(token)
        elif token.isalpha():
            last = upper + token
            array[-1] = [last]
        elif token.isdecimal():
            count = int(token)
            array[-1].extend(last*(int(token)-1))
        elif token in '([{': # or token == '[':
            array.append([])
        elif token in ')]}': #or token == ']':
            last = array.pop()
            array[-1].extend(last)
        else:
          raise ValueError("Unrecognized character in molecule")
    return dict(Counter(array[-1]))
#######################################################################################################

atomNumbers = numpy.transpose(numpy.array(numpy.loadtxt("data/A.csv", delimiter="\t", dtype=str)))

elName=[]
elRate=[]
if chemical == "A chemical formula":
  for key, value in parse_molecule(formula).items():
    elName.append(key)
    elRate.append(float(value))
if chemical == "A chemical composition":
  composition = []
  for key, value in parse_molecule(formula).items():
    elName.append(key)
    composition.append(float(st.sidebar.text_input("Percentage of the element **"+str(key)+"** in the composition", 50))/100)
    #composition.append(float(st.sidebar.text_input("Percent of the element **"+str(key)+"** in the composition", 50))/100)

elAtomNumber=[]
for j in elName:
  for i in range(len(atomNumbers)):
    if atomNumbers[i,0] == j:
      elAtomNumber.append(float(atomNumbers[i,1]))



Fs = numpy.array(numpy.loadtxt("data/F.csv", delimiter="\t", dtype=str))
Ss = numpy.array(numpy.loadtxt("data/S.csv", delimiter="\t", dtype=str))
X = numpy.array(numpy.delete(Fs[:,0],0), dtype=float)

st.sidebar.write("---")
st.sidebar.header("Energy and angle")
option = st.sidebar.selectbox("Which one you would like to use?",("An energy value and an angle value", "An energy value and multiple angle values", "Multiple energy values and an angle value","Multiple energy values and multiple angle values"))

# Set constants
r0=2.818e-13
m0c2=511
#Na=6.0221408e+23
Na=6.02e23
hc=12.4

if option == "An energy value and an angle value":
  OPT="op1"
  energy_single = float(st.sidebar.text_input("Energy (kev)", 500.0))
  theta_single = float(st.sidebar.text_input("Angle", 45.0))
elif option == "An energy value and multiple angle values":
  OPT="op2"
  energy_single = float(st.sidebar.text_input("Energy (kev)", 500.0))
  c1, c2, c3 = st.sidebar.columns(3)
  with c1:
    theta1 = float(st.text_input("Angle from", 45.0))
  with c2:
    theta2 = float(st.text_input("Angle to", 90.0))
  with c3:
    theta_skip = float(st.text_input("Angle skip", 5.0))
elif option == "Multiple energy values and an angle value":
  OPT="op3"
  c1, c2, c3 = st.sidebar.columns(3)
  with c1:
    energy1 = float(st.text_input("Energy from (kev)", 100.0))
  with c2:
    energy2 = float(st.text_input("Energy to (kev)", 500.0))
  with c3:
    energy_skip = float(st.text_input("Energy skip", 50.0))
  theta_single = float(st.sidebar.text_input("Angle", 45.0))

elif option == "Multiple energy values and multiple angle values":
  OPT="op4"
  c1, c2, c3 = st.sidebar.columns(3)
  cl1, cl2, cl3 = st.sidebar.columns(3)
  with c1:
    energy1 = float(st.text_input("Energy from (kev)", 100.0))
  with c2:
    energy2 = float(st.text_input("Energy to (kev)", 500.0))
  with c3:
    energy_skip = float(st.text_input("Energy skip", 50.0))
  with cl1:
    theta1 = float(st.text_input("Angle from", 45.0))
  with cl2:
    theta2 = float(st.text_input("Angle to", 90.0))
  with cl3:
    theta_skip = float(st.text_input("Angle skip", 5.0))

def decimal_range(start, stop, increment):
  while start < stop:
    yield start
    start+=increment

theta=[]
energy=[]
if OPT == "op1":
  energy.append(energy_single)
  energy=numpy.array(energy)
  theta.append(theta_single)
  theta=numpy.array(theta)
elif OPT == "op2":
  energy.append(energy_single)
  energy=numpy.array(energy)
  for i in decimal_range(theta1, theta2 + 0.001, theta_skip):
    theta.append(i)
  theta=numpy.array(theta)
elif OPT == "op3":
  theta.append(theta_single)
  theta=numpy.array(theta)
  for j in decimal_range(energy1, energy2 + 0.001, energy_skip):
    energy.append(j)
  energy=numpy.array(energy)
else:
  for j in decimal_range(energy1, energy2 + 0.001, energy_skip):
    energy.append(j)
  energy=numpy.array(energy)
  for i in decimal_range(theta1, theta2 + 0.001, theta_skip):
    theta.append(i)
  theta=numpy.array(theta)

result = {"Energy (kev)": [], "Angle": [], "Mom. Tra.": [], "KN": [], "Th": [], "TMMDSC": []}
for i in range(len(energy)):
  eng=energy[i]
  for j in range(len(theta)):
    tht=theta[j]
    alpha = eng/m0c2
    
    radians=math.radians
    cos=numpy.cos(radians(tht))
    sin=numpy.sin(radians(tht/2))
    
    lmb = hc/eng
    momTra=sin/lmb
    
    sigmaKN = 0.5*r0*r0*(1/((1+alpha*(1-cos))*(1+alpha*(1-cos))))*(1+cos*cos+(alpha*alpha*(1-cos)*(1-cos))/(1+alpha*(1-cos)))
    sigmaTh = 0.5*r0*r0*(1+numpy.cos(radians(tht))**2)
    
    elF=[]
    elS=[]
    for i in range(len(elName)):
      for j in range(len(Fs[0,:])):
        if Fs[0,j] == elName[i]:
          elNameF=numpy.array(numpy.delete(Fs[:,j],0), dtype=float)
          elNameS=numpy.array(numpy.delete(Ss[:,j],0), dtype=float)
          for k in range(0,len(X)):
            if X[k] == momTra:
              elF.append(elNameF[k])
              elS.append(elNameS[k])
            if X[k] < momTra < X[k+1]:
              elF.append((elNameF[k]*(numpy.log(X[k+1])-numpy.log(momTra))+elNameF[k+1]*(numpy.log(momTra)-numpy.log(X[k])))/(numpy.log(X[k+1])-numpy.log(X[k])))
              elS.append((elNameS[k]*(numpy.log(X[k+1])-numpy.log(momTra))+elNameS[k+1]*(numpy.log(momTra)-numpy.log(X[k])))/(numpy.log(X[k+1])-numpy.log(X[k])))
    
    wAF=0
    wAS=0
    if chemical == "A chemical formula":
      for i in range(len(elName)):
        wAF+=(elAtomNumber[i]*elRate[i]/numpy.sum(numpy.multiply(elAtomNumber,elRate)))/elAtomNumber[i]*elF[i]*elF[i]
        wAS+=(elAtomNumber[i]*elRate[i]/numpy.sum(numpy.multiply(elAtomNumber,elRate)))/elAtomNumber[i]*elS[i]
    if chemical == "A chemical composition":
      for i in range(len(elName)):
        wAF+=composition[i]/elAtomNumber[i]*elF[i]*elF[i]
        wAS+=composition[i]/elAtomNumber[i]*elS[i]
    
    TMMDSC=Na*sigmaTh*wAF+Na*sigmaKN*wAS

    result["Energy (kev)"].append(str(eng))
    result["Angle"].append(str(tht))
    result["Mom. Tra."].append(str(momTra))
    result["KN"].append(str(sigmaKN))
    result["Th"].append(str(sigmaTh))
    result["TMMDSC"].append(str(TMMDSC))

if OPT == "op1":
  if chemical == "A chemical formula":
    form = {"Formula": formula}
  elif chemical == "A chemical composition":
    form = {"Element": [], "Composition (%)": []}
    for j in range(len(elName)):
      form["Element"].append(elName[j])
      form["Composition (%)"].append(composition[j]*100)
  data = {"Energy (kev)": eng, "Angle": tht}
elif OPT == "op2":
  if chemical == "A chemical formula":
    form = {"Formula": formula}
  elif chemical == "A chemical composition":
    form = {"Element": [], "Composition (%)": []}
    for j in range(len(elName)):
      form["Element"].append(elName[j])
      form["Composition (%)"].append(composition[j]*100)
  data = {"Energy (kev)": eng, "Angle from": theta1, "Angle to": theta2, "Angle skip": theta_skip}
elif OPT == "op3":
  if chemical == "A chemical formula":
    form = {"Formula": formula}
  elif chemical == "A chemical composition":
    form = {"Element": [], "Composition (%)": []}
    for j in range(len(elName)):
      form["Element"].append(elName[j])
      form["Composition (%)"].append(composition[j]*100)
  data = {"Energy from (kev)": energy1, "Energy to (kev)": energy2, "Energy skip": energy_skip, "Angle": tht}
elif OPT == "op4":
  if chemical == "A chemical formula":
    form = {"Formula": formula}
  elif chemical == "A chemical composition":
    form = {"Element": [], "Composition (%)": []}
    for j in range(len(elName)):
      form["Element"].append(elName[j])
      form["Composition (%)"].append(composition[j]*100)
  data = {"Energy from (kev)": energy1, "Energy to (kev)": energy2, "Energy skip": energy_skip, "Angle from": theta1, "Angle to": theta2, "Angle skip": theta_skip}
  


if chemical == "A chemical formula":
  fTMMDSC = pd.DataFrame(form, index=[1])
elif chemical == "A chemical composition":
  fTMMDSC = pd.DataFrame(form)
  fTMMDSC.index = fTMMDSC.index + 1
dTMMDSC = pd.DataFrame(data, index=[1])

if cite == True:
  if chemical == "A chemical formula":
    st.write("""
    ### Chemical Formula 
    """)
    st.write(fTMMDSC)
  elif chemical == "A chemical composition":
    st.write("""
    ### Chemical Composition 
    """)
    st.write(fTMMDSC)
  st.write("""
  ### Variables 
  """)
  st.write(dTMMDSC)

result_print = pd.DataFrame(result)
result_print.index = result_print.index + 1


if cite == True:
  
  st.sidebar.header("Chemical formula or composition")
  st.write("""
  ### Results
  """)
  st.write(result_print)
  
  def convert_df(df):
     return df.to_csv(index=False, sep=',').encode('utf-8')
  csv = convert_df(result_print)
  st.download_button(
     "Download results",
     csv,
     "table.csv",
     "text/csv",
     key='download-csv'
   )
