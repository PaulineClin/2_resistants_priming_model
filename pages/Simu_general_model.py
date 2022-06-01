# -*- coding: utf-8 -*-
"""

Prevalence of the disease according to the number of host genotypes in the mixture. 

@author: clin
"""

# Importation des différents packages : 
import streamlit as st
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd

plt.rc('axes', labelsize=16) 


def app():
    
    st.markdown("## Prevalence of the disease at the equilibrium for all priming values. ")
    st.write(r""" Variation of the prevalence  $\mathcal{P}$ as a function of the proportion of resistant 2 in the mixture ($p$).""")
    
    # Test sur le cas sans priming : 
        
    # Paramètre : 
    R = st.slider('Transmission rate (R):', min_value=1, max_value=100, value = 7)
    c1 = st.slider('Virulence cost 1 (c1):', min_value=0.0, max_value=1.0, value = 0.1 )
    c2 = st.slider('Virulence cost 2 (c2):', min_value=0.0, max_value=1.0, value = 0.4 )
    rho = st.slider('Priming efficiency (rho):', min_value=0.0001, max_value=0.9999, value = 0.5)
    nu = st.slider('Lose of priming (nu):', min_value=1.0000, max_value=10.0000 , value = 1.0001)
    R1 = R*(1-c1)
    R2 = R*(1-c2)
    R3 = R*(1-c1)*(1-c2)
    p = np.linspace(0,1,500)
    
    # Equilibre singly 1 : 
    y1 = -(R1*p - R1 + 1)/R1
    y2 = 0
    z1 = 0
    z2 = 0
    m1 = 0
    m2 = (R1*p - R1 + 1)*p/(R1*p - R1 - nu + 1)
    
    # Matrice Jacobienne : 
    def jacobian_multiresistance(y1, y2, z1, z2, m1, m2, R, c1, c2, p, rho, nu):
        return np.array([[R1*(1 - p - m1 - y1 - z1) - R1*y1 + (1 - rho)*R1*m1 - 1, 0, -R1*y1, 0, -R1*y1 + (1 - rho)*R1*y1, 0], 
                        [0, R2*(p - m2 - y2 - z2) - R2*y2 + (1 - rho)*R2*m2 - 1, 0, -R2*y2, 0, -R2*y2 + (1 - rho)*R2*y2], 
                        [-R3*(z1 + z2), 0, R3*(1 - p - m1 - y1 - z1) - R3*(z1 + z2) + (1 - rho)*R3*m1 - 1, R3*(1 - p - m1 - y1 - z1) + (1 - rho)*R3*m1, -R3*(z1 + z2) + (1 - rho)*R3*(z1 + z2), 0],
                        [0, -R3*(z1 + z2), R3*(p - m2 - y2 - z2) + (1 - rho)*R3*m2, R3*(p - m2 - y2 - z2) - R3*(z1 + z2) + (1 - rho)*R3*m2 - 1, 0, -R3*(z1 + z2) + (1 - rho)*R3*(z1 + z2)],
                        [-R2*y2 - (1 - rho)*R1*m1, R2*(1 - p - m1 - y1 - z1), -R2*y2 - (1 - rho)*R3*m1, -(1 - rho)*R3*m1, -R2*y2 - (1 - rho)*R1*y1 - (1 - rho)*R3*(z1 + z2) - nu, 0], 
                        [R1*(p - m2 - y2 - z2), -R1*y1 - (1 - rho)*R2*m2, -(1 - rho)*R3*m2, -R1*y1 - (1 - rho)*R3*m2, 0, -R1*y1 - (1 - rho)*R2*y2 - (1 - rho)*R3*(z1 + z2) - nu]])
                        
    table = []
    
    for i in range(len(p)):
        ev = np.linalg.eig(jacobian_multiresistance(y1[i], y2, z1, z2, m1, m2[i], R, c1, c2, p[i], rho, nu))[0]
        if all( ev < 0):
            nature = 'Stable'
        else:
            nature = 'Unstable'
        dic = {'p': p[i], 'Prevalence': y1[i], 'Stability': nature}
        table.append(dic)
    
    data1 = pd.DataFrame(table)
    
    groups = data1.groupby("Stability")
    colors = {'Stable':'red', 'Unstable':'lightgrey'}
    
    fig1, ax1 = plt.subplots()
    
    for name, group in groups:
        plt.plot(group['p'][data1.Prevalence > 0], group['Prevalence'][data1.Prevalence > 0], marker="", markersize = 1, linestyle="-", label=name, color=colors[name])
    
    # Equilibre singly 2
    y1 = 0
    y2 = (R2*p - 1)/R2
    z1 = 0
    z2 = 0
    m1 = -(R2*p**2 - R2*p - p + 1)/(R2*p + nu - 1)
    m2 = 0 
        
    table = []
    
    for i in range(len(p)):
        ev = np.linalg.eig(jacobian_multiresistance(y1, y2[i], z1, z2, m1[i], m2, R, c1, c2, p[i], rho, nu))[0]
        if all( ev < 0):
            nature = 'Stable'
        else:
            nature = 'Unstable'
        dic = {'p': p[i], 'Prevalence': y2[i], 'Stability': nature}
        table.append(dic)
    
    data2 = pd.DataFrame(table)
    
    groups = data2.groupby("Stability")
    colors = {'Stable':'red', 'Unstable':'lightgrey'}
    
    for name, group in groups:
        plt.plot(group['p'][data2.Prevalence > 0], group['Prevalence'][data2.Prevalence > 0], marker="", markersize = 1, linestyle="-", label=name, color=colors[name])
    
    # Equilibre singly 1 and doubly :
    y1 = -(R1**2*p**2*rho - R1**2*p**2 - 2*R1**2*p*rho + R1*R3*p*rho + 2*R1**2*p + R1**2*rho - R1*R3*p - R1*R3*rho + R1*nu*p + R1*p*rho - R1**2 + R1*R3 - R1*nu - R1*p - R1*rho + R3*nu + R3*rho + R1 - R3)/(R1**2*p*rho - R1*R3*p*rho - R1**2*p - R1**2*rho + R1*R3*p + R1*R3*rho + R1**2 - R1*R3 + R1*nu - R3*nu - R1 + R3)
    y2 = 0
    z1 = (R1**2*R3*p**2*rho - R1**2*R3*p**2 - R1**2*R3*p*rho + R1**2*R3*p + R1*R3*nu*p + R1*R3*p*rho + R1**2*p - 2*R1*R3*p - R1**2 + R1*R3 - R1*nu + R3*nu + R1 - R3)/(R1*(R1**2*p*rho - R1*R3*p*rho - R1**2*p - R1**2*rho + R1*R3*p + R1*R3*rho + R1**2 - R1*R3 + R1*nu - R3*nu - R1 + R3))
    z2 = (R1**2*R3*p**2*rho - R1**2*R3*p**2 - R1**2*R3*p*rho + R1**2*R3*p + R1*R3*nu*p + R1*R3*p*rho + R1**2*p - 2*R1*R3*p - R1**2 + R1*R3 - R1*nu + R3*nu + R1 - R3)/(R3*R1*(R1*p*rho - R1*p - R1*rho + R1 + nu - 1))
    m1 = 0 
    m2 = -(R1*p - R1 + R3)/(R3*(R1*p*rho - R1*p - R1*rho + R1 + nu - 1))
    
    table = []
    
    for i in range(len(p)):
        ev = np.linalg.eig(jacobian_multiresistance(y1[i], y2, z1[i], z2[i], m1, m2[i], R, c1, c2, p[i], rho, nu))[0]
        if all( ev < 0):
            nature = 'Stable'
        else:
            nature = 'Unstable'
        dic = {'p': p[i], 'Prevalence': y1[i]+z1[i]+z2[i], 'Stability': nature}
        table.append(dic)
    
    data3 = pd.DataFrame(table)
    
    if all(data3.Stability == 'Unstable'):
        plt.plot(data3['p'],data3['Prevalence'], 'lightgrey')
    else :      
        plt.plot(data3['p'][data3.Stability == 'Stable'],data3['Prevalence'][data3.Stability == 'Stable'], 'gold')
     
       
    # Equilibre singly 2 and doubly :
    y1 = 0
    y2 = (R2**2*p**2*rho - R2**2*p**2 - R2*R3*p*rho + R2*R3*p - R2*nu*p - R2*p*rho + R2*p + R3*nu + R3*rho - R3)/(R2**2*p*rho - R2*R3*p*rho - R2**2*p + R2*R3*p - R2*nu + R3*nu + R2 - R3) 
    z1 = -(R2**2*R3*p**2*rho - R2**2*R3*p**2 - R2**2*R3*p*rho + R2**2*R3*p - R2*R3*nu*p - R2*R3*p*rho - R2**2*p + R2*R3*nu + 2*R2*R3*p + R2*R3*rho - R2*R3 - R2*nu + R3*nu + R2 - R3)/(R3*R2*(R2*p*rho - R2*p - nu + 1))
    z2 = -(R2**2*R3*p**2*rho - R2**2*R3*p**2 - R2**2*R3*p*rho + R2**2*R3*p - R2*R3*nu*p - R2*R3*p*rho - R2**2*p + R2*R3*nu + 2*R2*R3*p + R2*R3*rho - R2*R3 - R2*nu + R3*nu + R2 - R3)/(R2*(R2**2*p*rho - R2*R3*p*rho - R2**2*p + R2*R3*p - R2*nu + R3*nu + R2 - R3))
    m1 = -(R2*p - R3)/(R3*(R2*p*rho - R2*p - nu + 1)) 
    m2 = 0
    
    table = []
    
    for i in range(len(p)):
        ev = np.linalg.eig(jacobian_multiresistance(y1, y2[i], z1[i], z2[i], m1[i], m2, R, c1, c2, p[i], rho, nu))[0]
        if all( ev < 0):
            nature = 'Stable'
        else:
            nature = 'Unstable'
        dic = {'p': p[i], 'Prevalence': y2[i]+z1[i]+z2[i], 'Stability': nature}
        table.append(dic)
    
    data4 = pd.DataFrame(table)
    
    if all(data4.Stability == 'Unstable'):
        plt.plot(data4['p'],data4['Prevalence'], 'lightgrey')
    else :      
        plt.plot(data4['p'][data4.Stability == 'Stable'],data4['Prevalence'][data4.Stability == 'Stable'], 'gold')
     
        
    # Equilibre doubly :
    y1 = 0
    y2 = 0
    z1 = -(R3 - 1)*(p - 1)/R3
    z2 = p*(R3 - 1)/R3
    m1 = 0
    m2 = 0
       
    table = []
    
    for i in range(len(p)):
        ev = np.linalg.eig(jacobian_multiresistance(y1, y2, z1[i], z2[i], m1, m2, R, c1, c2, p[i], rho, nu))[0]
        if all( ev < 0):
            nature = 'Stable'
        else:
            nature = 'Unstable'
        dic = {'p': p[i], 'Prevalence': z1[i]+z2[i], 'Stability': nature}
        table.append(dic)
    
    data5 = pd.DataFrame(table)
    
    if all(data5.Stability == 'Unstable'):
        plt.plot(data5['p'],data5['Prevalence'], 'lightgrey')
    else :      
        plt.plot(data5['p'][data5.Stability == 'Stable'],data5['Prevalence'][data5.Stability == 'Stable'], 'green')
     
       
    
    # Equilibre singly 1 and 2 :
    y1 = (-((c1 - c2) * p - c1 + 1) * (-1 + c1) * (rho - 1) ** 2 * (p - 1) * R ** 2 + 2 * ((-1 + c1) * (p - 1) * rho / 2 + ((c1 + c2 / 2 - 3 / 2) * nu - 3 / 2* c1 + c2 / 2 + 1) * p - (-1 + c1) * (nu - 3 / 2)) * (rho - 1) * R + (-nu + 1) * rho + nu ** 2 - np.sqrt(((rho - 1) * (-1 + c1) * (p - 1) * R - nu + 1) ** 2 * (((p - 1) * c1 - p * c2 + 1) ** 2 * (rho - 1) ** 2 * R ** 2 - 2 * ((p - 1) * c1 - p * c2 + 1) * (rho + nu - 2) * (rho - 1) * R + rho ** 2 + (6 * nu - 4) * rho + (nu - 2) ** 2)) + nu - 2) / (-1 + c1) / (rho - 1) / (((c1 - c2) * p - c1 + 1) * (rho - 1) * R - 2 * nu + 2) / R / 2
    y2 = ((p * (rho - 1) * (c2 - 1) * R + nu - 1) * np.sqrt(((rho - 1) * (-1 + c1) * (p - 1) * R - nu + 1) ** 2 * (((p - 1) * c1 - p * c2 + 1) ** 2 * (rho - 1) ** 2 * R ** 2 - 2 * ((p - 1) * c1 - p * c2 + 1) * (rho + nu - 2) * (rho - 1) * R + rho ** 2 + (6 * nu - 4) * rho + (nu - 2) ** 2)) + ((rho - 1) * (-1 + c1) * (p - 1) * R - nu + 1) * (p * ((c1 - c2) * p - c1 + 1) * (rho - 1) ** 2 * (c2 - 1) * R ** 2 - (rho - 1) * (((c2 - 1) * rho + (nu + 1) * c1 + (2 * nu - 3) * c2 - 3 * nu + 2) * p - (nu + 1) * (-1 + c1)) * R + (nu - 1) * (nu - rho + 2))) / (rho - 1) / (c2 - 1) / (((c1 - c2) * p - c1 + 1) * (rho - 1) * R - 2 * nu + 2) / R / ((rho - 1) * (-1 + c1) * (p - 1) * R - nu + 1) / 2
    z1 = 0
    z2 = 0
    m1 = (-(rho - 1) ** 2 * (p - 1) ** 2 * R1 ** 2 + (rho - 1) * (p - 1) * (R2 * (rho - 1) * p - 2 * nu - rho + 3) * R1 + R2 * (rho - 1) * (nu + 2 * rho - 1) * p + (-3 * nu + 3) * rho - nu ** 2 + np.sqrt((((R1 - R2) * p - R1 + 1) ** 2 * rho ** 2 + (-2 * (R1 - R2) ** 2 * p ** 2 + 4 * (R1 + nu / 2 - 3/ 3) * (R1 - R2) * p + (-2 * R1 + 6) * nu - 2 * R1 ** 2 + 6 * R1 - 4) * rho + ((R1 - R2) * p - R1 - nu + 2) ** 2) * (R1 * (p - 1) * rho - R1 * p + R1 + nu - 1) ** 2) + 3 * nu - 2) / (rho - 1) / rho / R1 / ((rho - 1) * (p - 1) * R1 - R2 * (rho - 1) * p + 2 * nu - 2) / 2
    m2 = ((-p * (rho - 1) * (c2 - 1) * R - nu + 1) * np.sqrt(((rho - 1) * (-1 + c1) * (p - 1) * R - nu + 1) ** 2 * (((p - 1) * c1 - p * c2 + 1) ** 2 * (rho - 1) ** 2 * R ** 2 - 2 * ((p - 1) * c1 - p * c2 + 1) * (rho + nu - 2) * (rho - 1) * R + rho ** 2 + (6 * nu - 4) * rho + (nu - 2) ** 2)) + (p * ((c1 - c2) * p - c1 + 1) * (rho - 1) ** 2 * (c2 - 1) * R ** 2 + (((2 * c1 - c2 - 1) * rho + (nu - 1) * c1 + (-2 * nu + 3) * c2 + nu - 2) * p - (2 * rho + nu - 1) * (-1 + c1)) * (rho - 1) * R - (nu - 1) * (nu + 3 * rho - 2)) * ((rho - 1) * (-1 + c1) * (p - 1) * R - nu + 1)) / rho / (rho - 1) / (c2 - 1) / (((c1 - c2) * p - c1 + 1) * (rho - 1) * R - 2 * nu + 2) / R / ((rho - 1) * (-1 + c1) * (p - 1) * R - nu + 1) / 2
    
       
    table = []
    
    for i in range(len(p)):
        ev = np.linalg.eig(jacobian_multiresistance(y1[i], y2[i], z1, z2, m1[i], m2[i], R, c1, c2, p[i], rho, nu))[0]
        if all( ev < 0):
            nature = 'Stable'
        else:
            nature = 'Unstable'
        dic = {'p': p[i], 'Prevalence': y1[i]+y2[i], 'Stability': nature}
        table.append(dic)
    
    data6 = pd.DataFrame(table)
    
    groups = data6.groupby("Stability")
    colors = {'Stable':'blue', 'Unstable':'lightgrey'}
    
    if all(data6.Stability == 'Unstable'):
        for name, group in groups:
            plt.plot(group['p'][data6.Prevalence > 0], group['Prevalence'][data6.Prevalence > 0], marker="", markersize = 1, linestyle="-", label='Singly 1', color=colors[name])
    else: 
        plt.plot(data6['p'][data6.Stability == 'Stable'],data6['Prevalence'][data6.Stability == 'Stable'], 'blue')
    
    ax1.set_ylim(0, 1)
    ax1.set_xlim(0,1)
    ax1.set_xlabel('Proportion of resistant 2 (p)')
    ax1.set_ylabel('Disease prevalence (P)')

    # Show the pyplot figure in the app : 
    st.pyplot(fig1)
    
    st.caption(r""" Total equilibrium prevalence of the disease $\mathcal{P}$ as a function of the proportion of resistant 2 in the mixture ($p$).
               Color lines show the stable part of each equilibrium. Red lines correspond to singly 1 and sngly 2 equilibra. Yellow lines correspond to the two singly and doubly equilibria.
               Green line corresponds to the doubly equilibrium. Blue line corresponds to singly 1 and singly 2 equilbirum. The grey lines correspond to the unstable part of each equilibrium.
               For singly and doubly, and doubly equilibria, the unstable part are not shown.""")