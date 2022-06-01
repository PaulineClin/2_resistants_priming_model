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

    st.markdown("## Prevalence of the disease at the equilibrium without priming. ")
    st.write(r""" Variation of the prevalence  $\mathcal{P}$ as a function of the proportion of resistant 2 in the mixture ($p$).""")
    # Test sur le cas sans priming : 
        
    # Paramètre : 
    R = st.slider('Transmission rate (R):', min_value=1, max_value=100, value = 7)
    c1 = st.slider('Virulence cost 1 (c1):', min_value=0.0, max_value=1.0, value = 0.1 )
    c2 = st.slider('Virulence cost 2 (c2):', min_value=0.0, max_value=1.0, value = 0.4 )
    p = np.linspace(0,1,300)

    # Equilibre singly 1 : 
    y1 = -(R*c1*p - R*c1 - R*p + R - 1)/(R*(c1 - 1))
    y2 = 0
    z1 = 0
    z2 = 0

    # Matrice Jacobienne : 
    def jacobian_multiresistance(y1, y2, z1, z2, R, c1, c2, p):
        return np.array([[R * (-c1 + 1) * (1 - p - y1 - z1) - R * (-c1 + 1) * y1 - 1, 0, -R*(1 - c1)*y1, 0 ],
                        [0, R*(1 - c2)*(p - y2 - z2) - R*(1 - c2)*y2 - 1, 0, -R*(1 - c2)*y2],
                        [-R*(1 - c1)*(1 - c2)*(z1 + z2), 0, R*(1 - c1)*(1 - c2)*(1 - p - y1 - z1) - R*(1 - c1)*(1 - c2)*(z1 + z2) - 1, R*(1 - c1)*(1 - c2)*(1 - p - y1 - z1)], 
                        [0, -R*(1 - c1)*(1 - c2)*(z1 + z2), R*(1 - c1)*(1 - c2)*(p - y2 - z2), R*(1 - c1)*(1 - c2)*(p - y2 - z2) - R*(1 - c1)*(1 - c2)*(z1 + z2) - 1]])

    table = []

    for i in range(len(p)):
        ev = np.linalg.eig(jacobian_multiresistance(y1[i], y2, z1, z2, R, c1, c2, p[i]))[0]
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
    y2 = (R*c2*p - R*p + 1)/(R*(c2 - 1))
    z1 = 0
    z2 = 0
        
    table = []

    for i in range(len(p)):
        ev = np.linalg.eig(jacobian_multiresistance(y1, y2[i], z1, z2, R, c1, c2, p[i]))[0]
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
    y1 = (c2 - p)/c2
    y2 = 0
    z1 = -(R*c1*c2*p - R*c1*p - R*c2*p + R*p - c2)/((c1 - 1)*R*c2)
    z2 = (R*c1*c2*p - R*c1*p - R*c2*p + R*p - c2)/(R*(c1*c2 - c1 - c2 + 1)) 
       
    table = []

    for i in range(len(p)):
        ev = np.linalg.eig(jacobian_multiresistance(y1[i], y2, z1[i], z2[i], R, c1, c2, p[i]))[0]
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
    y2 = (c1 + p - 1)/c1
    z1 = -(R*c1*c2*p - R*c1*c2 - R*c1*p - R*c2*p + R*c1 + R*c2 + R*p - R + c1)/(R*(c1*c2 - c1 - c2 + 1))
    z2 = (R*c1*c2*p - R*c1*c2 - R*c1*p - R*c2*p + R*c1 + R*c2 + R*p - R + c1)/(c1*R*(c2 - 1))
       
    table = []

    for i in range(len(p)):
        ev = np.linalg.eig(jacobian_multiresistance(y1, y2[i], z1[i], z2[i], R, c1, c2, p[i]))[0]
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
    z1 = -(p - 1)*(R*c1*c2 - R*c1 - R*c2 + R - 1)/(R*(c1 - 1)*(c2 - 1)) 
    z2 = p*(R*c1*c2 - R*c1 - R*c2 + R - 1)/(R*(c1*c2 - c1 - c2 + 1))
       
    table = []

    for i in range(len(p)):
        ev = np.linalg.eig(jacobian_multiresistance(y1, y2, z1[i], z2[i], R, c1, c2, p[i]))[0]
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
    y1 = -(R*c1*p - R*c1 - R*p + R - 1)/(R*(c1 - 1))
    y2 = (R*c2*p - R*p + 1)/(R*(c2 - 1))
    z1 = 0
    z2 = 0
       
    table = []

    for i in range(len(p)):
        ev = np.linalg.eig(jacobian_multiresistance(y1[i], y2[i], z1, z2, R, c1, c2, p[i]))[0]
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