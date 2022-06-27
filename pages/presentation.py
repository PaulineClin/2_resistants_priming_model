import streamlit as st

def app():
    st.markdown("<h1 style='text-align: center;'> Host mixtures should not balance disease resistances, and should protect the one with the highest breaking cost </h1>", unsafe_allow_html=True)
  
    st.markdown("<h3 style='text-align: center;'> Pauline Clin, Frédéric Grognard, Didier Andrivon, Ludovic Mailleret, Frédéric Hamelin </h3>", unsafe_allow_html=True)
   
    st.markdown("This application is linked to our study submitted for publication. ")
    st.markdown("This application allows the user to produce the different figures presented in the manuscript for parameter values that can be chosen by the user.")
    st.markdown("Click on the left panel to choose the topic you are interested in:")
    
    st.write(r"""
    * **Primingless model:** 
      * Priming efficiency is null ($\rho=0$).  
    * **General model:**
      * Priming efficiency varies ($0<\rho<1$).
   """)

    
    
