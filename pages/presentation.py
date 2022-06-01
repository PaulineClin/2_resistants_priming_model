import streamlit as st

def app():
    st.markdown("<h1 style='text-align: center;'> The optimal ratio of disease-resistant plants in host mixtures is imbalanced and counter-intuitively biased </h1>", unsafe_allow_html=True)
  
    st.markdown("<h3 style='text-align: center;'> Clin Pauline, Grognard Frédéric, Andrivon Didier, Mailleret Ludovic, Hamelin Frédéric</h3>", unsafe_allow_html=True)
   
    st.markdown("This application is linked to our study submitted for publication to the journal Biology Letters.")
    st.markdown("This application allows to produce the different figures presented in the companion article for parameter values that can be chosen by the user.")
    st.markdown("Click on the left panel to choose the topic you are interested in:")
    
    st.write(r"""
    * **Prevalence of the disease without priming:** 
      * Priming efficiency is null ($\rho=0$).  
    * **Prevalence of the disease with full efficient priming:**
      * Priming efficiency is full ($\rho=1$).
    * **Prevalence of the disease :**
      * Priming efficiency varies ($0<\rho<1$).
   """)

    
    
