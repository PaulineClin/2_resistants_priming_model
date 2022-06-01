import streamlit as st

# Custom imports 
from multipage import MultiPage
from pages import presentation, Simu_primingless, Simu_fullpriming, Simu_general_model

# Create an instance of the app 
app = MultiPage()

# Title of the main page
st.title("")

# Add all your applications (pages) here
app.add_page("Presentation", presentation.app)
app.add_page("Primingless model", Simu_primingless.app)
app.add_page("Fullpriming model", Simu_fullpriming.app)
app.add_page("General model", Simu_general_model.app)

# The main app
app.run()