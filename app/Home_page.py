# app/app.py

# ------------------HOME PAGE------------------

import streamlit as st

st.set_page_config(page_title="Ionio Bioinformatics", layout="wide")
st.title("🧬 Interactive Machine Learning Platform")

col1, col2, col3 = st.columns([1, 3, 1])
with col2:
    # Replace it Gianni
    st.image("app/developers.jpg", caption="Developers: Matin & Giannis", use_container_width=True)