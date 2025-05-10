# app/pages/2_🔬_Run_Pipeline.py
import streamlit as st
from pipeline import run_pipeline

st.title("🔬 Run Analysis Pipeline")

if st.button("Run Pipeline"):
    st.info("Running pipeline, please wait...")
    run_pipeline()
