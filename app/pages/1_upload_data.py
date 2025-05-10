# app/pages/1_📁_Upload_Data.py
import streamlit as st
import pandas as pd

st.title("📁 Upload Your Data")

uploaded_file = st.file_uploader("Upload a CSV file", type=["csv"])
if uploaded_file:
    df = pd.read_csv(uploaded_file)
    st.success("Upload successful!")
    st.dataframe(df)
