"""Display Page for Proxidrugs"""

import streamlit as st

st.set_page_config(
    layout="wide",
    page_title="Proxidrugs",
    initial_sidebar_state="collapsed",
)

st.markdown(
    """
        <style>
            .block-container {
                padding-top: 1.5rem;
                padding-bottom: 1.5rem;
                padding-left: 5rem;
                padding-right: 5rem;
            }
            .stTabs [data-baseweb="tab-list"] button [data-testid="stMarkdownContainer"] p {font-size:1.3rem;}
        </style>
        """,
    unsafe_allow_html=True,
)  # .block-conatiner controls the padding of the page, .stTabs controls the font size of the text in the tabs

st.markdown(
    "<h1 style='text-align: center; color: #FFFFFF;'> PROXIDRUGS </h1> <br>",
    unsafe_allow_html=True,
)
st.markdown(
    "<h2 style='text-align: center; color: #1B3C67;'> This page is currently being updated. Please visit at a later time point. Thanks! </h2> <br>",
    unsafe_allow_html=True,
)
