"""App for Knowledge Graph Generator (KGG)"""

import base64
import datetime
import os
import uuid
import zipfile
import pickle
# from zipfile import ZipFile
from io import BytesIO

import pandas as pd
import streamlit as st
import streamlit_antd_components as sac
from pybel.struct.summary import supersummary as ss
from tabulate import tabulate

import kgg_utils

state = st.session_state

st.set_page_config(
    layout="wide",
    page_title="Knowledge Graph Generator (KGG)",
    page_icon="🕸️",
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
    "<h1 style='text-align: center; color: #149372;'> Knowledge Graph Generator (KGG)</h1> <br>",
    unsafe_allow_html=True,
)

if "active_kgg_tab" not in state:
    state["active_kgg_tab"] = None


kgg_tabs = sac.buttons(
    [
        "Description",
        "KG Generator",
        "Drug-likeness assessment",
        "KG Visualizer",
    ],
    size="25",
    direction="horizontal",
)
# tab1, tab2, tab3, tab4 = st.tabs(
#     ["Description", "KG Generator", "Drug-likeness assessment", "KG Visualizer"]
# )

if kgg_tabs != state["active_kgg_tab"]:
    # Reset the session state when the tab changes
    state.clear()
    state["active_kgg_tab"] = kgg_tabs


if kgg_tabs == "Description":
    st.markdown("### What is KGG?")
    st.markdown(
        """**K**nowledge **G**raph **G**enerator (KGG) is a fully automated workflow for creating disease-specific KGs. It is developed for a broad spectrum of researchers and scientists, especially for those who are into pre-clinical drug discovery, understanding disease mechanisms/comorbidity and drug-repurposing. The KGG embeds underlying schema of curated public databases to retrieve relevant knowledge which is regarded as the gold standard for high quality data. The KGG is leveraged on our previous contributions to the BY-COVID project where we developed workflows for identification of bio-active analogs for fragments identified in COVID-NMR studies (Berg, H et al.) and representation of Mpox biology (Karki, R et al.). The programmatic scripts and methods for KGG are written in python (version 3.10) and are available on [GitHub](https://github.com/Fraunhofer-ITMP/kgg).
        """
    )

    st.markdown(
        "Additional information about KGG can be found [here](https://www.infectious-diseases-toolkit.org/showcase/knowledge-graph-generator)"
    )

    kg_worklow = "./images/pages/KGG_workflow_updated.png"
    st.image(
        kg_worklow,
        caption="A schematic representation of the 3 phases of KGG workflow.",
    )

    st.markdown("### Examples and Usecases")

    option = st.selectbox(
        label="Please choose a disease by clicking on the grey area below",
        options=("Alzheimer", "Parkinson", "Depression", "COVID-19"),
        index=0,
    )

    kgg_utils.disease_figures(option)

if kgg_tabs == "KG Generator":
    st.markdown(
        """
        <meta http-equiv="refresh" content="0; url='https://kggapp.serve.scilifelab.se/'" />
        """,
        unsafe_allow_html=True,
    )

if kgg_tabs == "Drug-likeness assessment":
    st.markdown(
        """
        <meta http-equiv="refresh" content="0; url='https://druglikeness.serve.scilifelab.se/'" />
        """,
        unsafe_allow_html=True,
    )

if kgg_tabs == "KG Visualizer":
    st.markdown(
        """
        <meta http-equiv="refresh" content="0; url='https://kgviz.serve.scilifelab.se/'" />
        """,
        unsafe_allow_html=True,
    )


# footer with text and green background
current_year = datetime.datetime.now().year
st.markdown(
    f"<footer style='background-color: #149372; padding: 10px; border-radius: 10px;'>"
    f"<p style='color: white; text-align: center;'>Fraunhofer ITMP © {current_year}</p>"
    "<p style='color: white; text-align: center;'>This work has been conducted across several key projects in which ITMP has been actively involved.</p>"
    "</footer>",
    unsafe_allow_html=True,
)
