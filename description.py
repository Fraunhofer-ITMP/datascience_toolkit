"""Display Page for Description under KGG"""
import datetime
import uuid

import pandas as pd
import streamlit as st

import kgg_utils
import streamlit_antd_components as sac

st.set_page_config(
    layout="wide",
    page_title="KGG-Description",
    page_icon="🕸️",
    initial_sidebar_state="collapsed",
)

st.markdown(
    "<h1 style='text-align: center; color: #149372;'> Knowledge Graph Generator (KGG)</h1> <br>",
    unsafe_allow_html=True,
)

# ---------------- NAVIGATION ---------------- #
nav = sac.buttons(
    [
        "Description",
        "KG Generator",
        "Drug-likeness assessment",
        "KG Visualizer",
    ],
    index=0,  # highlighting "Drug-likeness assessment" (0=Description, 1=KG Generator, 2=Drug-likeness, 3=Visualizer)
    size="25",
    direction="horizontal",
)

if nav == "Description":
    pass
elif nav == "KG Generator":
    st.markdown(
        """<meta http-equiv="refresh" content="0; url='http://localhost:8505'" />""",
        unsafe_allow_html=True,
    )
elif nav == "Drug-likeness assessment":
    st.markdown(
        """<meta http-equiv="refresh" content="0; url='http://localhost:8503'" />""",
        unsafe_allow_html=True,
    )
elif nav == "KG Visualizer":
    st.markdown(
        """<meta http-equiv="refresh" content="0; url='http://localhost:8504'" />""",
        unsafe_allow_html=True,
    )

# ---------------- PAGE CONTENT ---------------- #
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


# footer with text and green background
current_year = datetime.datetime.now().year
st.markdown(
    f"<footer style='background-color: #149372; padding: 10px; border-radius: 10px;'>"
    f"<p style='color: white; text-align: center;'>Fraunhofer ITMP © {current_year}</p>"
    "<p style='color: white; text-align: center;'>This work has been conducted across several key projects in which ITMP has been actively involved.</p>"
    "</footer>",
    unsafe_allow_html=True,
)
