"""Landing page for all Fraunhofer apps."""

import os
import streamlit as st
import base64
from streamlit_card import card

st.set_page_config(
    layout="wide",
    page_title="ITMP Data Science Toolkits",
    page_icon="🧊",
    initial_sidebar_state="collapsed",
    menu_items={
        "Get Help": "https://www.extremelycoolapp.com/help",
        "Report a bug": "https://www.extremelycoolapp.com/bug",
        "About": "# This is a header. This is an *extremely* cool app!",
    },
)

st.markdown(
    "<h1 style='text-align: center; color: #149372;'>Fraunhofer ITMP Data Science Toolkits</h1>",
    unsafe_allow_html=True,
)

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))


def load_logo(filepath):
    with open(filepath, "rb") as f:
        data = f.read()
        encoded = base64.b64encode(data)
    data = "data:image/png;base64," + encoded.decode("utf-8")
    return data


tab1, tab2, tab3 = st.tabs(["About", "Tools", "Contact"])

with tab1:
    st.markdown("### About")
    st.markdown(
        """This website is a collection of data science tools developed by Fraunhofer ITMP SreeningPort across several key projects in which ITMP has been actively involved. These projects include \n:"""
    )

    st.write(
        """* **REMEDi4ALL**:  [REMEDi4ALL](https://remedi4all.org/) is an EU-funded initiative whose key mission is to make it easier and more reliable to find new medical uses for drugs we already know are safe and effective. We also aim to show that in many cases such “repurposed” medicines can be taken all the way into clinic faster and at the fraction of the cost of developing a completely new drug from scratch. 
        
        Project overview:
        - Start year: 2022
        - End year: 2027
        - Grant agreement ID: 101057442
        - Funding: 23M€
        - Funding source: HORIZON EUROPE
        """
    )

    st.markdown(
        """
        """
    )

    st.write(
        """* **BY-COVID**:  [BeYond-COVID](https://by-covid.org/) is a European initiative aimed at enhancing COVID-19 data sharing across sectors. It focuses on providing open access to data on COVID-19 variants, treatments, and other aspects of the pandemic. The project seeks to improve data infrastructure, enable effective research collaborations, and make data FAIR (Findable, Accessible, Interoperable, and Reusable). It involves various stakeholders, including public health institutions, researchers, and policymakers, to ensure rapid responses to current and future pandemics.
        
        Project overview:
        - Start year: 2021
        - End year: 2024
        - Grant agreement ID: 101046203
        - Funding: 12M€
        - Funding source: European Union
        """
    )

    st.write(
        """* **COMBINE**:  [COMBINE](https://amr-accelerator.eu/project/combine/) project is part of the IMI AMR Accelerator and focuses on antimicrobial resistance (AMR). It aims to enhance drug discovery and development by integrating and optimizing preclinical models, biomarkers, and clinical trial designs. The project combines expertise from various research institutions and pharmaceutical companies to improve the development of new treatments for infections caused by drug-resistant bacteria. COMBINE supports the fight against AMR by streamlining processes to accelerate the creation of effective therapies.
        
        Project overview:
        - Start year: 2019
        - End year: 2025
        - Grant agreement ID: 853967
        - Funding: 25M€
        - Funding source: Innovative Medicines Initiative 2
        """
    )

with tab2:
    st.markdown("## Tools")
    st.markdown("""The Data Science Toolkits include the following applications:""")

    st.markdown("### Knowledge Graph")

    col1, col2, col3 = st.columns(3)

    with col1:
        hasClicked = card(
            title="R4A Expertise Dashboard",
            text="Dashboard for the the Expertise knowledge graph",
            image=load_logo("images/r4a_logo.png"),
            url="https://r4a-expertise-dashboard.serve.scilifelab.se/",
        )

    with col2:
        hasClicked = card(
            title="KGG",
            text="Automated workflow for disease-specific KGs",
            image=load_logo("images/kg_2.png"),
            url="/KGGapp",
        )

    with col3:
        hasClicked = card(
            title="AntiMicrobial KG",
            text="Data warehouse of experimentally validated antibacterial chemicals",
            image=load_logo("images/Antimicrobial.webp"),
            url="https://antimicrobial-kg.serve.scilifelab.se/",
        )

    st.markdown("### Machine learning models")

    col1, col2, col3 = st.columns(3)

    with col1:
        hasClicked = card(
            title="AntiMicrobial Model",
            text="Model for predicting antimicorbial activity of small molecules",
            image=load_logo("images/Antimicrobial.webp"),
            url="https://antimicrobial-kg.serve.scilifelab.se/Model_Prediction",
        )

    st.markdown("### Screening data preprocessing")

    col1, col2, col3 = st.columns(3)

    with col1:
        hasClicked = card(
            title="Pareto Hit Analysis App",
            text="Some description",
            image=load_logo("images/pareto_front.png"),
            url="/Paretoapp",
        )

    with col2:
        hasClicked = card(
            title="ENSEMBL Gene Annotator App",
            text="Some description",
            image=load_logo("images/gene_icon.png"),
            url="/Ensemblapp",
        )

    with col3:
        hasClicked = card(
            title="Chemical Data Annotator",
            text="Workflow for annotation of compound libraries",
            image=load_logo("images/workflow.png"),
            url="/Ensemblapp",
        )

with tab3:
    st.markdown("### Contact")
    st.markdown("For any questions or feedback, please contact us.")

# footer with text and green background
st.markdown(
    "<footer style='background-color: #006c8b; padding: 10px; border-radius: 10px;'>"
    "<p style='color: white; text-align: center;'>Fraunhofer ITMP © 2024</p>"
    "<p style='color: white; text-align: center;'>This work has been conducted across several key projects in which ITMP has been actively involved.</p>"
    "</footer>",
    unsafe_allow_html=True,
)
