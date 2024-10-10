"""Landing page for all Fraunhofer apps."""

import os
import pandas as pd
import numpy as np
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
    "<h1 style='text-align: center; color: #149372;'>Fraunhofer ITMP Data Science Toolkit</h1>",
    unsafe_allow_html=True,
)

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))


def load_logo(filepath):
    with open(filepath, "rb") as f:
        data = f.read()
        encoded = base64.b64encode(data)
    data = "data:image/png;base64," + encoded.decode("utf-8")
    return data


tab1, tab2, tab3, tab4 = st.tabs(["About", "Tools", "Publications", "Team"])

with tab1:
    st.header("Fraunhofer ITMP ScreeningPort", anchor="about", divider="grey")
    st.markdown(
        """Our expertise lies in high-throughput drug research using high-quality substance and repurposing libraries (*in silico* and *in vitro* screening), which enables the identification of pharmacologically active substances. When investigating the mechanisms of action, an extensive portfolio of phenotypic and biochemical assays, as well as in vitro models based on induced pluripotent stem cells, are used. We also develop workflows to ensure the analysis of drug discovery data and the highest standards in FAIR data management, as well as algorithms and AI tools for the statistical analysis of patient cohorts from different medical indications. Our offering thus covers the broad field of medical data science."""
    )

    st.header("Our Data Science Projects", anchor="ds-projects", divider="grey")
    st.markdown(
        """This website hosts the collection of data science tools developed by Fraunhofer ITMP SreeningPort across several key projects. These projects include \n:"""
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
        """* **COMBINE**:  [COMBINE](https://amr-accelerator.eu/project/combine/) is part of the IMI AMR Accelerator and focuses on antimicrobial resistance (AMR). It aims to enhance drug discovery and development by integrating and optimizing preclinical models, biomarkers, and clinical trial designs. The project combines expertise from various research institutions and pharmaceutical companies to improve the development of new treatments for infections caused by drug-resistant bacteria. COMBINE supports the fight against AMR by streamlining processes to accelerate the creation of effective therapies.
        
        Project overview:
        - Start year: 2019
        - End year: 2025
        - Grant agreement ID: 853967
        - Funding: 25M€
        - Funding source: Innovative Medicines Initiative 2
        """
    )

    st.write(
        """* **PROXIDRUGS**:  [PROXIDRUGS](https://www.proxidrugs.de/) is a research initiative under Germany’s Clusters4Future program, focused on developing proximity-induced drugs for treating diseases with unmet medical needs. These drugs trigger the targeted degradation of disease-related proteins, opening new therapeutic avenues, especially for previously undruggable proteins. The project involves collaboration between academic institutions and pharmaceutical companies like AbbVie, Merck, and GSK. It aims to enhance drug development, support innovation, and promote regional collaboration in the pharmaceutical industry.
        
        Project overview:
        - Start year: 2021
        - End year: 2024
        - Grant agreement ID: 03ZU1109KB
        - Funding: 14M€
        - Funding source: Federal Ministry of Education and Research (BMBF)
        """
    )

    st.write(
        """* **IDERHA**:  [IDERHA](https://www.iderha.org/) (Improving Data-Driven Research in Health Applications) is a collaborative initiative aimed at enhancing clinical decision-making and patient access to healthcare innovations by optimizing the use of health data. Supported by the Innovative Health Initiative, the project focuses on developing a federated infrastructure for secure, cross-border health data sharing. It works across sectors like life sciences, healthcare, and policy-making to improve healthcare outcomes.
        
        Project overview:
        - Start year: 2023
        - End year: 2028
        - Grant agreement ID: 101112135
        - Funding: 19M€
        - Funding source: Innovative Health Initiative
        """
    )

    st.write(
        """* **SYNTHIA**:  [SYNTHIA](https://www.ihi-synthia.eu/) focuses on using synthetic data to advance personalized medicine. It aims to develop privacy-preserving, high-quality synthetic datasets for healthcare research, addressing patient privacy concerns while improving treatment options. The project involves collaboration across sectors like Medtech, Pharma, and Academia, providing a platform for generating and validating synthetic data to accelerate discoveries in medicine and enhance patient care. SYNTHIA's federated platform ensures data privacy and quality for research and innovation.

        Project overview:
        - Start year: 2024
        - End year: 2029
        - Grant agreement ID: 101172872
        - Funding: 0.26M€
        - Funding source: Innovative Health Initiative
        """
    )

    st.write(
        """* **FAIRplus**:  [FAIRplus](https://fairplus-project.eu/) focuses on improving the FAIR (Findable, Accessible, Interoperable, Reusable) status of life science data. It developed a reusable FAIRification framework, including tools and guidelines to help organizations manage their data more effectively. By promoting better data sharing and management, the project aims to boost innovation in health research.

        Project overview:
        - Start year: 2019
        - End year: 2022
        - Grant agreement ID: 802750
        - Funding: 7.8M€
        - Funding source: Innovative Medicines Initiative 2
        """
    )

# tools tab
with tab2:
    st.markdown("### Knowledge Graph")

    with st.expander(label="Graph based tools"):
        col1, col2, col3 = st.columns(3)

        with col1:
            hasClicked = card(
                title="R4A Expertise Dashboard",
                text="Dashboard for the the REMEDi4ALL expertise knowledge graph",
                image=load_logo("images/r4a_logo.png"),
                url="https://r4a-expertise-dashboard.serve.scilifelab.se/",
                styles={
                    "card": {
                        "border-radius": "10px",
                        "box-shadow": "0 0 10px rgba(0,0,0,0.5)",
                    }
                },
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
                image=load_logo("images/antimicrobial.webp"),
                url="https://antimicrobial-kg.serve.scilifelab.se/",
            )

    st.markdown("### Machine learning models")

    with st.expander(label="ML models"):
        col1, col2, col3 = st.columns(3)

        with col1:
            hasClicked = card(
                title="AntiMicrobial Model",
                text="Model for predicting antimicorbial activity of small molecules",
                image=load_logo("images/antimicrobial.webp"),
                url="https://antimicrobial-kg.serve.scilifelab.se/Model_Prediction",
            )

        with col2:
            hasClicked = card(
                title="E3 Ligase Model",
                text="Model for predicting speicificity of E3 ligase binders",
                image=load_logo("images/protac.webp"),
                url="https://github.com/Fraunhofer-ITMP/E3_binder_Model",
            )

    st.markdown("### Screening data preprocessing")

    with st.expander("Preprocessing tools"):
        col1, col2, col3 = st.columns(3)

        with col1:
            card(
                title="Dose-Response Analyses",
                text="Fit dose-response curves from normalised biochemical assay data",
                image=load_logo("images/drc.png"),
                url="https://hub.knime.com/fraunhoferitmp/spaces/Public/Dose_Response_Biochemical/DRCfit_biochemical_ECBD~6NLZB5Jkgn6j5a6Y/current-state",
            )

        with col2:
            hasClicked = card(
                title="Chemical annotator",
                text="Annotation of chemicals with biochemical and pharmaceutical information retrieved from ChEMBL",
                image=load_logo("images/workflow.png"),
                url="https://hub.knime.com/fraunhoferitmp/spaces/Public/Remedi4All/R4A_AnnotationTool_v1~RPNHTYMP7vUEoUwD/current-state",
            )

        with col3:
            hasClicked = card(
                title="Annotator Dashboard",
                text="Dashboard of biochemical and pharmaceutical information retrieved from ChEMBL",
                image=load_logo("images/workflow.png"),
                url="https://hub.knime.com/fraunhoferitmp/spaces/Public/Remedi4All/ChemblAssays~yQeFQzSCx6DSlZPd/current-state",
            )

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

    st.markdown("### Programming tools and workflows")

    with st.expander("Program tools"):
        col1, col2, col3 = st.columns(3)

        with col1:
            card(
                title="Patent Enrichment Tool PEMT",
                text="Python package for extracting patent for genes and compounds",
                image=load_logo("images/pemt.jpg"),
                url="https://github.com/Fraunhofer-ITMP/PEMT",
            )

        with col2:
            card(
                title="rSASC",
                text="R tool for Synthetic cohorts generator for longitudinal observational patient cohorts",
                image=load_logo("images/sasc.webp"),
                url="https://github.com/Fraunhofer-ITMP/SASC",
            )

        with col3:
            card(
                title="ALISTER",
                text="Application for lipid stability evaluation",
                image=load_logo("images/alister.png"),
                url="https://itmp.shinyapps.io/alister/",
            )

        col1, col2, col3 = st.columns(3)

        with col1:
            card(
                title="PySASC",
                text="Python tool for Synthetic cohorts generator for longitudinal observational patient cohorts",
                image=load_logo("images/sasc.webp"),
                url="https://github.com/Fraunhofer-ITMP/PySASC",
            )

        # with col3:
        #     card(
        #         title="ALISTER",
        #         text="Application for lipid stability evaluation",
        #         image=load_logo("images/alister.png"),
        #         url="https://itmp.shinyapps.io/alister/",
        #     )

    st.markdown("### Findable, Accessible, Reusable, Interoperable (FAIR) data tools")

    with st.expander(label="Endorsed FAIR tools"):
        col1, col2, col3 = st.columns(3)

        with col1:
            card(
                title="FAIR Cookbook",
                text="Resource for the Life Sciences with recipes that help you to make and keep data FAIR",
                image=load_logo("images/faircookbook.png"),
                url="https://faircookbook.elixir-europe.org/",
                styles={
                    "card": {
                        "border-radius": "10px",
                        "box-shadow": "0 0 10px rgba(0,0,0,0.5)",
                    }
                },
            )

        with col2:
            card(
                title="FAIRplus DSM",
                text="Tool to assess the FAIRness level of datasests",
                image=load_logo("images/fairplus_dsm.png"),
                url="https://fairdsm.biospeak.solutions/",
            )

        with col3:
            card(
                title="FAIRSharing.org",
                text="Resource on data and metadata standards, inter-related to databases and data policies",
                image=load_logo("images/fairsharing.jpg"),
                url="https://fairsharing.org/",
            )


# publications tab
with tab3:
    # create list with publication data in Chicago style
    df = pd.read_csv("data/publications.tsv", sep="\t", dtype=str)
    df = df.sort_values("year", ascending=False)

    for year in df["year"].unique():
        tmp = df[df["year"] == year]
        tmp.reset_index(drop=True, inplace=True)
        st.header(f"{year} Publications", anchor=year, divider="grey")

        for i, row in tmp.iterrows():
            st.markdown(
                f"{i+1}. {row['chicago_authors']} {row['chicago_title']} *{row['chicago_journal']}* {row['chicago_meta']} doi: {row['doi']}"
            )

# team tab
with tab4:
    st.markdown("### Meet the team!")
    st.markdown(
        "For any questions, feedback or collaborations, please contact our team members."
    )

    member_df = pd.read_csv("data/members.csv")

    member_chunk_df = [member_df[i : i + 3] for i in range(0, len(member_df), 3)]

    for sub_df in member_chunk_df:
        sub_df.reset_index(drop=True, inplace=True)
        columns = st.columns(3)
        for i, row in sub_df.iterrows():
            with columns[i]:
                st.write(f"### {row['member_name']}")
                st.markdown(
                    f"**Email**: {row['member_email']} <br> **ORCID**: https://orcid.org/{row['member_orcid']}",
                    unsafe_allow_html=True,
                )

# footer with text and green background
st.markdown(
    "<footer style='background-color: #149372; padding: 10px; border-radius: 10px;'>"
    "<p style='color: white; text-align: center;'>Fraunhofer ITMP © 2024</p>"
    "<p style='color: white; text-align: center;'>This work has been conducted across several key projects in which ITMP has been actively involved.</p>"
    "</footer>",
    unsafe_allow_html=True,
)
