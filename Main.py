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

# tools tab
with tab2:
    st.header("Expertise", anchor="expertise", divider="grey")
    st.markdown("Here are our areas of expertise:")

    st.header("Overview of Tools", anchor="tools", divider="grey")

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
                image=load_logo("images/Antimicrobial.webp"),
                url="https://antimicrobial-kg.serve.scilifelab.se/",
            )

    st.markdown("### Machine learning models")

    with st.expander(label="ML models"):
        col1, col2, col3 = st.columns(3)

        with col1:
            hasClicked = card(
                title="AntiMicrobial Model",
                text="Model for predicting antimicorbial activity of small molecules",
                image=load_logo("images/Antimicrobial.webp"),
                url="https://antimicrobial-kg.serve.scilifelab.se/Model_Prediction",
            )

        with col2:
            hasClicked = card(
                title="E3 Ligage Model",
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
                title="Patenet Enrichment Tool PEMT",
                text="Python package for extracting patent for genes and compounds",
                image=load_logo("images/pemt.jpg"),
                url="https://github.com/Fraunhofer-ITMP/PEMT",
            )

        with col2:
            card(
                title="rSASC",
                text="Synthetic cohorts generator for longitudinal observational patient cohorts",
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

# publications tab
with tab3:
    # create list with publication data in Chicago style
    st.header("2024 Publications", anchor="2024", divider="grey")
    st.markdown(
        """
        1. Gadiya, Yojana, Simran Shetty, Martin Hofmann-Apitius, Philip Gribbon, and Andrea Zaliani. "Exploring SureChEMBL from a drug discovery perspective." *Scientific Data* 11, no. 1 (2024): 507. https://doi.org/10.1038/s41597-024-03371-4
        
        2. Karki, Reagon, Yojana Gadiya, Simran Shetty, Philip Gribbon, and Andrea Zaliani. "Pharmacophore-based ML model to filter candidate E3 ligands and predict E3 Ligase binding probabilities." *Informatics in Medicine Unlocked* 44 (2024): 101424. https://doi.org/10.1016/j.imu.2023.101424
        """
    )

    st.header("2023 Publications", anchor="2023", divider="grey")
    st.markdown(
        """
        1. Gadiya, Yojana, Andrea Zaliani, Philip Gribbon, and Martin Hofmann-Apitius. "PEMT: a patent enrichment tool for drug discovery." *Bioinformatics* 39, no. 1 (2023): btac716. https://doi.org/10.1093/bioinformatics/btac716
        
        2. Karki, Reagon, Yojana Gadiya, Andrea Zaliani, and Philip Gribbon. "Mpox Knowledge Graph: a comprehensive representation embedding chemical entities and associated biology of Mpox." *Bioinformatics Advances* 3, no. 1 (2023): vbad045. https://doi.org/10.1093/bioadv/vbad045

        3. Karki, Reagon, Yojana Gadiya, Philip Gribbon, and Andrea Zaliani. "Pharmacophore-Based Machine Learning Model To Predict Ligand Selectivity for E3 Ligase Binders." *ACS omega* 8, no. 33 (2023): 30177-30185. https://doi.org/10.1021/acsomega.3c02803

        4. Gadiya, Yojana, Philip Gribbon, Martin Hofmann-Apitius, and Andrea Zaliani. "Pharmaceutical patent landscaping: A novel approach to understand patents from the drug discovery perspective." *Artificial Intelligence in the Life Sciences* 3 (2023): 100069. https://doi.org/10.1016/j.ailsci.2023.100069
        """
    )

    st.header("2022 Publications", anchor="2022", divider="grey")
    st.markdown(
        """
        1. Berg, Hannes, Maria A. Wirtz Martin, Nadide Altincekic, Islam Alshamleh, Jasleen Kaur Bains, Julius Blechar, Betül Ceylan et al. "Comprehensive fragment screening of the SARS‐CoV‐2 proteome explores novel chemical space for drug development." *Angewandte Chemie International Edition* 61, no. 46 (2022): e202205858. https://doi.org/10.1002/anie.202205858
        """
    )

    st.header("2021 Publications", anchor="2021", divider="grey")
    st.markdown(
        """
        1. Schultz, Bruce, Andrea Zaliani, Christian Ebeling, Jeanette Reinshagen, Denisa Bojkova, Vanessa Lage-Rupprecht, Reagon Karki *et al.* "A method for the rational selection of drug repurposing candidates from multimodal knowledge harmonization." *Scientific reports* 11, no. 1 (2021): 11049. https://doi.org/10.1038/s41598-021-90296-2

        2. Khorchani, Takoua, Yojana Gadiya, Gesa Witt, Delia Lanzillotta, Carsten Claussen, and Andrea Zaliani. "SASC: A simple approach to synthetic cohorts for generating longitudinal observational patient cohorts from COVID-19 clinical data." *Patterns* 3, no. 4 (2022). https://doi.org/10.1016/j.patter.2022.100453
        """
    )

# team tab
with tab4:
    st.markdown("### Meet the team!")
    st.markdown(
        "For any questions, feedback or collaborations, please contact our team members."
    )

    df = pd.read_csv("data/members.csv")

    for sub_df in np.array_split(df, 3):
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
