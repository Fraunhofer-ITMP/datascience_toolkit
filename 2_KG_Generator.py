"""App for Knowledge Graph Generator (KGG)"""

import base64
import datetime
import os
import uuid
import zipfile

# from zipfile import ZipFile
from io import BytesIO

import pandas as pd
import streamlit as st
from pybel.struct.summary import supersummary as ss
from tabulate import tabulate

import kgg_utils

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

tab1, tab2, tab3 = st.tabs(["Description", "KG Generator", "Drug-likeness assessment"])


with tab1:
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

with tab2:
    st.markdown(
        "👋🏻 Welcome to the KG Generator tool! Interested in building your own disease graph, lets get started..."
    )

    st.markdown("The workflow is divided into 3 main steps:")
    st.markdown(
        "1. **Disease search and Chemical filters** - Search for a disease and select the clinical trial phase of chemicals."
    )
    st.markdown(
        "2. **Generating the graph** - Generate the base knowledge graph with protein filter"
    )
    st.markdown(
        "3. **Finalizing the graph** - Add additional information to the graph and save the files."
    )

    st.header(
        "Disease search and Chemical filters",
        anchor="disease-search",
        divider="grey",
    )

    disease_name = st.text_input(
        "Enter the disease your are interested in generating a graph.",
        value=st.session_state.get("user_disease", "COVID-19"),
    )
    if disease_name == "":
        st.warning("Please enter a disease name.")
        st.stop()
    if "user_disease" not in st.session_state:
        st.session_state["user_disease"] = disease_name

    if st.session_state["user_disease"] != disease_name:
        st.session_state["user_disease"] = disease_name

    disease_df = kgg_utils.searchDisease(st.session_state["user_disease"])
    if disease_df.empty:
        st.write("No results found for the disease. Please try again.")
        st.stop()

    st.expander("Disease search results").dataframe(disease_df, hide_index=True)

    if "disease_df" not in st.session_state:
        st.session_state["disease_df"] = disease_df

    if not st.session_state["disease_df"].equals(disease_df):
        st.session_state["disease_df"] = disease_df

    col1, col2 = st.columns(2)

    with col1:
        st.markdown(
            "**Please select the disease of interest from the table above by entering the identifier value.** \n Each disease has a unique id, which fetches associated proteins using OpenTargets API. Double click twice on the id of interest to select, copy, and paste the id into the text box."
        )
        disease_id = st.text_input(
            "Please enter the identifier for disease of interest.",
            value=disease_df.iloc[0]["id"],
        )

        if "disease_id" not in st.session_state:
            st.session_state["disease_id"] = disease_id
            st.session_state["disease_name"] = disease_df[
                disease_df["id"] == disease_id
            ]["name"].values[0]

        elif disease_id != st.session_state["disease_id"]:
            st.session_state["disease_id"] = disease_id  # Update the session
            # st.session_state["disease_name"] = st.session_state["disease_df"][st.session_state["disease_df"]["id"] == disease_id]["name"].values[0]
            st.session_state["disease_name"] = disease_df[
                disease_df["id"] == disease_id
            ]["name"].values[0]

    with col2:
        st.markdown(
            """**Please enter the clinical trial phase of chemicals which should be used by the workflow.** \n Use a number between 1 (early phase) and 4 (FDA approved). Keep in mind, lower input values increase the number of identified chemicals and running time."""
        )

        ct_phase = st.number_input(
            label="Select your clinical trial phase", min_value=1, max_value=4, value=3
        )
        if "ct_phase" not in st.session_state:
            st.session_state["ct_phase"] = ct_phase
        elif ct_phase != st.session_state["ct_phase"]:
            st.session_state["ct_phase"] = ct_phase  # Update the session

    st.markdown(
        ":red[Selected disease:] "
        + st.session_state["disease_name"]
        + ":red[ with ID:] "
        + st.session_state["disease_id"]
        + ":red[ and clinical trial phase:] "
        + str(st.session_state["ct_phase"])
    )

    viral_prot = kgg_utils.GetViralProteins(st.session_state["user_disease"])
    if "viral_prot" not in st.session_state:
        st.session_state["viral_prot"] = viral_prot
    elif not st.session_state["viral_prot"] == viral_prot:
        # elif viral_prot != st.session_state["viral_prot"]:
        st.session_state["viral_prot"] = viral_prot

    # st.markdown(
    #     ":red[Selected disease:] "
    #     + st.session_state["disease_name"]
    #     + ":red[ with ID:] "
    #     + st.session_state["disease_id"]
    #     + ":red[ and clinical trial phase:] "
    #     + str(st.session_state["ct_phase"])
    # )
    st.write(st.session_state)

    st.header("Generating the graph", anchor="generate-graph", divider="grey")
    dis_name = (
        st.session_state["user_disease"].lower().replace(" ", "_").replace("-", "_")
    )
    kg_name = f"kgg_{dis_name}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
    if "kg_name" not in st.session_state:
        st.session_state["kg_name"] = kg_name
    elif kg_name != st.session_state["kg_name"]:
        st.session_state["kg_name"] = kg_name

    if "button_clicked" not in st.session_state:
        st.session_state.button_clicked = False

    def delete_old_cache_and_files():
        """
        This function deletes all graphs and photos that were just created. This function will be called when the user clicks the "Start over" button or updates the threshold score.
        """
        for key in st.session_state.keys():
            if key.startswith("graph_") or key.startswith("dis2prot_"):
                del st.session_state[key]
        st.cache_resource.clear()

    def callback():
        st.session_state.button_clicked = True

    if (
        st.button("Generate Base Knowledge Graph", on_click=callback)
        or st.session_state.button_clicked
    ):
        # st.write(st.session_state)
        st.write(st.session_state["user_disease"])
        st.write(st.session_state["disease_id"])
        drugs_df, dis2prot_df, dis2snp_df = kgg_utils.createInitialKG(
            disease_id=st.session_state["disease_id"],
            ct_phase=st.session_state["ct_phase"],
        )

        if "drugs_df" not in st.session_state:
            st.session_state["drugs_df"] = drugs_df
        elif not st.session_state["drugs_df"].equals(drugs_df):
            st.session_state["drugs_df"] = drugs_df

        if "dis2prot_df" not in st.session_state:
            st.session_state["dis2prot_df"] = dis2prot_df
        elif not st.session_state["dis2prot_df"].equals(dis2prot_df):
            st.session_state["dis2prot_df"] = dis2prot_df

        if "dis2snp_df" not in st.session_state:
            st.session_state["dis2snp_df"] = dis2snp_df
        elif not st.session_state["dis2snp_df"].equals(dis2snp_df):
            st.session_state["dis2snp_df"] = dis2snp_df

        #        st.session_state["dis2prot_df"] = dis2prot_df
        #        st.session_state["dis2snp_df"] = dis2snp_df
        kgg_utils.GetDiseaseAssociatedProteinsPlot(st.session_state["dis2prot_df"])

        score = st.number_input(
            "Enter threshold score (recommended > 0.3):",
            min_value=0.0,
            max_value=1.0,
            value=0.55,
            step=0.1,
        )
        if "protein_score" not in st.session_state:
            st.session_state["protein_score"] = score
        elif score != st.session_state["protein_score"]:
            st.session_state["protein_score"] = score

        # st.write(st.session_state)

        if st.button("Submit"):
            dis2prot_df = st.session_state["dis2prot_df"]
            disease_df = st.session_state["disease_df"]
            dis2snp_df = st.session_state["dis2snp_df"]

            filtered_df = dis2prot_df[dis2prot_df["Score"] >= score]

            st.warning(
                f"""ℹ️ Filter of protein score (> {score}) has been applied. This reduced the number of proteins from {len(dis2prot_df)} to {len(filtered_df)}."""
            )

            st.session_state.filtered_protein_df = filtered_df

            graph = kgg_utils.finalizeKG(filtered_df, session_inputs=st.session_state)

            st.session_state["graph"] = graph

            st.header("Graph summary", anchor="graph-summary", divider="grey")

            rv_base, rv_stats = kgg_utils.get_graph_summary(st.session_state["graph"])

            with st.container():
                st.markdown(
                    f"""\
                    <h3>Metadata</h3>
                    {tabulate(rv_base, tablefmt="html")}
                    <h3>Statistics</h3>
                    {tabulate(rv_stats, tablefmt="html")}
                    <h3>Nodes</h3>
                    {ss.functions_str(graph, examples=True, add_count=False, tablefmt="html")}
                    <h3>Namespaces</h3>
                    {ss.namespaces_str(graph, examples=True, add_count=False, tablefmt="html")}
                    <h3>Edges</h3>
                    {ss.edges_str(graph, examples=True, add_count=False, tablefmt="html")}""",
                    unsafe_allow_html=True,
                )

            kgg_utils.disease_figures(
                disease_name=st.session_state["user_disease"],
                graph=st.session_state["graph"],
            )
            if "graph" in st.session_state and st.session_state["graph"] is not None:
                zip_data = kgg_utils.create_zip()
                folder_name = f"{st.session_state.disease_name}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}.zip"

                st.download_button(
                    label="Download all files",
                    data=zip_data,
                    file_name=folder_name,
                    mime="application/zip",
                )

            else:
                st.error("No graph found. Please generate the graph first.")

            if st.button("Update threshold score"):
                delete_old_cache_and_files()
                st.session_state["graph"].nodes["protein"]["threshold_score"] = score
                st.success(f"Threshold score updated to {score}.")
                st.session_state["graph"].nodes["protein"]["threshold_score"] = score
                rv_base, rv_stats = kgg_utils.get_graph_summary(
                    st.session_state["graph"]
                )
                with st.container():
                    st.markdown(
                        f"""\
                        <h3>Metadata</h3>
                        {tabulate(rv_base, tablefmt="html")}
                        <h3>Statistics</h3>
                        {tabulate(rv_stats, tablefmt="html")}
                        <h3>Nodes</h3>
                        {ss.functions_str(graph, examples=True, add_count=False, tablefmt="html")}
                        <h3>Namespaces</h3>
                        {ss.namespaces_str(graph, examples=True, add_count=False, tablefmt="html")}
                        <h3>Edges</h3>
                        {ss.edges_str(graph, examples=True, add_count=False, tablefmt="html")}""",
                        unsafe_allow_html=True,
                    )
                    zip_data = kgg_utils.create_zip()
                    folder_name = f"{st.session_state.disease_name}.zip"
                    st.download_button(
                        label="Download updated files",
                        data=zip_data,
                        file_name=folder_name,
                        mime="application/zip",
                    )

        if "start_over" not in st.session_state:
            st.session_state.start_over = False
        if st.button("Start over"):
            delete_old_cache_and_files()
            for key in st.session_state.keys():
                #                if key.startswith("graph_") or key.startswith("dis2prot_"):
                del st.session_state[key]
            st.session_state["user_disease"] = disease_name
            st.session_state["disease_id"] = disease_id
            st.session_state["disease_name"] = disease_df[
                disease_df["id"] == disease_id
            ]["name"].values[0]
            st.session_state["ct_phase"] = ct_phase

            st.session_state.button_clicked = False
            st.session_state.start_over = True
            st.rerun()


with tab3:
    st.markdown("### Drug-likeness profile and beyond")
    st.markdown(
        """ This page allows users to get insights of chemicals and drugs associated with diseases, generated by the KGG. Alternatively, you can import 
         any file that has chemicals with SMILES representations.
        """
    )
    # File uploader
    uploaded_file = st.file_uploader(
        "Choose the CSV file named 'diseaseAssociatedDrugs' from a KGG output folder",
        type="csv",
    )

    if uploaded_file is not None:
        # Read the CSV file
        df = pd.read_csv(uploaded_file, index_col=0)
        # st.write(len(df))

        st.write(
            f"Total number of unique drugs is {len(set(df['drugId']))}. Please remember that a same drug can be in different phases of clinical trials."
        )

        calc_filters, unusedDrugs_df = kgg_utils.calculate_filters(df, "drugId")

        st.write(calc_filters)

        calc_filters = calc_filters.to_csv(index=False).encode("utf-8")

        unusedDrugs_df = unusedDrugs_df.to_csv(index=False).encode("utf-8")

        # st.write(unusedDrugs_df)

        st.download_button(
            label="Download druglikeness profile",
            data=calc_filters,
            file_name="druglikeness_df.csv",
            mime="text/csv",
        )

        st.write(
            "Some drugs may not have SMILES representation because their type is either antibody, protein or unknown. The unparsed drugs can be downloaded here."
        )

        st.download_button(
            label="Download unparsed drugs file",
            data=unusedDrugs_df,
            file_name="unparsedDrugs_df.csv",
            mime="text/csv",
        )


# footer with text and green background
st.markdown(
    "<footer style='background-color: #149372; padding: 10px; border-radius: 10px;'>"
    "<p style='color: white; text-align: center;'>Fraunhofer ITMP © 2024</p>"
    "<p style='color: white; text-align: center;'>This work has been conducted across several key projects in which ITMP has been actively involved.</p>"
    "</footer>",
    unsafe_allow_html=True,
)
