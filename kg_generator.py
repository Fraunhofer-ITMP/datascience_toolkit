"""Display Page for KG Generator under KGG"""

import datetime
import uuid
import pandas as pd
import streamlit as st
import streamlit_antd_components as sac
from pybel.struct.summary import supersummary as ss
from tabulate import tabulate

import kgg_utils

state = st.session_state

st.set_page_config(
    layout="wide",
    page_title="KGG-Knowledge Graph Generator",
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
    index=1,  # highlighting "Drug-likeness assessment" (0=Description, 1=KG Generator, 2=Drug-likeness, 3=Visualizer)
    size="25",
    direction="horizontal",
)

if nav == "Description":
    st.markdown(
        """<meta http-equiv="refresh" content="0; url='https://fraunhofer-itmp-ds-toolkit.serve.scilifelab.se/KGG'" />""",
        unsafe_allow_html=True,
    )
elif nav == "KG Generator":
    pass
elif nav == "Drug-likeness assessment":
    st.markdown(
        """<meta http-equiv="refresh" content="0; url='https://druglikeness.serve.scilifelab.se/'" />""",
        unsafe_allow_html=True,
    )
elif nav == "KG Visualizer":
    st.markdown(
        """<meta http-equiv="refresh" content="0; url='https://kgviz.serve.scilifelab.se/'" />""",
        unsafe_allow_html=True,
    )

# ---------------- PAGE CONTENT ---------------- #
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

def disease_input_and_search():
    if "previous_disease" in state and state.user_disease != state.previous_disease:
        keys_to_clear = [
            "disease_df",
            "drugs_df",
            "dis2prot_df",
            "dis2snp_df",
            "graph",
        ]
        for key in keys_to_clear:
            if key in state:
                del state[key]
        state.previous_disease = state.user_disease

    if "user_disease" not in state:
        state["user_disease"] = ""

    disease_name = st.text_input(
        "Enter the disease you are interested in generating a graph.",
        placeholder="e.g. AIDS",
        key="user_disease",
    )

    if "disease_name" in state and state["user_disease"] != state["disease_name"]:
        for key in list(state.keys()):
            if key != "user_disease":
                del state[key]
        state.clear()
        state["user_disease"] = disease_name

    if disease_name.strip():
        df = kgg_utils.searchDisease(disease_name)

        if df.empty:
            st.warning("No results found for the disease. Please try again.")
        else:
            with st.expander("Disease search results"):
                st.dataframe(df, hide_index=True)

        return df
    else:
        st.info("Please enter a disease name to search.")
        return pd.DataFrame()

def disease_input_component():
    if "user_disease" not in state:
        state["user_disease"] = "AIDS"

    disease_name = st.text_input(
        "Enter the disease you are interested in generating a graph. Please hit the **Search for disease** button to search.",
        placeholder="e.g. AIDS",
        key="user_disease",
    )
    if disease_name == "":
        st.warning("Please enter a disease name.")
        st.stop()

def disease_search_component():
    if "disease_df" not in state:
        state["disease_df"] = pd.DataFrame()
    search_button = st.button(
        "Search for disease",
        key="search_button",
    )
    if search_button:
        state["disease_df"] = kgg_utils.searchDisease(state["user_disease"])
        if state["disease_df"].empty:
            st.write("No results found for the disease. Please try again.")
            st.stop()
        else:
            with st.expander("Disease search results"):
                st.dataframe(state["disease_df"], hide_index=True)

    return state["disease_df"]

#    disease_input_component()
#    disease_df = disease_search_component()
disease_df = disease_input_and_search()

col1, col2 = st.columns(2)

with col1:
    st.markdown(
        "**Please select the disease of interest from the table above by entering the identifier value.** \n Each disease has a unique id, which fetches associated proteins using OpenTargets API. Double click twice on the id of interest to select, copy, and paste the id into the text box."
    )
    if disease_df.empty:
        st.stop()
    else:
        #            st.write(f"Disease ID at this point: {disease_df.iloc[0]['id']}")
        if "disease_id" not in state:
            state["disease_id"] = disease_df.iloc[0]["id"]

        disease_id = st.text_input(
            "Please enter the identifier for disease of interest.",
            value=disease_df.iloc[0]["id"],
            key="disease_id",
        )
        
        disease_id_clean = disease_id.strip() if disease_id else ""
        valid_disease_ids = disease_df["id"].tolist()
        
        if disease_id_clean not in valid_disease_ids:
            st.error(f"Invalid disease ID: '{disease_id}'. Please select a valid ID from the search results table above.")
            st.stop()        
        if disease_id_clean != disease_id:
            state["disease_id"] = disease_id_clean

with col2:
    st.markdown(
        """**Please enter the clinical trial phase of chemicals which should be used by the workflow.** \n Use a number between 1 (early phase) and 4 (FDA approved). Keep in mind, lower input values increase the number of identified chemicals and running time."""
    )

    if "ct_phase" not in state:
        state["ct_phase"] = 3
    ct_phase = st.number_input(
        label="Select your clinical trial phase",
        min_value=1,
        max_value=4,
        value=state["ct_phase"],
        key="ct_phase",
    )

st.markdown(
    ":red[Selected disease:] "
    + str(state["user_disease"])
    + ":red[ with ID:] "
    + str(state["disease_id"])
    + ":red[ and clinical trial phase:] "
    + str(state["ct_phase"])
)

try:
    viral_proteins = kgg_utils.GetViralProteins(state["user_disease"])
    state["viral_prot"] = viral_proteins if viral_proteins is not None else []
except Exception as e:
    st.error(f"Error in viral protein identification: {str(e)}")
    state["viral_prot"] = []

# if "viral_prot" not in st.session_state:
#     st.session_state["viral_prot"] = viral_prot
# elif not st.session_state["viral_prot"] == viral_prot:
#     # elif viral_prot != st.session_state["viral_prot"]:
#     st.session_state["viral_prot"] = viral_prot

#    st.write(state)

st.header("Generating the graph", anchor="generate-graph", divider="grey")

dis_name = state["user_disease"].lower().replace(" ", "_").replace("-", "_")
kg_name = f"kgg_{dis_name}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
if "kg_name" not in state:
    state["kg_name"] = kg_name
elif kg_name != state["kg_name"]:
    state["kg_name"] = kg_name

if "button_clicked" not in state:
    state["button_clicked"] = False

def callback():
    state["button_clicked"] = True

if st.button("Generate Base Knowledge Graph", on_click=callback) or state.get(
    "button_clicked", False
):
    # st.write(state)
    #        st.write(state["user_disease"])
    #        st.write(state["disease_id"])
    state["drugs_df"] = pd.DataFrame()
    state["dis2prot_df"] = pd.DataFrame()
    state["dis2snp_df"] = pd.DataFrame()
    state["drugs_df"], state["dis2prot_df"], state["dis2snp_df"] = (
        kgg_utils.createInitialKG(
            _ct_phase=state["ct_phase"],
        )
    )

    #        state["dis2prot_df"] = dis2prot_df
    #        state["dis2snp_df"] = dis2snp_df
    #        if not state["dis2prot_df"].empty:
    #            st.write("Current protein data preview:")
    #            st.dataframe(state["dis2prot_df"].head())
    #            st.write(
    #                f"Current disease in session state: {st.session_state.get('user_disease')}"
    #            )
    kgg_utils.GetDiseaseAssociatedProteinsPlot(state["dis2prot_df"])

    def threshold_input_component():
        score_input = st.text_input(
            "Enter threshold score (recommended > 0.3):",
            placeholder="e.g., 0.3",
            key="protein_score_input",
        )

        if not score_input.strip():
            st.info("Please enter a threshold score and hit Enter.")
            st.stop()

        try:
            score = float(score_input)
            if not (0.0 <= score <= 1.0):
                st.warning("Please enter a score between 0.0 and 1.0.")
                st.stop()
        except ValueError:
            st.warning("Please enter a valid numeric value.")
            st.stop()

        return score

    def apply_threshold_and_update(score):
        if "filtered_df" not in state:
            state["filtered_df"] = pd.DataFrame()
        state["filtered_df"] = state["dis2prot_df"][
            state["dis2prot_df"]["Score"] >= score
        ]

        st.warning(
            f"""ℹ️ Filter of protein score (> {score}) has been applied. 
            This reduced the number of proteins from {len(state["dis2prot_df"])} to {len(state["filtered_df"])}."""
        )

        if "graph" not in state:
            state["graph"] = None
        state["graph"] = kgg_utils.finalizeKG(
            state["filtered_df"], session_inputs=state
        )

        st.header("Graph summary", anchor="graph-summary", divider="grey")
        if "rv_base" not in state:
            state["rv_base"] = ""
        if "rv_stats" not in state:
            state["rv_stats"] = ""
        state["rv_base"], state["rv_stats"] = kgg_utils.get_graph_summary(
            state["graph"]
        )

        with st.container():
            st.markdown(
                f"""\n<h3>Metadata</h3>
                {tabulate(state["rv_base"], tablefmt="html")}
                <h3>Statistics</h3>
                {tabulate(state["rv_stats"], tablefmt="html")}
                <h3>Nodes</h3>
                {ss.functions_str(state["graph"], examples=True, add_count=False, tablefmt="html")}
                <h3>Namespaces</h3>
                {ss.namespaces_str(state["graph"], examples=True, add_count=False, tablefmt="html")}
                <h3>Edges</h3>
                {ss.edges_str(state["graph"], examples=True, add_count=False, tablefmt="html")}""",
                unsafe_allow_html=True,
            )

    state["protein_score"] = threshold_input_component()
    apply_threshold_and_update(state["protein_score"])

    kgg_utils.disease_figures(
        disease_name=state["user_disease"],
        graph=state["graph"],
    )

#        state["button_clicked"] = False

if "graph" in state and state["graph"] is not None:
    if "zip_data" not in state:
        state["zip_data"] = ""
    if "folder_name" not in state:
        state["folder_name"] = ""
    state["zip_data"] = kgg_utils.create_zip()
    state["folder_name"] = (
        f"{state['user_disease']}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:4]}.zip"
    )
    #                st.balloons()

    st.download_button(
        label="Download all files",
        data=state["zip_data"],
        file_name=state["folder_name"],
        mime="application/zip",
        on_click="ignore",
        help="Click to download the zip file containing all graphs and CSVs. This will also download the CSV file *diseaseAssociatedDrugs.csv* that you can utilize for Drug-likeness assessment on the next tab.",
    )

    if st.button(
        "Start Over",
        help="This button usually takes a while. Please have patience.",
    ):
        #            kgg_utils.createInitialKG.clear()
        state.clear()
        state["button_clicked"] = False
        st.rerun()

else:
    st.warning("KG will be generated once all inputs are given.")
