"""Display Page for KG Visualizer under KGG"""
import datetime
import uuid

import pickle
import pandas as pd
import streamlit as st

import kgg_utils
import streamlit_antd_components as sac

st.set_page_config(
    layout="wide",
    page_title="KGG-KG Visualizer",
    page_icon="🕸️",
    initial_sidebar_state="collapsed",
)

st.markdown(
    "<h1 style='text-align: center; color: #149372;'> Knowledge Graph Generator (KGG)</h1> <br>",
    unsafe_allow_html=True,
)

nav = sac.buttons(
    [
        "Description",
        "KG Generator",
        "Drug-likeness assessment",
        "KG Visualizer",
    ],
    index=3,  # highlighting "Drug-likeness assessment" (0=Description, 1=KG Generator, 2=Drug-likeness, 3=Visualizer)
    size="25",
    direction="horizontal",
)

if nav == "Description":
    st.markdown(
        """<meta http-equiv="refresh" content="0; url='https://fraunhofer-itmp-ds-toolkit.serve.scilifelab.se/KGG'" />""",
        unsafe_allow_html=True,
    )
elif nav == "KG Generator":
    st.markdown(
        """<meta http-equiv="refresh" content="0; url='https://kggapp.serve.scilifelab.se/'" />""",
        unsafe_allow_html=True,
    )
elif nav == "Drug-likeness assessment":
    st.markdown(
        """<meta http-equiv="refresh" content="0; url='https://druglikeness.serve.scilifelab.se/'" />""",
        unsafe_allow_html=True,
    )
elif nav == "KG Visualizer":
    pass

# ---------------- PAGE CONTENT ---------------- #
st.markdown("### KG Visualizer")
st.markdown(
    """ This page allows users to get insights of a PyBEL Knowledge Graph...
    """
)
st.info(
    "**Important:** Please only upload the pickle file of a graph. Other types of files (CSV,XLSX) are not supported.",
    icon="ℹ️",
)

# File uploader
uploaded_file = st.file_uploader(
    "Please upload your pickle file.",
    type="pkl",
    
)

if uploaded_file is not None:
    query_graph = kgg_utils.load_pickle_file(uploaded_file=uploaded_file)
    if query_graph is None:
        st.error(
            "Uploaded file is either corrupt or the file format is wrong. Please upload the right file."
        )
        st.stop()
    kgg_utils.query_graph_info(query_graph)
    st.markdown("### Display of your graph:")
    graph_subset_html = kgg_utils.display_interactive_belgraph(query_graph)
    kgg_utils.download_interactive_belgraph(graph_subset_html)
elif uploaded_file is None:
    query_graph = pickle.load(
        open(
            "data/kgg_initial_loadfiles/drugRepurposing.pkl",
            "rb",
        )
    )
    st.markdown("Below is a sample graph for Parkinson's disease.")
    kgg_utils.query_graph_info(query_graph)
    st.markdown("### Display of the initial graph:")
    graph_subset_html = kgg_utils.display_interactive_belgraph(query_graph)
    kgg_utils.download_interactive_belgraph(graph_subset_html)


# footer with text and green background
current_year = datetime.datetime.now().year
st.markdown(
    f"<footer style='background-color: #149372; padding: 10px; border-radius: 10px;'>"
    f"<p style='color: white; text-align: center;'>Fraunhofer ITMP © {current_year}</p>"
    "<p style='color: white; text-align: center;'>This work has been conducted across several key projects in which ITMP has been actively involved.</p>"
    "</footer>",
    unsafe_allow_html=True,
)
