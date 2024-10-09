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
        'Get Help': 'https://www.extremelycoolapp.com/help',
        'Report a bug': "https://www.extremelycoolapp.com/bug",
        'About': "# This is a header. This is an *extremely* cool app!"
    }
)

# Customize sidebar
# markdown = """
# **Info**: This is a collection of data science tools developed by Fraunhofer ITMP SreeningPort.

# **External links**:
# * [GitHub](https://github.com/Fraunhofer-ITMP/streamlit_ITMP)
# * [Website](https://www.itmp.fraunhofer.de/)
# """

# st.sidebar.title("About")
# st.sidebar.markdown(markdown)
# st.sidebar.image("images/fraunhofer_ITMP-logo_900p.jpg")

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

st.markdown("### About")
st.markdown(
    "This is a collection of data science tools developed by Fraunhofer ITMP SreeningPort."
)

tab1, tab2, tab3 = st.tabs(3)



col1, col2, col3  = st.columns(3)

# col3, col4  = st.columns(2)


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
        title="KGG",
        text="Automated workflow for disease-specific KGs",
        image=load_logo("images/kg_2.png"),
        url="/kgg",
    )

# with col4:
#     hasClicked = card(
#         title="KGG_user",
#         text="Test",
#         image=load_logo("images/kg_1.png"),
#         url="/kggUser",
#     )


# footer with text and green background
st.markdown(
    "<footer style='background-color: #006c8b; padding: 10px; border-radius: 10px;'>"
    "<p style='color: white; text-align: center;'>Fraunhofer ITMP © 2021</p>"
    "<p style='color: white; text-align: center;'>This work has been conducted across several key projects in which ITMP has been actively involved.</p>"
    "</footer>",
    unsafe_allow_html=True,
)
