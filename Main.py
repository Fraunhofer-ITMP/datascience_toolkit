"""Landing page for all Fraunhofer apps."""

import os
import streamlit as st
import base64
from streamlit_card import card

st.set_page_config(layout="wide")

# Customize sidebar
markdown = """
**Info**: This is a collection of data science tools developed by Fraunhofer ITMP SreeningPort.

**External links**:
* [GitHub](https://github.com/Fraunhofer-ITMP/streamlit_ITMP)
* [Website](https://www.itmp.fraunhofer.de/)
"""

st.sidebar.title("About")
st.sidebar.markdown(markdown)
st.sidebar.image("images/fraunhofer_ITMP-logo_900p.jpg")

st.markdown(
    "<h1 style='text-align: center; color: #149372;'>Fraunhofer ITMP ScreeningPort <br> Data Science Toolkits</h1>",
    unsafe_allow_html=True,
)

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))


def load_logo(filepath):
    with open(filepath, "rb") as f:
        data = f.read()
        encoded = base64.b64encode(data)
    data = "data:image/png;base64," + encoded.decode("utf-8")
    return data


# Add cards for each app
#col1, col2, col3, col4, col5, col6 = st.columns(6)

col1, col2  = st.columns(2)

col3, col4  = st.columns(2)

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

with col4:
    hasClicked = card(
        title="KGG_user",
        text="Test",
        image=load_logo("images/clinical_ds.png"),
        url="/kggUser",
    )

# with col5:
#     hasClicked = card(
#         title="KGG workflow test 1",
#         text="Test2Jupyter",
#         image=load_logo("images/kg_2.png"),
#         url="/kgg_user_2",
#     )

# with col6:
#     hasClicked = card(
#         title="KGG workflow test 2",
#         text="Test",
#         image=load_logo("images/kg_1.png"),
#         url="/kgg_user_3",
#     )
