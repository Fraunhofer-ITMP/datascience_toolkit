"""Display Page for OpenScreen Impulse"""

import streamlit as st
from PIL import Image

st.set_page_config(
    layout="wide",
    page_title="EU OpenScreen-Impulse",
    page_icon="🇪🇺",
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
    "<h1 style='text-align: center; color: #FF8E1C;'> EU OpenScreen-Impulse </h1> <br>",
    unsafe_allow_html=True,
)
st.markdown(
    "<h2 style='text-align: center; color: #1B3C67;'> European Infrastructure of Open Screening Platforms for Chemical Biology </h2> <br>",
    unsafe_allow_html=True,
)
bg_image = Image.open("images/pages/openscreen_impulse_background.png")
resized_bg = bg_image.resize((1800, 300))
st.image(resized_bg, use_container_width=True)

st.markdown("## IMPULSE Initiative")
st.markdown(
    "**IMPULSE** seeks to enhance EU-OPENSCREEN's role as a premier platform for chemical biology and early drug discovery in Europe "
    "with regards to the key areas of **enhancing service catalogue**, **improving standards**, **fueling community engagement**, and **ensuring sustainability**"
)
st.markdown(
    "Additional information about **EU-OpenScreen IMPULSE** can be found [here.](https://www.eu-openscreen.eu/impulse-overview.html)"
)
st.markdown(
    """
    ### Key Focus Areas

    **Validating pharmacology using advanced disease models & genetic screens:** Open call for proposal to develop new models in uncovered disease areas or improve existing models using better assay technology.

    **New chemical modalities:** Demonstrator projects started, covering fluorescent probes, traceless proximity probes and quenched activity-based probes, targeted proteins/RNA degraders, chemically diverse covalent and allosteric protein modulators and multitarget compounds.

    **Data Reproducibility & Operational Standards:** Mapping of data standard across EU-OPENSCREEN partner sites ongoing. Preparation of an interlaboratory comparison to assess and improve reproducibility and consistency across EU-OPENSCREEN screening partners.

    **ML/AI for prediction of modes of action**  with cell painting and small molecule data.
    """
)

st.markdown("### Cataloguing Screening")
st.image("images/pages/openscreen_cataloguing_screening.png")
st.caption("Cataloguing Screening Visualized based on Site, Types and Assays.")
