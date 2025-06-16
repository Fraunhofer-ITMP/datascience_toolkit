"""App for CityCycling"""

import streamlit as st
import streamlit.components.v1 as components
import datetime
import base64

st.set_page_config(
    layout="wide",
    page_title="CityCycling",
    page_icon=":bike:",
    initial_sidebar_state="collapsed",
)

st.markdown(
    """
    <style>
        button[title^=Exit]+div 
    </style>
    """,
    unsafe_allow_html=True,
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
    "<h1 style='text-align: center; color: #149372;'> City Cycling 2025</h1> <br>",
    unsafe_allow_html=True,
)

st.markdown("**And.... we are back again representing :fast_forward: Fraunhofer ITMP Hamburg - Pedaling for Progress!!** :bike: :woman-biking: :man-biking: :right-facing_fist::left-facing_fist:")

st.markdown('**Cycling for a better climate and health:** Why not when it contributes to lesser emissions of CO<sub>2</sub> and improves your quality of life. Last year 15 of us covered ' \
'a total distance of **3007** kms from 376 rides. We were placed **247** out of ~1800 teams and contributed to avoiding 500 kg of CO<sub>2</sub> emission. The numbers in a bigger picture ' \
' are actually amazing. A total of 1.2 million cyclists covered a joint distance of 200 million kms thereby avoiding CO<sub>2</sub> emission of 35k tons. Wow!!! ' \
'This year 10k cyclists from Hamburg including **12** ITMPians will participate in the campaign. Our live data is displayed below (right).',unsafe_allow_html=True)

col1, col2 =st.columns(2)

gif_path = "./images/social_media/cycle.gif"

with col1:

    with open(gif_path, "rb") as f:
        contents = f.read()
        data_url = base64.b64encode(contents).decode("utf-8")

    st.markdown(
    f'<img src="data:image/gif;base64,{data_url}" alt="gif">',
    unsafe_allow_html=True,
)

with col2:


#<div style="display: flex; justify-content: center; align-items: center; width: 100%; min-height: 65vh;">

    itmp_rad = """

    <div style="width: auto !important; min-width: 375px; max-width: 415px; height: 415px;">
    <iframe style="width: 100%; height: 100%;" frameborder="0" scrolling="no" src="https://login.stadtradeln.de/specials/radelmeter/team/74776?L=2"></iframe>
    </div>"""

    st.markdown(itmp_rad,unsafe_allow_html=True)

# st.image(
#     './images/social_media/itmp_group.jpg',width=800
# )
st.write('')
st.write('')
st.write('Our group photo from 2024 :arrow_down:')

st.markdown(
    """
    <style>
        [data-testid=stImage]{
            display: block;
            margin-left: auto;
            margin-right: auto;
        }
    </style>
    """, unsafe_allow_html=True
)
st.image("./images/social_media/itmp_group.jpg")