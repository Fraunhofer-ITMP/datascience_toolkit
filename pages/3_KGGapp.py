"""App for Knowledge Graph Generator (KGG)"""

import streamlit as st
import streamlit.components.v1 as components
# from PIL import Image
# import pandas as pd
# import numpy as np
# from bokeh.plotting import figure
# from bokeh.models import ColumnDataSource, HoverTool
# from rdkit import Chem
# from rdkit.Chem import Draw
# import io
# import base64
# from openpyxl import Workbook
# from openpyxl.utils.dataframe import dataframe_to_rows
import os
import glob
#from tqdm import tqdm
from utils import *
#from kg_gen_4 import *

def load_kg(path):
    infile = open(path, 'rb')
    kg = pickle.load(infile)
    infile.close()
    return(kg)

def disease_figures(disease_name):

    st.write('You selected:',disease_name)

    if disease_name == 'Alzheimer':
        disName = 'ad_'
    if disease_name == 'Parkinson':
        disName = 'pd_'
    if disease_name == 'Depression':
        disName = 'dep_'
    if disease_name == 'COVID-19':
        disName = 'covid_'

    path = os.getcwd()

    #locate kg pkl file with disease sub-string match on filename
    kg = glob.glob(path + '\\kgs\\' + disName +'*pkl')

    kg = load_kg(kg[0])
    node_num = kg.number_of_nodes()
    edge_num = kg.number_of_edges()

    st.write('The ', option, ' KG consists of ', node_num, ' nodes and ', edge_num, ' edges.')

    st.write('Now let\'s take a look into various figures and tables related to', option)

    fig_summary = './images/kg_summary_demo.png'
    fig_ns = glob.glob(path + '\\kgs\\' + disName + 'namespace.png')[0]
    fig_ct = glob.glob(path + '\\kgs\\' + disName + 'CTphase.png')[0]
    fig_dt = glob.glob(path + '\\kgs\\' + disName + 'drugType.png')[0]
    fig_dlk = glob.glob(path + '\\kgs\\' + disName + 'drugLikeness.png')[0]

    st.subheader('1. Types and numbers of KG entities and their relationships')
    st.image(fig_summary)

    with open(fig_summary, "rb") as file:
        btn = st.download_button(
            label="Download image",
            data=file,
            file_name="kg_summary.png",
            mime="image/png",
        )

    st.subheader('2. Summary of namespaces')
    st.image(fig_ns)

    with open(fig_ns, "rb") as file:
        btn = st.download_button(
            label="Download image",
            data=file,
            file_name="namespace.png",
            mime="image/png",
        )

    st.subheader('3. Summary of namespaces (Interactive Version)')

    path_ns = glob.glob(path + '\\kgs\\' + disName + 'namespace.html')[0]
    HtmlFile = open(path_ns, 'r', encoding='utf-8')
    source_code = HtmlFile.read()
    #hit and trial for height adjustment == 400
    components.html(source_code,height=450,width=800)

    st.subheader('4. Pie charts for drug types and their clinical trial phases')
    # align images one after another
    img_1, img_2 = st.columns(2)

    img_1.image(fig_dt, use_column_width=True)

    with open(fig_dt, "rb") as file:
        btn = st.download_button(
            label="Download image",
            data=file,
            file_name="drugType.png",
            mime="image/png",
        )

    img_2.image(fig_ct, use_column_width=True)

    # image download button
    with open(fig_ct, "rb") as file:
        btn = st.download_button(
            label="Download image",
            data=file,
            file_name="CTphase.png",
            mime="image/png",
        )

    st.subheader('5. Pie chart for types of drugs')

    path_dt = glob.glob(path + '\\kgs\\' + disName + 'drugType.html')[0]
    HtmlFile = open(path_dt, 'r', encoding='utf-8')
    source_code = HtmlFile.read()
    # hit and trial for height adjustment == 400
    components.html(source_code, height=400, width=800)

    st.subheader('6. Pie chart for clinical trial phases of drugs')

    path_ct = glob.glob(path + '\\kgs\\' + disName + 'CTphase.html')[0]
    HtmlFile = open(path_ct, 'r', encoding='utf-8')
    source_code = HtmlFile.read()
    # hit and trial for height adjustment == 400
    components.html(source_code, height=400, width=800)

    st.subheader('7. Assessment of druglikeness using various filters')

    st.image(fig_dlk)

    with open(fig_dlk, "rb") as file:
        btn = st.download_button(
            label="Download image",
            data=file,
            file_name="druglikeness.png",
            mime="image/png",
        )

    st.subheader('8. Interactive parallel co-ordinate plot for KG drugs')

    path_pcp = glob.glob(path + '\\kgs\\' + disName + '*PCP*')[0]
    HtmlFile = open(path_pcp, 'r', encoding='utf-8')
    source_code = HtmlFile.read()
    # hit and trial for height adjustment == 400
    components.html(source_code, height=400)


st.set_page_config(layout="wide")

st.markdown(
    "<h1 style='text-align: center; color: #149372;'>Knowledge Graph Generator (KGG)</h1>",
    unsafe_allow_html=True,
)

st.subheader('This page provides various figures and plots for diseases that were depicted in the KGG manuscript')

url_git = 'https://github.com/Fraunhofer-ITMP/kgg'
url_idkit = 'https://www.infectious-diseases-toolkit.org/showcase/knowledge-graph-generator'
st.write('KGG is a fully automated workflow for creating disease-specific KGs. '
             'It is developed for a broad spectrum of researchers and scientists, especially for those who are into pre-clinical drug discovery, understanding disease mechanisms/comorbidity and drug-repurposing. '
             'The KGG embeds underlying schema of curated public databases to retrieve relevant knowledge which is regarded as the gold standard for high quality data. '
         'The KGG is leveraged on our previous contributions to the BY-COVID project where we developed workflows for identification of bio-active analogs for fragments identified in COVID-NMR studies (Berg, H et al.) and representation of Mpox biology (Karki, R et al.). '
         'The programmatic scripts and methods for KGG are written in python (version 3.10) and are available on [GitHub](%s).' % url_git, 'Additional information about KGG can be found [here](%s).' % url_idkit)


kg_worklow = './images/KGG_workflow_updated.png'
st.image(kg_worklow, caption = 'A schematic representation of the 3 phases of KGG workflow')

option = st.selectbox('Please choose a disease by clicking on the grey area below',('None','Alzheimer','Parkinson','Depression','COVID-19'))

if option == 'None':

    #st.write('Waiting for a selection!')
    components.html('<iframe src="https://giphy.com/embed/MgRKCBGvlpqTENUzWk" width="480" height="346" style="" frameBorder="0" class="giphy-embed" allowFullScreen></iframe><p><a href="https://giphy.com/gifs/fallontonight-jimmy-fallon-hurry-up-tonightshow-MgRKCBGvlpqTENUzWk">via GIPHY</a></p>',height = 400)

if option == 'Alzheimer':
    disease_figures(option)

if option == 'Parkinson':
    disease_figures(option)

if option == 'Depression':
    disease_figures(option)

if option == 'COVID-19':
    disease_figures(option)

# if option == 'Alzheimer':
#
#     st.write('You selected:', option)
#
#     kg = load_kg('./kgs/ad_aug_2024.pkl')
#     node_num = kg.number_of_nodes()
#     edge_num = kg.number_of_edges()
#
#     # st.write(kg.number_of_nodes())
#
#     st.write('The ', option, ' KG consists of ', node_num, ' nodes and ', edge_num, ' edges.')
#
#     st.write('Now let\'s take a look into various figures and tables related to', option)
#
#     fig_summary = './images/kg_summary_demo.png'
#     #fig_ns = './images/ad_aug_2024_namespace.png'
#     fig_ns = './images/ad_namespace.png'
#     fig_ct = './images/ad_CTphase.png'
#     fig_dt = './images/ad_drugType.png'
#     fig_dlk = './images/ad_druglikeness.png'
#
#     st.subheader('1. Types and numbers of KG entities and their relationships')
#     st.image(fig_summary)
#
#     with open(fig_summary, "rb") as file:
#         btn = st.download_button(
#             label="Download image",
#             data=file,
#             file_name="kg_summary.png",
#             mime="image/png",
#         )
#
#     st.subheader('2. Summary of namespaces')
#     st.image(fig_ns)
#
#     with open(fig_ns, "rb") as file:
#         btn = st.download_button(
#             label="Download image",
#             data=file,
#             file_name="namespace.png",
#             mime="image/png",
#         )
#
#     st.subheader('3. Summary of namespaces (Interactive Version)')
#
#     HtmlFile = open("images/ad_namespace.html", 'r', encoding='utf-8')
#     source_code = HtmlFile.read()
#     #hit and trial for height adjustment == 400
#     components.html(source_code,height=450,width=800)
#
#
#     st.subheader('4. Pie charts for drug types and their clinical trial phases')
#     # align images one after another
#     img_1, img_2 = st.columns(2)
#
#     img_1.image(fig_dt, use_column_width=True)
#
#     with open(fig_dt, "rb") as file:
#         btn = st.download_button(
#             label="Download image",
#             data=file,
#             file_name="drugType.png",
#             mime="image/png",
#         )
#
#     img_2.image(fig_ct, use_column_width=True)
#
#     # image download button
#     with open(fig_ct, "rb") as file:
#         btn = st.download_button(
#             label="Download image",
#             data=file,
#             file_name="CTphase.png",
#             mime="image/png",
#         )
#
#     st.subheader('5. Pie chart for types of drugs')
#
#     HtmlFile = open("images/ad_drugType.html", 'r', encoding='utf-8')
#     source_code = HtmlFile.read()
#     #hit and trial for height adjustment == 400
#     components.html(source_code,height=400,width=800)
#
#     st.subheader('6. Pie chart for clinical trial phases of drugs')
#     HtmlFile = open("images/ad_CTphase.html", 'r', encoding='utf-8')
#     source_code = HtmlFile.read()
#     #hit and trial for height adjustment == 400
#     components.html(source_code,height=400,width=800)
#
#     st.subheader('7. Assessment of druglikeness using various filters')
#
#     st.image(fig_dlk)
#
#     with open(fig_dlk, "rb") as file:
#         btn = st.download_button(
#             label="Download image",
#             data=file,
#             file_name="druglikeness.png",
#             mime="image/png",
#         )
#
#     st.subheader('8. Interactive parallel co-ordinate plot for KG drugs')
#
#     HtmlFile = open("images/ad_druglikeness_PCP.html", 'r', encoding='utf-8')
#     source_code = HtmlFile.read()
#     #hit and trial for height adjustment == 400
#     components.html(source_code,height=400)



    #images = [img_1, img_2]
    #st.image(images, use_column_width=True, caption=["some generic text"] * len(images))

    #st.html('images/ad_drugType.html', unsafe_allow_html=True)


