import pickle
import pandas as pd
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
from stqdm import stqdm
import datetime
import streamlit as st
import streamlit.components.v1 as components
import io
import xlsxwriter
import json
from urllib.request import urlopen

#import plotly #pip install -U kaleido required
import plotly.graph_objects as go
import plotly.express as px

st.set_page_config(
    layout="wide",
    page_title="PROXIDRUGS Dashboard",
    page_icon=":pill:",
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
    "<h1 style='text-align: center; color: #149372;'> PROXIDRUGS Dashboard</h1> <br>",
    unsafe_allow_html=True,
)

st.markdown('#### Description')
st.write('PROXIDRUGS aims at improving this novel class of drugs in multiple ways: the strategy comprises the identification of suitable functional subunits for molecule engineering; the directed optimization of pharmacological properties; and the transfer to clinical phases. In line with this, the cluster interlinks partners from academia and industry, covering the entire value chain from basic to translational research. Academic participants are based at Goethe University Frankfurt, TU Darmstadt, University of Heidelberg, the MPI of Biophysics, and the Fraunhofer Institute for Translational Medicine and Pharmacology. They are closely working together with major pharma partners, namely AbbVie Deutschland, Merck Healthcare, and GlaxoSmithKline, as well as Revvity as a technology provider.')

st.write("This dashboard provides information about partners and individuals involved in sub-projects of PROXIDRUGS.")

# st.header(
#     "PROXIDRUGS partners in Germany",
#     divider="gray",
#     help="This section allows you know more about the PROXIDRUGS project and see the basic information surrounding the sub-projects",
# )

st.markdown('#### PROXIDRUGS partners and their locations')

# Get map data
with urlopen(
    #map of EU
    "https://raw.githubusercontent.com/eurostat/Nuts2json/master/pub/v2/2021/4326/20M/nutsrg_1.json"
    #state level
    #"https://raw.githubusercontent.com/isellsoap/deutschlandGeoJSON/main/2_bundeslaender/1_sehr_hoch.geo.json"
    # #region level
    # "https://raw.githubusercontent.com/isellsoap/deutschlandGeoJSON/refs/heads/main/4_kreise/1_sehr_hoch.geo.json"
) as response:
    map_eu = json.load(response)

#state level
map_data = pd.read_csv("./data/proxidrugs_data/partner_location.csv",encoding='utf-8-sig')


custom_color_map = {

    "Hessen": "Yellow",
    "Hamburg": "Red",
    "London" : "Green",
    "Baden-Württemberg" : "Blue",
    "Rheinland-Pfalz" : "Orange",

}

default_color = "#CCCCCC"  # Light gray

# Create a new column for coloring
map_data['color'] = map_data['name'].map(lambda x: x if x in custom_color_map else 'Other')

# Geographic Map
fig = px.choropleth_mapbox(
    map_data,
    geojson=map_eu,
    locations="name",
    color="color",
    #color_discrete_map=custom_color_map,
    color_discrete_map={**custom_color_map, "Other": default_color},
    #color_continuous_scale="Viridis",
    mapbox_style="open-street-map",
    zoom=4,
    center={"lat": 51.0057, "lon": 13.7274},
    opacity=0.5,
    labels={"partner":"Institute","counts": "People"},
    hover_data=["counts", "partner"],
    featureidkey="properties.na",
)
#fig.update_layout(coloraxis_showscale=False)
fig.update_layout(showlegend=False,margin={"r": 0, "t": 0, "l": 0, "b": 0},height=500)
# Update color axis to show discrete legend
#fig.update_layout(coloraxis_showscale=False)

st.plotly_chart(fig, use_container_width=True)

st.write("")

st.markdown('#### Overview of individuals per sub-project')

col = st.columns(2,gap='large')

with col[0]:

    #st.write("**Overview of number of individuals per sub-project**")

    prx = pd.read_csv('./data/proxidrugs_data/proxidrugs_partners.csv')

    prx = prx.groupby(['Project'])['Project'].count().reset_index(name='count')
    fig = px.pie(prx, values='count', names='Project')
    fig.update_layout(title_text='Proportion of individuals per project', title_x=0.18)

    st.plotly_chart(fig)

    # path_ns = './images/subProject.html'
    # HtmlFile = open(path_ns, 'r', encoding='utf-8')
    # source_code = HtmlFile.read()
    # # hit and trial for height adjustment == 400
    # components.html(source_code, height=400, width=600)

with col[1]:

    st.write("**Data champions for each sub-project**")
    data_champ = pd.read_excel('./data/proxidrugs_data/data_champion.xlsx')
    st.write(data_champ)

st.markdown("#### Overview of partners")

df = pd.read_csv('./data/proxidrugs_data/proxidrugs_partners.csv', encoding='utf-8-sig')
# st.write(df)

# Create a dropdown menu for institute selection
institutes = sorted(df['Institute/Company'].unique())
selected_institute = st.selectbox("Select an Institute/Company", institutes)

# Filter the dataframe based on the selected institute
filtered_df = df[df['Institute/Company'] == selected_institute]

# Count the number of people per project
project_counts = filtered_df['Project'].value_counts()
# project_counts = pd.DataFrame(project_counts,columns=['Name','Number'])
# data = dict(project_counts)


col = st.columns(2,gap='large')

with col[0]:
    # Create a pie chart using Plotly
    fig = px.pie(
        values=project_counts.values,
        names=project_counts.index,
        # title=f"Project Distribution for {selected_institute}",
        labels={'names': 'Project', 'values': 'Number of People'}
    )

    # Update layout
    fig.update_layout(
        title_text="Proportion of partners per each sub-project in PROXIDRUGS"
        # Add any other layout customizations here
    )
    # Update the hover template to use the new labels
    fig.update_traces(hovertemplate='<b>%{label}</b><br>Number of People: %{value}<extra></extra>')

    # Display the pie chart
    st.plotly_chart(fig)


with col[1]:
    st.write("**Number of individuals per sub-project**")

    st.write(project_counts)


# Optional: Display the filtered dataframe
st.write("**Details about individuals**")
st.dataframe(filtered_df)

st.markdown("#### Overview of sub-projects")

# Load your data
#df = pd.read_csv('your_data.csv')  # Replace with your actual data source

# Create the Streamlit app
#st.title("Institute Distribution by Project")

# Create a dropdown menu for project selection
projects = sorted(df['Project'].unique())
selected_project = st.selectbox("Select a Project", projects)

# Filter the dataframe based on the selected project
filtered_df = df[df['Project'] == selected_project]

# Count the number of people per institute
institute_counts = filtered_df['Institute/Company'].value_counts()

col = st.columns(2)

with col[0]:

    # Create a pie chart using Plotly
    fig = px.pie(
        values=institute_counts.values,
        names=institute_counts.index,
        #title=f"Institute(s) in {selected_project}",
        labels={'names': 'Institute/Company', 'values': 'Number of People'}
    )

    fig.update_layout(title_text=f"Institute(s) in {selected_project}", title_x=0.25)

    # Update the hover template
    fig.update_traces(hovertemplate='<b>%{label}</b><br>Number of People: %{value}<extra></extra>')

    # Display the pie chart
    st.plotly_chart(fig)

with col[1]:
    # Optional: Display the filtered dataframe
    st.write("**Details about individuals**")
    st.dataframe(filtered_df)




