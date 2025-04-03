"""App for Proxidrugs"""

import pandas as pd
import streamlit as st
import json
from urllib.request import urlopen

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

st.header("Description", anchor="description", divider="gray")
st.write(
    "PROXIDRUGS aims at improving this novel class of drugs in multiple ways: the strategy comprises the identification of suitable functional subunits for molecule engineering; the directed optimization of pharmacological properties; and the transfer to clinical phases. In line with this, the cluster interlinks partners from academia and industry, covering the entire value chain from basic to translational research. Academic participants are based at Goethe University Frankfurt, TU Darmstadt, University of Heidelberg, the MPI of Biophysics, and the Fraunhofer Institute for Translational Medicine and Pharmacology. They are closely working together with major pharma partners, namely AbbVie Deutschland, Merck Healthcare, and GlaxoSmithKline, as well as Revvity as a technology provider."
)

st.write(
    "This dashboard provides information about partners and individuals involved in sub-projects of PROXIDRUGS."
)

st.header(
    "Our partners and their locations",
    divider="gray",
    help="This section allows you know more about the PROXIDRUGS project and see the basic information surrounding the sub-projects",
)


@st.cache_data
def build_map():
    # Get map data - https://raw.githubusercontent.com/eurostat/Nuts2json/master/pub/v2/2021/4326/20M/nutsrg_1.json
    with open("./data/proxidrugs_data/europe.json", "r") as f:
        map_germany = json.load(f)

    # state level
    map_data = pd.read_csv(
        "./data/proxidrugs_data/partner_location.csv", encoding="utf-8-sig"
    )

    custom_color_map = {
        "Hessen": "Yellow",
        "Hamburg": "Red",
        "London": "Green",
        "Baden-Württemberg": "Blue",
        "Rheinland-Pfalz": "Orange",
    }
    default_color = "#CCCCCC"  # Grey

    # Create a new column for coloring
    map_data["color"] = map_data["name"].apply(
        lambda x: x if x in custom_color_map else "Other"
    )

    custom_color_map["Other"] = default_color
    # Geographic Map
    fig = px.choropleth_mapbox(
        map_data,
        geojson=map_germany,
        locations="name",
        color="color",
        color_discrete_map={**custom_color_map},
        mapbox_style="open-street-map",
        zoom=4,
        center={"lat": 51.0057, "lon": 13.7274},
        opacity=0.5,
        labels={"partner": "Institute", "counts": "People"},
        featureidkey="properties.na",
        hover_name="name",
        hover_data={"counts": True, "partner": True, "name": False, "color": False},
    )
    fig.update_layout(
        showlegend=False,
        margin={"r": 0, "t": 0, "l": 0, "b": 0},
        height=500,
        hovermode="closest",
    )
    st.plotly_chart(fig, use_container_width=True)


build_map()

st.header("Overview of individuals", divider="gray")

prx = pd.read_csv("./data/proxidrugs_data/proxidrugs_partners.csv")

col = st.columns(2, gap="large")
with col[0]:
    prx_projects = prx.groupby("Project")["Project"].count().reset_index(name="count")
    fig = px.pie(prx_projects, values="count", names="Project")
    fig.update_layout(title_text="Proportion of individuals per project", title_x=0.18)
    st.plotly_chart(fig)

with col[1]:
    st.write("**Data champions for each sub-project**")
    data_champ = pd.read_excel("./data/proxidrugs_data/data_champion.xlsx")
    st.dataframe(data_champ, use_container_width=True, hide_index=True)

st.header("Overview of partners", divider="gray")

# Create a dropdown menu for institute selection
institutes = sorted(prx["Institute/Company"].unique())
selected_institute = st.selectbox("Select an Institute/Company", institutes)

# Filter the dataframe based on the selected institute
filtered_df = prx[prx["Institute/Company"] == selected_institute]

# Count the number of people per project
project_counts = filtered_df["Project"].value_counts().to_frame().reset_index()

col = st.columns(2, gap="large")
with col[0]:
    # Create a pie chart using Plotly
    fig = px.pie(
        values=project_counts["count"],
        names=project_counts["Project"],
        title=f"Project Distribution for {selected_institute}",
        labels={"names": "Project", "values": "Number of People"},
    )

    # Update the hover template to use the new labels
    fig.update_traces(
        hovertemplate="<b>%{label}</b><br>Number of People: %{value}<extra></extra>"
    )
    st.plotly_chart(fig)

with col[1]:
    st.write("**Number of individuals per sub-project**")
    st.dataframe(project_counts, use_container_width=True, hide_index=True)

    st.write("**Need to contact the partner?**")
    csv_file = filtered_df.to_csv(index=False).encode("utf-8")
    st.download_button(
        "Press to Download",
        csv_file,
        f"{institutes}_file.csv",
        "text/csv",
        key="download-csv",
    )

st.header("Overview of sub-projects", divider="gray")

# Create a dropdown menu for project selection
projects = sorted(prx["Project"].unique())
selected_project = st.selectbox("Select a Project", projects)

# Filter the dataframe based on the selected project
project_filtered_df = prx[prx["Project"] == selected_project]

# Count the number of people per institute
institute_counts = (
    project_filtered_df["Institute/Company"].value_counts().to_frame().reset_index()
)

col = st.columns(2)

with col[0]:
    # Create a pie chart using Plotly
    fig = px.pie(
        values=institute_counts["count"],
        names=institute_counts["Institute/Company"],
        # title=f"Institute(s) in {selected_project}",
        labels={"names": "Institute/Company", "values": "Number of People"},
    )

    fig.update_layout(title_text=f"Institute(s) in {selected_project}", title_x=0.25)

    # Update the hover template
    fig.update_traces(
        hovertemplate="<b>%{label}</b><br>Number of People: %{value}<extra></extra>"
    )

    # Display the pie chart
    st.plotly_chart(fig)

with col[1]:
    # Optional: Display the filtered dataframe
    st.write("**Details about individuals**")
    project_filtered_df["name"] = (
        project_filtered_df["First Name"] + " " + project_filtered_df["Last Name"]
    )
    subset = project_filtered_df[["name", "Institute/Company", "Email"]]
    st.dataframe(subset, use_container_width=True, hide_index=True)
