"""Display Page for OpenScreen Impulse"""

import streamlit as st
from PIL import Image
import pandas as pd
#import matplotlib.pyplot as plt
import plotly.express as px
from streamlit_agraph import agraph, Node, Edge, Config
import re
import plotly.graph_objects as go


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

    **🎯 Validating pharmacology using advanced disease models & genetic screens:** Open call for proposal to develop new models in uncovered disease areas or improve existing models using better assay technology.

    **🎯 New chemical modalities:** Demonstrator projects started, covering fluorescent probes, traceless proximity probes and quenched activity-based probes, targeted proteins/RNA degraders, chemically diverse covalent and allosteric protein modulators and multitarget compounds.

    **🎯 Data Reproducibility & Operational Standards:** Mapping of data standard across EU-OPENSCREEN partner sites ongoing. Preparation of an interlaboratory comparison to assess and improve reproducibility and consistency across EU-OPENSCREEN screening partners.

    **🎯 ML/AI for prediction of modes of action**  with cell painting and small molecule data.

    """,

    unsafe_allow_html=True,
)

# st.markdown("### Cataloging Screening")
# st.image("images/pages/openscreen_cataloguing_screening.png")
# st.caption("Cataloging Screening Visualized based on Site, Types and Assays.")

# df = pd.read_csv('data/Cataloguing screening.csv')

# df = df.loc[df['Site'].notna()]

# st.write(df)

tab1, tab2, tab3 = st.tabs(
    [
        "Cataloging Screening",
        "KG overview",
        "Data standard heatmap"
    ]
)

with tab1:

    st.header("Overview of partners", divider="gray")

    df = pd.read_csv('data/Cataloguing screening.csv')

    df = df.loc[df['Type'].notna()]

    # Create a dropdown menu for institute selection
    institutes = sorted(df["Site"].unique())
    selected_institute = st.selectbox("Select a partner institute", institutes,key='overview')

    # Filter the dataframe based on the selected institute
    filtered_df = df[df["Site"] == selected_institute]

    #st.dataframe(filtered_df)

    #st.dataframe(filtered_df[['Type','Model organism','Cell type','Assay format']], use_container_width=True, hide_index=True)

    # Count the number of people per project
    screen_counts = filtered_df["Type"].value_counts().to_frame().reset_index()

    #st.dataframe(screen_counts)

    col = st.columns(2, gap="large")
    with col[0]:
        # Create a pie chart using Plotly
        fig = px.pie(
            values=screen_counts["count"],
            names=screen_counts["Type"],
            title=f"Types of screening for {selected_institute}",
            labels={"names": "Type", "values": "Number of screens"},
        )

        # Update the hover template to use the new labels
        fig.update_traces(
            hovertemplate="<b>%{label}</b><br>Number of screens: %{value}<extra></extra>"
        )
        st.plotly_chart(fig)

    with col[1]:
        st.write("**Details of screening types**")
        st.dataframe(filtered_df[['Type','Model organism','Cell type','Assay format']], use_container_width=True, hide_index=True)
        # st.write("**Need to contact the partner?**")
        # csv_file = filtered_df.to_csv(index=False).encode("utf-8")
        # st.download_button(
        #     "Press to Download",
        #     csv_file,
        #     f"{institutes}_file.csv",
        #     "text/csv",
        #     key="download-csv",
        # )

    st.subheader("Overview of screening details")

    # # Prepare the summary counts (can adjust as needed)
    # pivot = df.groupby(['Site', 'Assay format']).size().unstack(fill_value=0)

    # # Create the matplotlib figure and axis
    # fig, ax = plt.subplots(figsize=(12, 6))
    # pivot.plot(kind='bar', stacked=False, ax=ax)
    # plt.title('Distribution of Assay Format Usage per Site')
    # plt.xlabel('Site')
    # plt.ylabel('Number of Experiments')
    # plt.xticks(rotation=45)
    # plt.legend(title='Assay format')
    # plt.tight_layout()

    # # Show the plot in Streamlit
    # st.pyplot(fig)


    # Load or receive your df here
    # For example:
    # df = pd.read_csv("your_data.csv")

    df = filtered_df[['Type','Model organism','Cell type','Assay format']]

    # Create a dropdown menu for institute selection
    df_cols = sorted(['Model organism','Cell type','Assay format'])
    #st.write(temp_col)

    #break 
    selected_column = st.selectbox("Select an attribute", df_cols)

    # Filter the dataframe based on the selected institute
    #filtered_df = df[df["Site"] == selected_institute]

    # Assuming df is already loaded with your data that has a 'Type' column

    type_counts = df[selected_column].value_counts().reset_index()
    type_counts.columns = [selected_column, 'Count']

    fig1 = px.bar(type_counts,
        x=selected_column,
        y='Count',
        title=f'Number of {selected_column.lower()}s per screening',
        labels={'Count': f'Number of {selected_column.lower()}s', 'Type': 'Screening Type'},
        text='Count',
        height=600)

    fig1.update_traces(textposition='outside')
    fig1.update_layout(uniformtext_minsize=8, uniformtext_mode='hide', xaxis_tickangle=-45)

    st.plotly_chart(fig1, use_container_width=True)

    grouped = df.groupby(['Type', selected_column]).size().reset_index(name='count')

    fig2 = px.bar(grouped,
        x='Type',
        y='count',
        color=selected_column,
        barmode='group',
        title=f'{selected_column} Distribution per Type',
        labels={'count': f'Number of {selected_column}s', 'Type': 'Type', selected_column: selected_column})

    fig2.update_layout(xaxis_tickangle=-45)

    st.plotly_chart(fig2, use_container_width=True)


    total_counts = df[selected_column].value_counts().reset_index()
    total_counts.columns = [selected_column, 'count']

    fig = px.pie(total_counts, 
                values='count', 
                names=selected_column, 
                title=f'Overall {selected_column} Distribution',
                hole=0.3)

    st.plotly_chart(fig, use_container_width=True)

with tab2:

    # Set Page Configuration
    st.set_page_config(layout="wide")

    st.title("IMPULSE partner overview")

    df = pd.read_csv('data/Cataloguing screening.csv')

    df = df.loc[df['Type'].notna()]

    cols = ['Site','Type','Cell type','Assay format','Model organism']

    df=df[cols]

    mappings = {'Site':'red','Type':'blue','Cell type':'green','Assay format':'purple','Model organism':'cyan'}

    nodes_2,edges_2 = [],[]

    node_unique_2 = []

    for item in cols:

        if item != 'Cell type':
            
            for obj in df[item]:

                if obj not in node_unique_2:
                    node_unique_2.append(obj) 
                    nodes_2.append(Node(id=obj,label=obj,color=mappings[item]))
                #st.write(Node(id=obj,label=obj,color=mappings[item]))
                #break
        #break

    #st.write(nodes)

    
    for site,type,cell,assay,org in df.values:

        edges_2.append(Edge(source=site,target=type))
        edges_2.append(Edge(source=type,target=assay))
        #edges.append(Edge(source=cell,target=assay))
        edges_2.append(Edge(source=assay,target=org))


    # Configure the graph
    config = Config(
        width=1500,
        height=900,
        directed=True,
        nodeHighlightBehavior=True,
        highlightColor="#F7A7A6",
        collapsible=False,
        node={'labelProperty': 'label'},
        link={'labelProperty': 'label', 'renderLabel': False},
        hierarchical=True,
        hierarchicalSorting=True,
        layout={
            "hierarchical": {
                "enabled": True,
                "levelSeparation": 150,
                "nodeSpacing": 100,
                "treeSpacing": 200,
                "direction": "UD",  # UD for top to bottom
                "sortMethod": "directed"
            }
        },
            # Set initial zoom level
        zoom=1.2  # Adjust as needed
    )


    # Display the graph
    agraph(nodes=nodes_2, edges=edges_2, config=config)

    st.subheader('Overview per partner')

    # Create a dropdown menu for institute selection
    partner = sorted(df["Site"].unique())
    selected_partner = st.selectbox("Select a partner institute", partner, key='kg')

    # Filter the dataframe based on the selected institute
    filtered_df = df[df["Site"] == selected_partner]

    nodes,edges = [],[]

    node_unique = []

    for item in cols:
        for obj in filtered_df[item]:

            if obj not in node_unique:
                node_unique.append(obj) 
                nodes.append(Node(id=obj,label=obj,color=mappings[item]))

    for site,type,cell,assay,org in filtered_df.values:

        edges.append(Edge(source=site,target=type,type="CURVE_SMOOTH"))
        edges.append(Edge(source=type,target=assay,type="CURVE_SMOOTH"))
        edges.append(Edge(source=assay,target=cell))
        edges.append(Edge(source=cell,target=org,type="CURVE_SMOOTH"))

    # Configure the graph
    config = Config(
        width=2000,
        height=1080,
        directed=True,
        nodeHighlightBehavior=True,
        highlightColor="#F7A7A6",
        collapsible=False,
        node={'labelProperty': 'label'},
        link={'labelProperty': 'label', 'renderLabel': False},
        hierarchical=False,
        hierarchicalSorting=True,
        layout={
            "hierarchical": {
                "enabled": True,
                "levelSeparation": 150,
                "nodeSpacing": 100,
                "treeSpacing": 200,
                "direction": "UD",  # UD for top to bottom
                "sortMethod": "directed"
            }
        },
            # Set initial zoom level
        zoom=1.2  # Adjust as needed
    )


    # Display the graph
    agraph(nodes=nodes, edges=edges, config=config)

with tab3:

    #st.header('Heatmap of IMPULSE survey')

    st.markdown('Following is a preliminary survey result that was conducted among the partners and institutes in the IMPULSE project. A detailed overview of the survery  will be provided soon.')

    data = pd.read_excel('data/Task 5.1 - Mapping data standards survey_cs.xlsx',sheet_name='Survey',header=1)

    data = data[data['Unnamed: 72'] != 'not filled']
    #st.write(data.head(5))

    #st.write(data.columns)

    cols = ['Name','Data Management Plan','Laboratory Information Management System (LIMS) [yes/no]',
            'Electronic Laboratory Notebook (ELN) [yes/no]','Lab data steward [yes/no]','Are your data FAIR? [yes/no]','BioAssay Ontology (BAO) [yes/no]','Relevant [yes/no].3',
            'Relevant [yes/no].4','Relevant [yes/no].5','Relevant [yes/no].6','Relevant [yes/no].7','Relevant [yes/no].8','Scripting [yes/no]','Version control system [yes/no]']
    
    data = data[cols]

    data = data.reset_index(drop=True)

    data = data.fillna("no")

    #st.write(data.head(2))

    def yes_no_replace(val):
        if isinstance(val, str):
            # If val starts with 'yes' (case-insensitive), replace with 'yes'
            if re.match(r'^yes\b', val, flags=re.IGNORECASE):
                return 'yes'
            else:
                # String that does NOT start with yes
                return 'no'
        else:
            # For non-string cells, also replace with 'no'
            return 'no'

        # Select columns excluding 'Name'
    cols_to_replace = [col for col in data.columns if col != 'Name']

    # Apply only on selected columns
    data[cols_to_replace] = data[cols_to_replace].applymap(yes_no_replace)


    cols = ['Name','DMP','LIMS',
            'ELN','Data steward','FAIR','BAO','Protein',
            'Nucleic acid','Organism','Analytical data','Imaging data','OMICs data','Scripting','Version control']
    
    data.columns = cols

    #st.write(data)

    # Prepare feature data (exclude 'Name' column)
    feature_data = data.drop('Name', axis=1)

    # Map 'yes' to 1, 'no' to 0 for heatmap
    z = feature_data.replace({'yes': 1, 'no': 0}).values.T 

    # Define colors: 0 => red, 1 => green
    colorscale = [[0, 'blue'], [1, 'cyan']]

    fig = go.Figure(
        go.Heatmap(
            z=z, 
            x=data['Name'], 
            y=feature_data.columns,
            colorscale=colorscale,
            colorbar=dict(tickvals=[0, 1], ticktext=['No', 'Yes']),
            showscale=True,
        )
    )

    fig.update_layout(
        title="Heat map of IMPULSE survey based on availability of data and resources",
        xaxis_title="Institutes",
        yaxis_title="Survey answers",
        xaxis_tickangle=45,
        height = 700
    )

    st.plotly_chart(fig, use_container_width=True)


