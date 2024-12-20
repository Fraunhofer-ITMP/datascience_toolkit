"""App for Chemical space analysis of compounds."""

import io
import os
import base64
import streamlit as st
import pandas as pd
import numpy as np
from umap import UMAP
from bokeh.plotting import figure
from bokeh.palettes import Category10
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.transform import factor_cmap

import tempfile
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from sklearn.manifold import TSNE
from bokeh.io import export_svgs

st.set_page_config(
    layout="wide", page_title="Chemical space exploration tool", page_icon="🔎"
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
            [data-testid="stExpander"] details:hover summary{background-color: #d0e9e2;}
        </style>
        """,
    unsafe_allow_html=True,
)  # .block-conatiner controls the padding of the page, .stTabs controls the font size of the text in the tabs


def plot_molecule(smiles):
    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is not None:
            img = Draw.MolToImage(mol, size=(200, 200))
            return img
        else:
            return None
    except Exception:
        return None


def get_molecule_image_src(smiles):
    img = plot_molecule(smiles)
    if img:
        buffered = io.BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode()
        return f"data:image/png;base64,{img_str}"
    return ""


def get_svg_download_link(fig):
    """Create a download link for an SVG file."""
    with tempfile.NamedTemporaryFile(suffix=".svg", delete=False) as temp_file:
        fig.output_backend = "svg"
        export_svgs(fig, filename=temp_file.name)

    with open(temp_file.name, "rb") as file:
        svg_content = file.read()

    os.unlink(temp_file.name)  # Delete the temporary file

    b64 = base64.b64encode(svg_content).decode()
    return f"data:image/svg+xml;base64,{b64}"


st.markdown(
    "<h1 style='text-align: center; color: #149372;'> Chemical Space Analysis </h1> <br>",
    unsafe_allow_html=True,
)


def safe_calc(mol):
    if mol is None:
        return [None] * len(desc_functions)
    try:
        return [func(mol) for _, func in desc_functions]
    except:
        return [None] * len(desc_functions)


@st.fragment  # Prevents the app from running the code below all the time
def load_dataset(file: str = None):
    # File type and parameters
    file_type = file.name.split(".")[-1]

    if file_type == "csv":
        df = pd.read_csv(file, sep=",")
    elif file_type == "xlsx":
        sheet_name = st.text_input("Enter sheet name (leave blank for first sheet)", "")
        if not sheet_name:
            sheet_name = 0

        df = pd.read_excel(file, sheet_name=sheet_name, engine="openpyxl")
    else:
        st.error("Please upload a CSV or Excel file.")
        return

    df.columns = [i.lower() for i in df.columns]

    # Convert SMILES column to string if it exists
    if "smiles" in df.columns:
        df["smiles"] = df["smiles"].astype(str)
        df["molecule_img"] = df["smiles"].apply(get_molecule_image_src)

    return df


tab_1, tab_2 = st.tabs(
    [
        "Chemical space analysis",
        "UMAP Plot",
    ]
)

with tab_1:
    st.header("🔎 2D chemical space presentation", anchor="chem-space", divider="gray")
    st.markdown(
        "*This app* allows user to upload **any** :blue-background[CSV file] containing a column named ***SMILES*** and a label column named ***Type***. It generates a [t-SNE plot](https://en.wikipedia.org/wiki/T-distributed_stochastic_neighbor_embedding) or [UMAP plot](https://umap-learn.readthedocs.io/) which shows the chemical space calculated on a selection of most common RDKIT 2D descriptors. Plot colors are dependent on ***Type*** column indicating the groups."
    )

    st.header("Data loader", anchor="data-loader", divider="gray")

    # File uploader
    uploaded_file = st.file_uploader(
        "Choose a CSV file (comma separated) or excel sheet",
        type="csv",
        key="chem-space",
    )

    if uploaded_file is None:
        st.info(
            "Please upload a CSV or Excel file to proceed! Loading an example file",
            icon="ℹ️",
        )
    else:
        chem_df = load_dataset(uploaded_file)

        st.info("Reading the CSV file...")
        # Convert SMILES to RDKit molecules
        mols = [Chem.MolFromSmiles(smiles) for smiles in chem_df["smiles"]]

        # Convert SMILES column to string if it exists
        if "smiles" in chem_df.columns:
            chem_df["smiles"] = chem_df["smiles"].astype(str)
            chem_df["molecule_img"] = chem_df["smiles"].apply(get_molecule_image_src)

        st.header("Data plotter", anchor="data-plot", divider="gray")
        st.markdown("### Select plot type and calculate 2D RDKit descriptors")

        plot_type = st.radio("Plotting type", ["t-SNE", "UMAP"], horizontal=True)

        if st.button("Calculate Descriptors and Plot"):
            # Calculate 2D RDKit descriptors
            desc_names = [
                "MolWt",
                "FractionCSP3",
                "HeavyAtomCount",
                "NumAliphaticCarbocycles",
                "NumAliphaticHeterocycles",
                "NumAliphaticRings",
                "NumAromaticCarbocycles",
                "NumAromaticHeterocycles",
                "NumAromaticRings",
                "NumHAcceptors",
                "NumHDonors",
                "NumHeteroatoms",
                "NumRotatableBonds",
                "NumSaturatedCarbocycles",
                "NumSaturatedHeterocycles",
                "NumSaturatedRings",
                "RingCount",
                "MolLogP",
                "MolMR",
            ]

            desc_functions = [(name, getattr(Descriptors, name)) for name in desc_names]
            calc = lambda m: [func(m) for _, func in desc_functions]
            descriptors = [calc(mol) for mol in mols]
            descriptors_df = pd.DataFrame(descriptors, columns=desc_names)

            if plot_type == "t-SNE":
                # Perform t-SNE
                tsne = TSNE(
                    n_components=2,
                    random_state=42,
                    init="random",
                    perplexity=30,
                    n_iter=500,
                    learning_rate="auto",
                )
                results = tsne.fit_transform(descriptors_df)
            elif plot_type == "UMAP":
                # Perform UMAP
                umap_model = UMAP(
                    n_components=2, random_state=42, n_neighbors=15, metric="cosine"
                )
                results = umap_model.fit_transform(descriptors_df)

            # Create a Bokeh ColumnDataSource
            source = ColumnDataSource(
                data=dict(
                    x=results[:, 0],
                    y=results[:, 1],
                    smiles=df["smiles"],
                    type=df["type"],
                    molecule_img=df["molecule_img"],
                )
            )

            # Create Bokeh figure
            p = figure(
                width=800, height=600, title=f"{plot_type} Plot of 2D RDKit Descriptors"
            )

            # Create a color mapper
            colors = Category10[10][: len(df["type"].unique())]
            color_mapper = factor_cmap(
                "type", palette=colors, factors=df["type"].unique()
            )

            # Create the scatter plot
            scatter = p.scatter(
                "x",
                "y",
                source=source,
                color=color_mapper,
                legend_field="type",
                size=10,
                alpha=0.8,
            )

            # Create tooltip HTML
            tooltip_html = """
            <div>
                <div>
                    <img src="@molecule_img" height="200" alt="@molecule_img" width="200">
                </div>
                <div>
                    <span style="font-size: 12px; color: #666;">Type: @type</span>
                </div>
                <div>
                    <span style="font-size: 12px; color: #666;">SMILES: @smiles</span>
                </div>
            </div>
            """

            hover = HoverTool(renderers=[scatter], tooltips=tooltip_html)
            p.add_tools(hover)

            # Show the plot in Streamlit
            st.bokeh_chart(p, use_container_width=True)

            if plot_type == "t-SNE":
                file_name = "tsne_plot.svg"
            elif plot_type == "UMAP":
                file_name = "umap_plot.svg"

            # Add download button
            if st.button("Download Plot as SVG"):
                svg_href = get_svg_download_link(p)
                st.markdown(
                    f'<a href="{svg_href}" download=f"{file_name}">Download Plot as SVG</a>',
                    unsafe_allow_html=True,
                )


with tab_2:
    st.header("🔎 UMAP plot", anchor="umap", divider="gray")
    st.markdown(
        "Uniform Manifold Approximation and Projection (UMAP) is a dimension reduction technique that can be used for visualisation similarly to t-SNE, but also for general non-linear dimension reduction. It visualises the closest topological structure of the data."
    )
    st.markdown(
        "*This Streamlit app* allows user to upload **any** :blue-background[CSV file] containing a column named ***SMILES*** and a label column named ***Type***. It generates a [UMAP plot](https://umap-learn.readthedocs.io/) which shows the chemical space calculated on a selection of most common RDKIT 2D descriptors. Plot colors are dependent on ***Type*** column indicating the groups. A user can upload another CSV file with only a column called SMILES and the plot will show with label ***USER*** the relative molecule within the plot with a :red[slightly bigger red dot]."
    )

    st.header("Data loader", anchor="data-loader", divider="gray")

    # File uploaders
    uploaded_file = st.file_uploader(
        "Choose a reference CSV file with labels (Column name = Types)", type="csv"
    )
    user_file = st.file_uploader(
        "Upload your compounds (CSV with SMILES column)", type="csv"
    )

    if uploaded_file is not None and user_file is not None:
        # Read the CSV file
        ref_df = pd.read_csv(uploaded_file)
        ref_df.columns = [i.lower() for i in ref_df.columns]

        # Read the user's CSV file
        user_df = pd.read_csv(user_file)
        user_df.columns = [i.lower() for i in user_df.columns]

        type_col = [i for i in ref_df.columns if "type" in i][0]
        smile_col = [i for i in ref_df.columns if "smile" in i][0]

        # Check if SMILES column exists
        if type_col not in ref_df.columns or smile_col not in ref_df.columns:
            st.error(
                f"The user's CSV file must contain a {smile_col} and {type_col} column."
            )

        if smile_col not in user_df.columns:
            st.error(f"The user's CSV file must contain a {smile_col} column.")

        if type_col not in user_df.columns:
            user_df["type"] = "User"

        full_df = pd.concat([ref_df, user_df], ignore_index=True)
        full_df = full_df.dropna(subset=[smile_col])

        # Continue with the analysis
        st.success("File uploaded successfully. Generating UMAP plot...")

        # Filter out invalid molecules
        # Convert SMILES to RDKit molecules
        mols = [Chem.MolFromSmiles(smiles) for smiles in full_df[smile_col]]
        mols = [mol for mol in mols if mol is not None]

        # Convert SMILES column to string if it exists
        full_df["smiles"] = full_df[smile_col].astype(str)
        full_df["molecule_img"] = full_df["smiles"].apply(get_molecule_image_src)

        # Calculate 2D RDKit descriptors
        desc_names = [
            "MolWt",
            "FractionCSP3",
            "HeavyAtomCount",
            "NumAliphaticCarbocycles",
            "NumAliphaticHeterocycles",
            "NumAliphaticRings",
            "NumAromaticCarbocycles",
            "NumAromaticHeterocycles",
            "NumAromaticRings",
            "NumHAcceptors",
            "NumHDonors",
            "NumHeteroatoms",
            "NumRotatableBonds",
            "NumSaturatedCarbocycles",
            "NumSaturatedHeterocycles",
            "NumSaturatedRings",
            "RingCount",
            "MolLogP",
            "MolMR",
        ]

        desc_functions = [(name, getattr(Descriptors, name)) for name in desc_names]
        calc = lambda m: [func(m) for _, func in desc_functions]
        descriptors = [safe_calc(mol) for mol in mols]
        descriptors_df = pd.DataFrame(descriptors, columns=desc_names)
        descriptors_df = descriptors_df.dropna()  # Remove rows with None values

        # Perform UMAP
        umap_model = UMAP(
            n_components=2,
            random_state=42,
            n_neighbors=30,
            metric="cosine",
            n_epochs=300,
            min_dist=0.3,
        )
        umap_results = umap_model.fit_transform(descriptors_df)
        # Ensure df and umap_results are aligned
        df = full_df.loc[descriptors_df.index].reset_index(drop=True)

        source = ColumnDataSource(
            data=dict(
                x=umap_results[:, 0],
                y=umap_results[:, 1],
                smiles=df["smiles"],
                type=df[type_col],
                molecule_img=df["molecule_img"],
            )
        )

        # Create Bokeh figure
        p = figure(width=800, height=600, title="UMAP Plot of 2D RDKit Descriptors")

        # Create a color mapper
        colors = Category10[10][: len(df[type_col].unique())]
        color_mapper = factor_cmap(
            "type", palette=colors, factors=df[type_col].unique()
        )

        # Create the scatter plot
        scatter = p.scatter(
            "x",
            "y",
            source=source,
            color=color_mapper,
            legend_field="type",
            size=10,
            alpha=0.8,
        )

        # Highlight user compounds
        if "User" in df[type_col].unique():
            user_source = ColumnDataSource(
                data=dict(
                    x=umap_results[df[type_col] == "User", 0],
                    y=umap_results[df[type_col] == "User", 1],
                    smiles=df[df[type_col] == "User"]["smiles"],
                    type=df[df[type_col] == "User"][type_col],
                    molecule_img=df[df[type_col] == "User"]["molecule_img"],
                )
            )
            p.scatter(
                "x",
                "y",
                source=user_source,
                color="red",
                size=15,
                alpha=0.8,
                legend_label="User Compounds",
            )
        # Create tooltip HTML
        tooltip_html = """
        <div>
            <div>
                <img src="@molecule_img" height="200" alt="@molecule_img" width="200">
            </div>
            <div>
                <span style="font-size: 12px; color: #666;">Type: @type</span>
            </div>
            <div>
                <span style="font-size: 12px; color: #666;">SMILES: @smiles</span>
            </div>
        </div>
        """

        hover = HoverTool(renderers=[scatter], tooltips=tooltip_html)
        p.add_tools(hover)

        # Show the plot in Streamlit
        st.bokeh_chart(p, use_container_width=True)

# footer with text and green background
st.markdown(
    "<footer style='background-color: #149372; padding: 10px; border-radius: 10px;'>"
    "<p style='color: white; text-align: center;'>Fraunhofer ITMP © 2024</p>"
    "<p style='color: white; text-align: center;'>This work has been conducted across several key projects in which ITMP has been actively involved.</p>"
    "</footer>",
    unsafe_allow_html=True,
)
