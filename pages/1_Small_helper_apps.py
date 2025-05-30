"""App for Pareto analysis of multi-objective optimization problems."""

import io
import os
import re
import requests
import base64
import time
import streamlit as st
import pandas as pd
import numpy as np
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool

import tempfile
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows
from rdkit import Chem
from rdkit.Chem import Draw
from bokeh.io import export_svgs
import datetime

st.set_page_config(layout="wide", page_title="Small helper apps", page_icon="🧪")

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


# Functions
def is_numeric(col):
    return pd.api.types.is_numeric_dtype(col)


def pareto_efficient(costs):
    is_efficient = np.ones(costs.shape[0], dtype=bool)
    for i, c in enumerate(costs):
        if is_efficient[i]:
            is_efficient[is_efficient] = np.any(
                costs[is_efficient] < c, axis=1
            )  # Keep any point with a lower cost
            is_efficient[i] = True  # And keep self
    return is_efficient


def create_excel():
    """Create Excel file to download."""
    wb = Workbook()
    ws = wb.active
    ws.title = "Pareto Analysis Results"

    for r in dataframe_to_rows(df, index=False, header=True):
        ws.append(r)

    excel_file = io.BytesIO()
    wb.save(excel_file)
    excel_file.seek(0)
    return excel_file


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


def get_ensembl_info(ensembl_list):
    base_url = "https://rest.ensembl.org"
    result_list = []

    progress_bar = st.progress(0)
    status_text = st.empty()

    for i, ensembl_id in enumerate(ensembl_list):
        try:
            # Fetch gene information
            gene_url = f"{base_url}/lookup/id/{ensembl_id}?expand=1"
            response = requests.get(
                gene_url, headers={"Content-Type": "application/json"}
            )
            response.raise_for_status()
            gene_data = response.json()

            # Fetch cross-references (including UniProt)
            xref_url = f"{base_url}/xrefs/id/{ensembl_id}?external_db=UniProt/SWISSPROT"
            xref_response = requests.get(
                xref_url, headers={"Content-Type": "application/json"}
            )
            xref_response.raise_for_status()
            xref_data = xref_response.json()

            uniprot_id = next(
                (
                    x["primary_id"]
                    for x in xref_data
                    if x["dbname"] == "UniProt/SWISSPROT"
                ),
                "Not found",
            )

            result_list.append(
                {
                    "Ensembl_ID": ensembl_id,
                    "Gene_Symbol": gene_data.get("display_name", "Not found"),
                    "Gene_Description": gene_data.get("description", "Not found"),
                    "Chromosome": gene_data.get("seq_region_name", "Not found"),
                    "Start": gene_data.get("start", "Not found"),
                    "End": gene_data.get("end", "Not found"),
                    "Strand": gene_data.get("strand", "Not found"),
                    "Biotype": gene_data.get("biotype", "Not found"),
                    "UniProt_ID": uniprot_id,
                }
            )

        except requests.exceptions.RequestException as e:
            st.error(f"Error processing {ensembl_id}: {str(e)}")
            result_list.append(
                {
                    "Ensembl_ID": ensembl_id,
                    "Gene_Symbol": "Error",
                    "Gene_Description": "Error",
                    "Chromosome": "Error",
                    "Start": "Error",
                    "End": "Error",
                    "Strand": "Error",
                    "Biotype": "Error",
                    "UniProt_ID": "Error",
                }
            )

        # Update progress
        progress = (i + 1) / len(ensembl_list)
        progress_bar.progress(progress)
        status_text.text(f"Processed {i+1}/{len(ensembl_list)} IDs")

        time.sleep(0.1)  # To avoid overwhelming the Ensembl server

    return pd.DataFrame(result_list)


def get_uniprot_info(gene_list):
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    result_list = []

    progress_bar = st.progress(0)
    status_text = st.empty()

    for i, gene in enumerate(gene_list):
        query = f"(gene:{gene}) AND (organism_id:9606)"
        params = {"query": query, "format": "json"}

        url = f"{base_url}?{requests.compat.urlencode(params)}"

        try:
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()

            if data["results"]:
                entry = data["results"][0]

                recommended_name = (
                    entry["proteinDescription"]["recommendedName"]["fullName"]["value"]
                    if "proteinDescription" in entry
                    and "recommendedName" in entry["proteinDescription"]
                    else "Not found"
                )

                function = next(
                    (
                        comment["texts"][0]["value"]
                        for comment in entry.get("comments", [])
                        if comment["commentType"] == "FUNCTION"
                    ),
                    "Not found",
                )

                subcellular_location = next(
                    (
                        loc["location"]["value"]
                        for comment in entry.get("comments", [])
                        if comment["commentType"] == "SUBCELLULAR LOCATION"
                        for loc in comment.get("subcellularLocations", [])
                    ),
                    "Not found",
                )

                structural_family = next(
                    (
                        comment["texts"][0]["value"]
                        for comment in entry.get("comments", [])
                        if comment["commentType"] == "SIMILARITY"
                    ),
                    "Not found",
                )

                result_list.append(
                    {
                        "Gene": gene,
                        "UniProt_ID": entry["primaryAccession"],
                        "Recommended_name": recommended_name,
                        "Subcellular_location": subcellular_location,
                        "Function": function,
                        "Structural_family": structural_family,
                    }
                )
            else:
                result_list.append(
                    {
                        "Gene": gene,
                        "UniProt_ID": "Not found",
                        "Recommended_name": "Not found",
                        "Subcellular_location": "Not found",
                        "Function": "Not found",
                        "Structural_family": "Not found",
                    }
                )
        except requests.exceptions.RequestException as e:
            st.error(f"Error processing {gene}: {str(e)}")
            result_list.append(
                {
                    "Gene": gene,
                    "UniProt_ID": "Error",
                    "Recommended_name": "Error",
                    "Subcellular_location": "Error",
                    "Function": "Error",
                    "Structural_family": "Error",
                }
            )

        # Update progress
        progress = (i + 1) / len(gene_list)
        progress_bar.progress(progress)
        status_text.text(f"Processed {i+1}/{len(gene_list)} genes")

        time.sleep(0.1)  # To avoid overwhelming the UniProt server

    return pd.DataFrame(result_list)


st.markdown(
    "<h1 style='text-align: center; color: #149372;'>Small helper applications</h1>",
    unsafe_allow_html=True,
)

tab1, tab2 = st.tabs(
    [
        "Pareto Analysis",
        "Ensemble grounder",
    ]
)

with tab1:
    col1, col2 = st.columns(
        [5, 2],
        vertical_alignment="center",
    )
    with col1:
        st.header(
            "📈 Pareto Analysis App",
            anchor="pareto-analysis-app",
        )
        st.markdown(
            "*This app* allows user to upload **any** :blue-background[CSV file] or **any** :blue-background[EXCEL file] containing a column named ***SMILES*** and some numerical columns whose data have to be minimized or maximize during the [Pareto front calculation](https://en.wikipedia.org/wiki/Pareto_front). The app generates a subselection of compounds which calculated using the constraints/requirements indicated as: (1) user selection of numerical columns and (2) maximization/minimization of those columns."
        )
    with col2:
        st.image("images/pages/pareto_front.png", width=250)

    st.header(
        "Pareto parameterization", anchor="pareto-parameterization", divider="gray"
    )
    # File upload
    uploaded_file = st.file_uploader("Choose a CSV or Excel file", type=["csv", "xlsx"])

    if uploaded_file is None:
        st.info(
            "Please upload a CSV or Excel file to proceed! Currently, using example file",
            icon="ℹ️",
        )

    @st.fragment  # Prevents the app from running the code below all the time
    def load_dataset():
        # File type and parameters
        file_type = uploaded_file.name.split(".")[-1]

        if file_type == "csv":
            sep = st.text_input("Enter CSV separator", ",")
            df = pd.read_csv(uploaded_file, sep=sep)
        else:  # Excel
            sheet_name = st.text_input(
                "Enter sheet name (leave blank for first sheet)", ""
            )
            if not sheet_name:
                sheet_name = 0

            df = pd.read_excel(uploaded_file, sheet_name=sheet_name, engine="openpyxl")

        # Convert SMILES column to string if it exists
        if "SMILES" in df.columns:
            df["SMILES"] = df["SMILES"].astype(str)
            df["molecule_img"] = df["SMILES"].apply(get_molecule_image_src)

        return df

    # Load dataset
    st.header("Dataset loading and preview", anchor="dataset-loading", divider="gray")
    try:
        if uploaded_file is None:
            st.info(
                "Please upload a CSV or Excel file to proceed! Currently, using example file",
                icon="ℹ️",
            )
            df = pd.read_excel("data/SigUp-v-SigDown_July2024.xlsx")
        else:
            df = load_dataset()

        st.write("Data Preview:")
        st.dataframe(df.head())

    except Exception as e:
        st.error(f"An error occurred: {str(e)}")
        st.error("Please check your file format and try again.")

    # Select numerical columns
    st.header(
        "Pareto parameterization", anchor="pareto-parameterization", divider="gray"
    )
    numeric_cols = [col for col in df.columns if is_numeric(df[col])]
    selected_cols = st.multiselect(
        "Select numerical columns for Pareto analysis",
        numeric_cols,
        max_selections=2,
    )

    # Optimization direction
    directions = {}
    for col in selected_cols:
        direction = st.radio(f"Optimize {col}", ["Minimize", "Maximize"])
        directions[col] = direction

    # Perform Pareto analysis
    st.header("Pareto analysis", anchor="pareto-analysis", divider="gray")
    if st.button("Perform Pareto Analysis"):
        costs = df[selected_cols].values
        for i, col in enumerate(selected_cols):
            if directions[col] == "Maximize":
                costs[:, i] = -costs[:, i]  # Invert for maximization

        is_efficient = pareto_efficient(costs)

        # Add Pareto efficiency column to the DataFrame
        df["Pareto_Efficient"] = is_efficient

        # Interactive plot with molecule visualization
        st.write("Hover over points to see molecule structures:")

        # Round the selected columns to 3 decimal places
        for col in selected_cols:
            df[col] = df[col].round(3)

        df_efficient = df[is_efficient]
        df_inefficient = df[~is_efficient]

        source_efficient = ColumnDataSource(df_efficient)
        source_inefficient = ColumnDataSource(df_inefficient)

        p = figure(width=800, height=600, title="Interactive Pareto Front")

        # Plot inefficient points
        p.circle(
            x=selected_cols[0],
            y=selected_cols[1],
            size=10,
            color="gray",
            alpha=0.5,
            source=source_inefficient,
        )

        # Plot efficient points
        p.circle(
            x=selected_cols[0],
            y=selected_cols[1],
            size=15,
            color="red",
            alpha=0.8,
            source=source_efficient,
        )

        # Create tooltip HTML
        tooltip_html = """
        <div>
            <div>
                <img src="@molecule_img" height="200" alt="@molecule_img" width="200">
            </div>
        """

        for col in selected_cols:
            tooltip_html += f"""
            <div>
                <span style="font-size: 12px; color: #666;">{col}: @{{{col}}}</span>
            </div>
            """

        tooltip_html += """
            <div>
                <span style="font-size: 12px; color: #666;">SMILES: @SMILES</span>
            </div>
        </div>
        """

        hover = HoverTool(tooltips=tooltip_html)

        p.add_tools(hover)
        p.xaxis.axis_label = selected_cols[0]
        p.yaxis.axis_label = selected_cols[1]

        st.bokeh_chart(p, use_container_width=True)

        # Download button
        excel_file = create_excel()
        st.download_button(
            label="Download Excel file",
            data=excel_file,
            file_name="pareto_analysis_results.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        )


with tab2:
    st.header(
        "🧬 Ensembl Gene Information Retriever",
        anchor="ensemble-app",
    )

    st.header("Grounding genes", anchor="tsne-plot", divider="gray")

    # Input for gene IDs or symbols
    input_text = st.text_area(
        "Enter Ensembl IDs or HGNC symbols (one per line):",
        "ENSG00000234186\nENSG00000184385\nIL10\nHDC\nENSG00000180251\nFAM216B\nDRD3\nENSG00000229590",
    )

    if st.button("Retrieve Information"):
        input_list = [id.strip() for id in input_text.split("\n") if id.strip()]

        if input_list:
            # Separate ENSG IDs and HGNC symbols
            ensembl_ids = [id for id in input_list if re.match(r"ENSG\d+", id)]
            hgnc_symbols = [id for id in input_list if not re.match(r"ENSG\d+", id)]

            ensembl_results = None
            uniprot_results = None

            if ensembl_ids:
                with st.spinner("Retrieving information for Ensembl IDs..."):
                    ensembl_results = get_ensembl_info(ensembl_ids)
                    st.success(
                        f"Retrieved information for {len(ensembl_ids)} Ensembl IDs"
                    )
                    st.dataframe(ensembl_results)

                    # Option to download Ensembl results as CSV
                    csv_ensembl = ensembl_results.to_csv(index=False)
                    st.download_button(
                        label="Download Ensembl data as CSV",
                        data=csv_ensembl,
                        file_name="ensembl_info.csv",
                        mime="text/csv",
                    )

            if hgnc_symbols:
                with st.spinner("Retrieving information for HGNC symbols..."):
                    uniprot_results = get_uniprot_info(hgnc_symbols)
                    st.success(
                        f"Retrieved information for {len(hgnc_symbols)} HGNC symbols"
                    )
                    st.dataframe(uniprot_results)

                    # Option to download UniProt results as CSV
                    csv_uniprot = uniprot_results.to_csv(index=False)
                    st.download_button(
                        label="Download UniProt data as CSV",
                        data=csv_uniprot,
                        file_name="uniprot_info.csv",
                        mime="text/csv",
                    )

            if not ensembl_ids and not hgnc_symbols:
                st.warning("No valid Ensembl IDs or HGNC symbols found.")
        else:
            st.warning("Please enter at least one Ensembl ID or HGNC symbol.")


# footer with text and green background
current_year = datetime.datetime.today().year
st.markdown(
    f"<footer style='background-color: #149372; padding: 10px; border-radius: 10px;'>"
    f"<p style='color: white; text-align: center;'>Fraunhofer ITMP © {current_year}</p>"
    "<p style='color: white; text-align: center;'>This work has been conducted across several key projects in which ITMP has been actively involved.</p>"
    "</footer>",
    unsafe_allow_html=True,
)
