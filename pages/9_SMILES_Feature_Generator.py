import streamlit as st

import pandas as pd
import datetime


st.set_page_config(
    layout="wide",
    page_title="SMILES Feature Generator",
    page_icon="🧪",
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

        .stTabs [data-baseweb="tab-list"] button
        [data-testid="stMarkdownContainer"] p {
            font-size: 1.3rem;
        }
    </style>
    """,
    unsafe_allow_html=True,
)

st.markdown(
    """
    <h1 style="text-align: center; color: #149372;"> SMILES Feature Generator</h1> <br>""",
    unsafe_allow_html=True,
)










st.header(
    "SMILES preprocessing and molecular feature generation",
    anchor="smiles-feature-generator",
    divider="gray",
)

st.write(
    """
    This tool cleans and standardizes molecular SMILES representations
    and generates molecular descriptors and fingerprints for machine-learning
    applications.

    Users can upload a CSV or Excel file, select the column containing the
    SMILES representations and download the generated ML-ready feature tables.
    """
)

st.info(
    "Uploaded structures are processed for feature generation and are not "
    "saved by this application.",
    icon="ℹ️",
)


st.header("Upload your data", anchor="upload-data", divider="gray")
# File upload
uploaded_file = st.file_uploader("Choose a CSV or Excel file", type=["csv", "xlsx"])
st.info('Please make sure your SMILES column is named "SMILES".', icon="ℹ️")


@st.fragment # Prevents the app from running the code below all the time
def load_dataset():
    """Load the uploaded CSV or Excel dataset."""

    file_type = uploaded_file.name.split(".")[-1]

    if file_type == "csv":
        sep = st.text_input("Enter CSV separator", ",")
        df = pd.read_csv(uploaded_file,sep=sep)

    else:  # Excel
        sheet_name = st.text_input("Enter sheet name (leave blank for first sheet)","")

        if not sheet_name:
            sheet_name = 0

        df = pd.read_excel(uploaded_file,sheet_name=sheet_name,engine="openpyxl")

    return df

#Load dataset
@st.fragment
def load_preview():
    """Preview the dataset and display feature options."""

    st.header("Dataset loading and preview", anchor="dataset-loading", divider="gray")
    try:
        df = load_dataset()

        # Check for the required SMILES column
        if "SMILES" not in df.columns:
            st.error('The uploaded dataset must contain a column named "SMILES".')
            return

        # Ensure SMILES values are strings
        df["SMILES"] = df["SMILES"].astype(str)

        st.write("Data Preview:")
        st.dataframe(df.head(), use_container_width=True)

        st.header("Feature selection", anchor="feature-selection", divider="gray")

        selected_features = st.multiselect(
            "Select the molecular features to generate",
            options=[
                "Cleaned and standardized SMILES",
                "Physicochemical descriptors",
                "Morgan fingerprints",
                "RDKit fingerprints",
                "MACCS fingerprints",
            ],
            default=[
                "Cleaned and standardized SMILES",
                "Physicochemical descriptors",
            ],
        )

        if st.button(
            "Generate molecular features",
            type="primary",
        ):
            st.info(
                "The interface is working. Your SMILES-processing "
                "functions will be connected here next.",
                icon="ℹ️",
            )

            # This will be added later:
            #
            # results = process_smiles_dataframe(
            #     df,
            #     selected_features,
            # )
            #
            # st.dataframe(results)
            #
            # results_csv = results.to_csv(
            #     index=False
            # ).encode("utf-8")
            #
            # st.download_button(
            #     label="Download molecular features",
            #     data=results_csv,
            #     file_name="smiles_features.csv",
            #     mime="text/csv",
            # )

    except Exception as e:
        st.error(f"An error occurred: {str(e)}")
        st.error("Please check your file format and try again.")


if uploaded_file is not None:
    load_preview()


# Footer can remain the same
current_year = datetime.datetime.today().year
st.markdown(
    f"<footer style='background-color: #149372; padding: 10px; border-radius: 10px;'>"
    f"<p style='color: white; text-align: center;'>Fraunhofer ITMP © {current_year}</p>"
    "<p style='color: white; text-align: center;'>This work has been conducted across several key projects in which ITMP has been actively involved.</p>"
    "</footer>",
    unsafe_allow_html=True,
)