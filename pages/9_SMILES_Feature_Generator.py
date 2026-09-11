# 1. Imports

import streamlit as st
import pandas as pd
import datetime
from smiles_utils import (
    molToInChIkey,
    calculate_descriptors,
    generate_fingerprints,
    fingerprint_Formatting,
    fingerprint_methods,
)


# 2. Page configuration and styling
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



# INTRODUCTION
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

# FILE UPLOAD
st.header("Upload your data", anchor="upload-data", divider="gray")
# File upload
uploaded_file = st.file_uploader("Choose a CSV, TSV or Excel file", type=["csv", "tsv", "xlsx"],)
st.info("After uploading your file, select the column containing " "the SMILES representations.", icon="ℹ️",)
# ============================================================
# LOAD THE DATASET
# ============================================================

# Run this section independently during Streamlit interactions
@st.fragment
def load_dataset():
    """Load the uploaded CSV, TSV, or Excel dataset."""

    file_type = uploaded_file.name.split(".")[-1].lower()

    # Return the file-reading position to the beginning
    uploaded_file.seek(0)

    if file_type in ["csv", "tsv"]:

        # Automatically detect comma, semicolon, tab, or another separator
        df = pd.read_csv(
            uploaded_file,
            sep=None,
            engine="python",
        )

    else:  # Excel
        sheet_name = st.text_input(
            "Enter sheet name (leave blank for first sheet)",
            "",
        )

        if not sheet_name:
            sheet_name = 0

        uploaded_file.seek(0)

        df = pd.read_excel(
            uploaded_file,
            sheet_name=sheet_name,
            engine="openpyxl",
        )

    return df

# PREVIEW AND FEATURE GENERATION
@st.fragment
def load_preview():
    """Preview the dataset and display feature options."""

    st.header("Dataset loading and preview", anchor="dataset-loading", divider="gray",)

    try:

        df = load_dataset()

        # Show file information
        file_size_kb = uploaded_file.size / 1024

        st.write(f"File size: {file_size_kb:.2f} KB")
        st.write(f"Number of rows: {len(df)}")
        st.write(f"Number of columns: {len(df.columns)}")

        # Show five rows
        st.write("Data Preview:")

        st.dataframe(df.head(5), use_container_width=True,)

        # Create the column dropdown
        smiles_column = st.selectbox("Select the column containing the SMILES", options=list(df.columns),)

        # Convert the selected column to text
        df[smiles_column] = df[smiles_column].astype(str)

        st.success(f'Column "{smiles_column}" will be used as the SMILES column.')

        # Feature-selection code continues here

        # ====================================================
        # SMILES CLEANING
        # ====================================================
        st.header("SMILES cleaning", anchor="smiles-cleaning", divider="gray",)

        if st.button("Clean SMILES", type="primary",):
            with st.spinner("Cleaning SMILES..."):

                # Pass the column selected by the user
                cleaned_df = molToInChIkey(df=df.copy(), smiles_col=smiles_column, clean_smiles_col="SMILES_clean", inchi_col="InChI", inchikey_col="InChIKey",)

            # Check whether valid molecules were found
            if cleaned_df.empty:
                st.error("No valid SMILES were found in the selected column.")
                return
             # Save the cleaned data for the next Streamlit step
            st.session_state["cleaned_df"] = cleaned_df
            st.session_state["smiles_column"] = smiles_column

            # Calculate how many rows were removed
            removed_rows = len(df) - len(cleaned_df)

            st.success(f"Successfully processed " f"{len(cleaned_df)} molecule(s).")

            if removed_rows > 0:
                st.warning(f"{removed_rows} invalid molecule(s) " "could not be processed.")

            # ====================================================
            # CLEANED DATA PREVIEW
            # ====================================================

            st.subheader("Cleaned data preview")

            st.dataframe(cleaned_df.head(5), use_container_width=True,)

            # ====================================================
            # CLEANED DATA DOWNLOAD
            # ====================================================

            cleaned_csv = cleaned_df.to_csv(index=False).encode("utf-8")

            st.download_button(label="Download cleaned SMILES", data=cleaned_csv, file_name="cleaned_smiles.csv", mime="text/csv", key="cleaned_download",)



        # ====================================================
        # FEATURE SELECTION
        # ====================================================

        if "cleaned_df" in st.session_state:

            cleaned_df = st.session_state["cleaned_df"]
            original_smiles_column = st.session_state["smiles_column"]

            st.header("Feature selection", anchor="feature-selection", divider="gray",)

            selected_features = st.multiselect(
                "Select additional molecular features to generate",
                options=[
                    "Physicochemical descriptors",
                    "Morgan fingerprints",
                    "RDKit fingerprints",
                    "MACCS fingerprints",
                    "MHFP fingerprints",
                ],
                default=[
                    "Physicochemical descriptors",
                ],
            )

            if st.button(
                "Generate selected features",
                type="primary",
            ):

                if not selected_features:
                    st.warning("Please select at least one feature type.")
                    return

                with st.spinner(
                    "Generating selected molecular features..."
                ):


                    # PHYSICOCHEMICAL DESCRIPTORS


                    if (
                        "Physicochemical descriptors"
                        in selected_features
                    ):

                        descriptor_df = calculate_descriptors(
                            df=cleaned_df,
                            smiles_col="SMILES_clean",
                            keep_cols=[
                                original_smiles_column,
                                "SMILES_clean",
                                "InChI",
                                "InChIKey",
                            ],
                        )

                        st.subheader(
                            "Physicochemical descriptors"
                        )

                        st.dataframe(
                            descriptor_df.head(5),
                            use_container_width=True,
                        )

                        descriptor_csv = descriptor_df.to_csv(
                            index=False
                        ).encode("utf-8")

                        st.download_button(
                            label="Download descriptors",
                            data=descriptor_csv,
                            file_name="smiles_descriptors.csv",
                            mime="text/csv",
                            key="descriptor_download",
                        )

                    # ============================================
                    # FINGERPRINT METHODS
                    # ============================================

                    fingerprint_selection = {
                        "Morgan fingerprints": "ECFP",
                        "RDKit fingerprints": "RDKit",
                        "MACCS fingerprints": "MACCS",
                        "MHFP fingerprints": "MHFP",
                    }

                    # ============================================
                    # GENERATE SELECTED FINGERPRINTS
                    # ============================================

                    for (
                        display_name,
                        method_name,
                    ) in fingerprint_selection.items():

                        if display_name not in selected_features:
                            continue

                        fingerprint_df = generate_fingerprints(
                            df=cleaned_df,
                            smiles_col="SMILES_clean",
                            fingerprint_methods={
                                method_name:
                                    fingerprint_methods[method_name]
                            },
                        )

                        formatted_df = fingerprint_Formatting(
                            df=fingerprint_df,
                            fingerprint_col=method_name,
                            keep_cols=[
                                original_smiles_column,
                                "SMILES_clean",
                                "InChIKey",
                            ],
                        )

                        st.subheader(display_name)

                        st.dataframe(
                            formatted_df.head(5),
                            use_container_width=True,
                        )

                        fingerprint_csv = formatted_df.to_csv(
                            index=False
                        ).encode("utf-8")

                        st.download_button(
                            label=f"Download {display_name}",
                            data=fingerprint_csv,
                            file_name=(
                                f"{method_name.lower()}"
                                "_fingerprints.csv"
                            ),
                            mime="text/csv",
                            key=(
                                f"{method_name.lower()}"
                                "_download"
                            ),
                        )
    except Exception as e:
        st.error(
            f"An error occurred: {str(e)}"
        )
        st.error(
            "Please check your file format and try again."
        )


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