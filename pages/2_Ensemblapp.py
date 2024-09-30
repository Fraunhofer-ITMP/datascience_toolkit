import streamlit as st
import requests
import pandas as pd
import time
from PIL import Image
import re


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


def main():
    # Add the logo
    logo = Image.open("fraunhofer_ITMP-logo_900p.jpg")
    st.image(logo, width=200)  # Adjust the width as needed

    st.markdown(
        "<h1 style='text-align: center; color: green;'>Ensembl Gene Information Retriever</h1>",
        unsafe_allow_html=True,
    )

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


if __name__ == "__main__":
    main()
