"""App for CDBReg Number Generator"""

import io
import pickle
import pandas as pd
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
import datetime
import streamlit as st

st.set_page_config(
    layout="wide",
    page_title="CBDREG Number Generator",
    page_icon="🔢",
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
    "<h1 style='text-align: center; color: #149372;'> CBDREGNum Generator</h1> <br>",
    unsafe_allow_html=True,
)


def create_inchikey_dict(node_dict):
    items = [
        "PROTAC",
        "E3 Binder",
        "Linker",
        "Degrader",
        "Shuttle Protein",
        "Molecular Glue",
        "Warhead",
    ]
    temp_list = []
    for item in items:
        for obj in node_dict[item]:
            temp = [
                obj,
                node_dict[item][obj]["SMILES"],
                node_dict[item][obj]["InChI Key"],
                item,
            ]
            temp_list.append(temp)

    regnoDF = pd.DataFrame(
        temp_list, columns=["CBDREGNO", "SMILES", "InChI Key", "Type"]
    )

    regnoList = list(regnoDF["CBDREGNO"])
    ikey_org = list(regnoDF["InChI Key"])
    dict_reg_ikey = dict(zip(ikey_org, regnoList))
    return dict_reg_ikey


def removeSalt_getIkey(mol):
    remover = SaltRemover()
    stripped = remover.StripMol(mol)
    smiles = Chem.MolToSmiles(stripped)

    if "." in smiles:
        temp = smiles.split(".")
        smiles = max(temp, key=len)
        smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
        ikey = Chem.inchi.MolToInchiKey(Chem.MolFromSmiles(smiles))

        return (smiles, ikey)

    else:
        smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
        ikey = Chem.inchi.MolToInchiKey(Chem.MolFromSmiles(smiles))
        return (smiles, ikey)


def cbdregs(df):
    infile = open("./data/cbreg_data/cbdregDict.pkl", "rb")
    cbdregDict = pickle.load(infile)
    infile.close()

    temp_smiles, temp_ikey = [], []

    for smiles in list(df["SMILES"]):

        try:
            clean_smiles, ikey = removeSalt_getIkey(Chem.MolFromSmiles(smiles))
            temp_smiles.append(clean_smiles)
            temp_ikey.append(ikey)

        except:
            temp_smiles.append("Error parsing SMILES")
            temp_ikey.append("Invalid SMILES")

    df["Canonical_SMILES"] = temp_smiles
    df["InChI Key"] = temp_ikey

    temp_reg = []

    for ikey in list(df["InChI Key"]):

        if ikey != "Invalid SMILES":

            if ikey in cbdregDict:
                # st.write(ikey)
                temp_reg.append(cbdregDict[ikey])
                # st.write(cbdregDict[ikey])

            else:
                max_key = max(cbdregDict, key=cbdregDict.get)
                newRegNo = cbdregDict[max_key] + 1

                # add reg to list
                temp_reg.append(newRegNo)

                # update_dict
                cbdregDict.update({ikey: newRegNo})

        else:
            temp_reg.append(None)

    df["CBDREGNO"] = temp_reg

    # get current date and time
    now = datetime.datetime.now()
    now = now.strftime("%Y-%m-%d_%H-%M")

    filename = f"./data/cdbdregDict_{now}.pkl"

    outfile = open(filename, "wb")
    pickle.dump(cbdregDict, outfile)
    outfile.close()
    return df


st.header("CBDREGNum Generator", anchor="cbdregnum-generator", divider="gray")
st.write(
    """Welcome to the CBDREGNum Generator. This is a pre-requisite for submitting data in PROXIDRUGSDB to avoid duplicates and redundancy. This tool parses chemical SMILES representations in 2 phases. Firstly, it will remove salts from SMILES and generate canonical SMILES and InChI Key. Secondly, it checks the presence of InChI Key in our compound libary. If present, the existing CBDREGNum will be assigned to it. Else, a new CBDREGNum will be assigned. If you have your chemicals in other formats such as SDF or InChI and need help converting to SMILES, please feel free to reach us out."""
)


st.info(
    "Note: We assure you that we do not save any of your SMILES representations. We only store the InChI Key with which it is impossible to "
    "convert to the original structure.",
    icon="ℹ️",
)

st.header("Upload your data", anchor="upload-data", divider="gray")
# File upload
uploaded_file = st.file_uploader("Choose a CSV or Excel file", type=["csv", "xlsx"])
st.info('Please make sure your column is named "SMILES".', icon="ℹ️")


@st.fragment  # Prevents the app from running the code below all the time
def load_dataset():
    # File type and parameters
    file_type = uploaded_file.name.split(".")[-1]

    if file_type == "csv":
        sep = st.text_input("Enter CSV separator", ",")
        df = pd.read_csv(uploaded_file, sep=sep)
    else:  # Excel
        sheet_name = st.text_input("Enter sheet name (leave blank for first sheet)", "")
        if not sheet_name:
            sheet_name = 0

        df = pd.read_excel(uploaded_file, sheet_name=sheet_name, engine="openpyxl")

    # # Convert SMILES column to string if it exists
    # if "SMILES" in df.columns:
    #     df["SMILES"] = df["SMILES"].astype(str)
    #     df["molecule_img"] = df["SMILES"].apply(get_molecule_image_src)

    return df


# Load dataset
@st.fragment
def load_preview():
    # df = load_dataset()

    st.header("Dataset loading and preview", anchor="dataset-loading", divider="gray")
    try:
        df = load_dataset()
        st.write("Data Preview:")
        st.dataframe(df.head())

        parse_df = cbdregs(df)
        st.markdown("### CBDREGNums assigned. Please check last column to see it.")
        st.write(parse_df)

        parse_df_csv = parse_df.to_csv(index=False).encode("utf-8")

        st.download_button(
            label="Download data as CSV",
            data=parse_df_csv,
            file_name="cdbreg_df.csv",
            mime="text/csv",
        )
        # buffer to use for excel writer
        buffer = io.BytesIO()
        # download button 2 to download dataframe as xlsx
        with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
            # Write each dataframe to a different worksheet.
            parse_df.to_excel(writer, sheet_name="Sheet1", index=False)

            download2 = st.download_button(
                label="Download data as Excel",
                data=buffer,
                file_name="cdbreg_df.xlsx",
                mime="application/vnd.ms-excel",
            )

    except Exception as e:
        st.error(f"An error occurred: {str(e)}")
        st.error("Please check your file format and try again.")


if uploaded_file is not None:
    load_preview()


# footer with text and green background
current_year = datetime.datetime.today().year
st.markdown(
    f"<footer style='background-color: #149372; padding: 10px; border-radius: 10px;'>"
    f"<p style='color: white; text-align: center;'>Fraunhofer ITMP © {current_year}</p>"
    "<p style='color: white; text-align: center;'>This work has been conducted across several key projects in which ITMP has been actively involved.</p>"
    "</footer>",
    unsafe_allow_html=True,
)
