# ============================================================
# IMPORTS
# ============================================================

import pandas as pd
import numpy as np

from tqdm import tqdm

from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem import MACCSkeys
from rdkit.Chem import inchi
from rdkit.Chem.SaltRemover import SaltRemover

from mhfp.encoder import MHFPEncoder


# ============================================================
# 1. REMOVE SALTS
# ============================================================

def removeSalt(mol):

    remover = SaltRemover()

    stripped = remover.StripMol(mol)

    smiles = Chem.MolToSmiles(stripped)

    # If several fragments remain,
    # keep the largest one
    if "." in smiles:

        fragments = smiles.split(".")

        smiles = max(fragments, key=len)

    return smiles


# ============================================================
# 2. CLEAN SMILES + CREATE InChI + InChIKey
# ============================================================

def molToInChIkey(
    df,
    smiles_col,
    clean_smiles_col,
    inchi_col,
    inchikey_col
):

    smiles_list = []
    inchi_list = []
    inchikey_list = []

    for item in tqdm(df[smiles_col]):

        try:

            # SMILES → RDKit molecule
            mol = Chem.MolFromSmiles(item)

            # Remove salts
            clean_smiles = removeSalt(mol)

            # Clean SMILES → molecule again
            mol = Chem.MolFromSmiles(clean_smiles)

            # Create InChI
            inchi_value = inchi.MolToInchi(mol)

            # Create InChIKey
            inchikey_value = inchi.MolToInchiKey(mol)

            smiles_list.append(clean_smiles)
            inchi_list.append(inchi_value)
            inchikey_list.append(inchikey_value)

        except:

            smiles_list.append(None)
            inchi_list.append(None)
            inchikey_list.append(None)

    # Add new columns
    df[clean_smiles_col] = smiles_list
    df[inchi_col] = inchi_list
    df[inchikey_col] = inchikey_list

    # Remove molecules that could not be processed
    df = df.dropna(subset=[inchikey_col])

    # Reset row numbers
    df = df.reset_index(drop=True)

    return df

# ============================================================
# 3. CALCULATE RDKit DESCRIPTORS
# ============================================================

def calculate_descriptors(
    df,
    smiles_col,
    keep_cols=None
):

    descriptor_list = []

    for smi in tqdm(df[smiles_col]):

        mol = Chem.MolFromSmiles(smi)

        if mol is not None:

            descriptors = Descriptors.CalcMolDescriptors(mol)

        else:

            descriptors = {}

        descriptor_list.append(descriptors)

    # Convert descriptors into dataframe
    descriptor_df = pd.DataFrame(descriptor_list)

    # Keep important original columns
    if keep_cols is not None:

        descriptor_df = pd.concat(
            [
                df[keep_cols].reset_index(drop=True),
                descriptor_df
            ],
            axis=1
        )

    return descriptor_df


# ============================================================
# 4. FINGERPRINT GENERATORS
# ============================================================

# ECFP4
# radius 2 = diameter 4
mfpgen = rdFingerprintGenerator.GetMorganGenerator(
    radius=2,
    fpSize=1024
)

# RDKit fingerprint
rdkgen = rdFingerprintGenerator.GetRDKitFPGenerator(
    fpSize=1024
)

# MHFP
mhfp_encoder = MHFPEncoder()


# ============================================================
# 5. DIFFERENT FINGERPRINT METHODS
# ============================================================

fingerprint_methods = {

    "ECFP": lambda mol, smi:
        mfpgen.GetFingerprint(mol),

    "RDKit": lambda mol, smi:
        rdkgen.GetFingerprint(mol),

    "MACCS": lambda mol, smi:
        MACCSkeys.GenMACCSKeys(mol), # type: ignore

    "MHFP": lambda mol, smi:
        mhfp_encoder.encode(smi)
}


# ============================================================
# 6. GENERATE FINGERPRINTS
# ============================================================

def generate_fingerprints(
    df,
    smiles_col,
    fingerprint_methods):

    # Copy dataframe
    fingerprint_df = df.copy()

    # Create empty list for every fingerprint type
    fingerprint_lists = {
        name: []
        for name in fingerprint_methods
    }

    # Go through every molecule
    for smi in tqdm(fingerprint_df[smiles_col]):

        # Canonical SMILES
        try:

            can_smiles = Chem.CanonSmiles(smi)

        except:

            can_smiles = smi

        # SMILES → molecule
        mol = Chem.MolFromSmiles(can_smiles)

        # If molecule is invalid
        if mol is None:

            for name in fingerprint_methods:

                fingerprint_lists[name].append(None)

            continue

        # Generate every fingerprint type
        for name, method in fingerprint_methods.items():

            try:

                fp = method(
                    mol,
                    can_smiles
                )

            except:

                fp = None

            fingerprint_lists[name].append(fp)

    # Add fingerprints to dataframe
    for name in fingerprint_methods:

        fingerprint_df[name] = fingerprint_lists[name]

    return fingerprint_df


# ============================================================
# 7. FORMAT FINGERPRINTS INTO FEATURE COLUMNS
# ============================================================

def fingerprint_Formatting(
    df,
    fingerprint_col,
    keep_cols=None):

    # Remove rows where fingerprint failed
    fp_df = df.dropna(
        subset=[fingerprint_col])

    data = []

    for fp_object in tqdm(fp_df[fingerprint_col]):

        # MHFP is already a NumPy array
        if isinstance(fp_object, np.ndarray):

            fp = np.array(fp_object)

        # RDKit fingerprints
        else:

            fp = np.zeros(
                (0,),
                dtype=int)

            DataStructs.ConvertToNumpyArray(
                fp_object,
                fp)

        # One molecule = one row
        t = pd.DataFrame(fp).T

        # bit0, bit1, bit2...
        t.rename(
            columns=lambda x: f"bit{x}",
            inplace=True)

        data.append(t)

    # Combine all molecules
    fingerprint_dataframe = pd.concat(
        data,
        ignore_index=True)

    # Add SMILES / Type / metadata
    if keep_cols is not None:

        fingerprint_dataframe = pd.concat(
            [
                fp_df[keep_cols].reset_index(drop=True),
                fingerprint_dataframe
            ],
            axis=1)

    return fingerprint_dataframe


# ============================================================
# 8. PREPARE DATA FOR MACHINE LEARNING
# ============================================================

def prepare_ML_data(
    df,
    target_col,
    exclude_cols=None,
    feature_prefix=None):

    if exclude_cols is None:

        exclude_cols = []

    # y = what we want to predict
    y = df[target_col]

    # Fingerprint case:
    # choose bit0, bit1, bit2...
    if feature_prefix is not None:

        X = df.filter(
            regex=f"^{feature_prefix}")

    # Descriptor case:
    # remove SMILES, Type, IDs...
    else:

        X = df.drop(
            columns=[target_col] + exclude_cols)

    return X, y