import streamlit as st
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from umap import UMAP
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool, CustomJS
from bokeh.palettes import Category10
from bokeh.transform import factor_cmap
import io
import base64

# Streamlit app title
st.title("UMAP Plot of Molecules")

st.markdown("Uniform Manifold Approximation and Projection (UMAP) is a dimension reduction technique that can be used for visualisation similarly to t-SNE, but also for general "
            "non-linear dimension reduction. It visualises the closest topological structure of the data.")
st.markdown("*This Streamlit app* allows user to upload **any** :blue-background[CSV file] containing a column named ***SMILES*** and a label column named ***Type***. "
            "It generates a UMAP plot (https://umap-learn.readthedocs.io/) which shows the chemical space calculated on a selection of most common RDKIT 2D descriptors."
            " Plot colors are dependent on ***Type*** column indicating the groups. A user can upload another CSV file with only a column called SMILES and the plot will "
            "show with label ***USER*** the relative molecule within the plot with a :red[slightly bigger red dot].")

def plot_molecule(smiles):
    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is not None:
            img = Draw.MolToImage(mol, size=(200, 200))
            return img
        else:
            return None
    except:
        return None

def get_molecule_image_src(smiles):
    img = plot_molecule(smiles)
    if img:
        buffered = io.BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode()
        return f"data:image/png;base64,{img_str}"
    return ""
    
def safe_calc(mol):
    if mol is None:
        return [None] * len(desc_functions)
    try:
        return [func(mol) for _, func in desc_functions]
    except:
        return [None] * len(desc_functions)


# File uploaders
uploaded_file = st.file_uploader("Choose a reference CSV file with labels (Column name = Types)", type="csv")
user_file = st.file_uploader("Upload your compounds (CSV with SMILES column)", type="csv")

if uploaded_file is not None:
    # Read the CSV file
    df = pd.read_csv(uploaded_file)
    
    # Check if required columns exist
    if 'SMILES' not in df.columns or 'Type' not in df.columns:
        st.error("The CSV file must contain 'SMILES' and 'Type' columns.")
    else:
        # Process user file if uploaded
        if user_file is not None:
            # Read the user's CSV file
            user_df = pd.read_csv(user_file)
        
        # Check if SMILES column exists
        if 'SMILES' not in user_df.columns:
            st.error("The user's CSV file must contain a 'SMILES' column.")
        else:
            # Add 'Type' column with 'user' label
            user_df['Type'] = 'User'
            
            # Concatenate user_df with the main df
            df = pd.concat([df, user_df], ignore_index=True)
            df = df.dropna(subset = 'SMILES' )
            
            st.success("User compounds added successfully.")

        # Continue with the analysis
        st.success("File uploaded successfully. Generating UMAP plot...")

        # Filter out invalid molecules
        # Convert SMILES to RDKit molecules
        mols = [Chem.MolFromSmiles(smiles) for smiles in df['SMILES']]
        mols = [mol for mol in mols if mol is not None]
        
        

        # Convert SMILES column to string if it exists
        if 'SMILES' in df.columns:
            df['SMILES'] = df['SMILES'].astype(str)
            df['molecule_img'] = df['SMILES'].apply(get_molecule_image_src)

        # Calculate 2D RDKit descriptors
        desc_names = ['MolWt', 'FractionCSP3', 'HeavyAtomCount', 'NumAliphaticCarbocycles',
                      'NumAliphaticHeterocycles', 'NumAliphaticRings', 'NumAromaticCarbocycles',
                      'NumAromaticHeterocycles', 'NumAromaticRings', 'NumHAcceptors',
                      'NumHDonors', 'NumHeteroatoms', 'NumRotatableBonds', 'NumSaturatedCarbocycles',
                      'NumSaturatedHeterocycles', 'NumSaturatedRings', 'RingCount',
                      'MolLogP', 'MolMR']

        desc_functions = [(name, getattr(Descriptors, name)) for name in desc_names]
        calc = lambda m: [func(m) for _, func in desc_functions]
        descriptors = [safe_calc(mol) for mol in mols]
        descriptors_df = pd.DataFrame(descriptors, columns=desc_names)
        descriptors_df = descriptors_df.dropna()  # Remove rows with None values

        # Perform UMAP
        umap_model = UMAP(n_components=2, random_state=42, n_neighbors=30, metric='cosine', n_epochs=300, min_dist=0.3)
        umap_results = umap_model.fit_transform(descriptors_df)
        # Ensure df and umap_results are aligned
        df = df.loc[descriptors_df.index].reset_index(drop=True)

        source = ColumnDataSource(data=dict(
            x=umap_results[:, 0],
            y=umap_results[:, 1],
            smiles=df['SMILES'],
            type=df['Type'],
            molecule_img=df['molecule_img']
        ))


        # Create Bokeh figure
        p = figure(width=800, height=600, title="UMAP Plot of 2D RDKit Descriptors")

        # Create a color mapper
        colors = Category10[10][:len(df['Type'].unique())]
        color_mapper = factor_cmap('type', palette=colors, factors=df['Type'].unique())

        # Create the scatter plot
        scatter = p.scatter('x', 'y', source=source, color=color_mapper, legend_field='type', size=10, alpha=0.8)

        # Highlight user compounds
        if 'User' in df['Type'].unique():
            user_source = ColumnDataSource(data=dict(
                x=umap_results[df['Type'] == 'User', 0],
                y=umap_results[df['Type'] == 'User', 1],
                smiles=df[df['Type'] == 'User']['SMILES'],
                type=df[df['Type'] == 'User']['Type'],
                molecule_img=df[df['Type'] == 'User']['molecule_img']
            ))
            p.scatter('x', 'y', source=user_source, color='red', size=15, alpha=0.8, legend_label='User Compounds')
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