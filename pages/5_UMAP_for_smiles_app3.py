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
st.markdown("*This Streamlit app* allows user to upload **any** :blue-background[CSV file] containing a column named ***SMILES*** and a label column named ***Type***. It generates a UMAP plot (https://umap-learn.readthedocs.io/) which shows the chemical space calculated on a selection of most common RDKIT 2D descriptors. Plot colors are dependent on ***Type*** column indicating the groups.")

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

# File uploader
uploaded_file = st.file_uploader("Choose a CSV file", type="csv")

if uploaded_file is not None:
    # Read the CSV file
    df = pd.read_csv(uploaded_file)
    
    # Check if required columns exist
    if 'SMILES' not in df.columns or 'Type' not in df.columns:
        st.error("The CSV file must contain 'SMILES' and 'Type' columns.")
    else:
        # Continue with the analysis
        st.success("File uploaded successfully. Generating UMAP plot...")

        # Convert SMILES to RDKit molecules
        mols = [Chem.MolFromSmiles(smiles) for smiles in df['SMILES']]

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
        descriptors = [calc(mol) for mol in mols]
        descriptors_df = pd.DataFrame(descriptors, columns=desc_names)

        # Perform UMAP
        umap_model = UMAP(n_components=2, random_state=42, n_neighbors = 15, metric = 'cosine')
        umap_results = umap_model.fit_transform(descriptors_df)

        # Create a Bokeh ColumnDataSource
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

        # Create tooltip HTML
        tooltip_html = """
        <div>
            <div>
                <img src="@molecule_img" height="200" alt="@molecule_img" width="200">
            </div>
            <div>
                <span style="font-size: 12px; color: #666;">Type: @type</span>
            </div>
        </div>
        """

        hover = HoverTool(renderers=[scatter], tooltips=tooltip_html)
        p.add_tools(hover)

        # Show the plot in Streamlit
        st.bokeh_chart(p, use_container_width=True)

        # Function to get SVG download link
        import tempfile
        import os
        from bokeh.io import export_svgs

        def get_svg_download_link(fig):
            with tempfile.NamedTemporaryFile(suffix=".svg", delete=False) as temp_file:
                fig.output_backend = "svg"
                export_svgs(fig, filename=temp_file.name)
                
            with open(temp_file.name, "rb") as file:
                svg_content = file.read()
            
            os.unlink(temp_file.name)  # Delete the temporary file
            
            b64 = base64.b64encode(svg_content).decode()
            return f'data:image/svg+xml;base64,{b64}'

        # Add download button
        if st.button("Download Plot as SVG"):
            svg_href = get_svg_download_link(p)
            st.markdown(f'<a href="{svg_href}" download="umap_plot.svg">Download Plot as SVG</a>', unsafe_allow_html=True)