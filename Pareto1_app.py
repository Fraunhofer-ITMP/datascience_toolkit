import streamlit as st
import pandas as pd
import numpy as np
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool
from rdkit import Chem
from rdkit.Chem import Draw
import io
import base64
from PIL import Image
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows



def is_numeric(col):
    return pd.api.types.is_numeric_dtype(col)

def pareto_efficient(costs):
    is_efficient = np.ones(costs.shape[0], dtype = bool)
    for i, c in enumerate(costs):
        if is_efficient[i]:
            is_efficient[is_efficient] = np.any(costs[is_efficient]<c, axis=1)  # Keep any point with a lower cost
            is_efficient[i] = True  # And keep self
    return is_efficient

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

def main():

    logo = Image.open("fraunhofer_ITMP-logo_900p.jpg")
    st.image(logo, width=200)  # Adjust the width as needed

    st.markdown("<h1 style='text-align: center; color: green;'>Pareto Analysis App</h1>", unsafe_allow_html=True)
    #st.title("Pareto Analysis App")

    # File upload
    uploaded_file = st.file_uploader("Choose a CSV or Excel file", type=["csv", "xlsx"])
    
    if uploaded_file is not None:
        # File type and parameters
        file_type = uploaded_file.name.split('.')[-1]
        
        try:
            if file_type == 'csv':
                sep = st.text_input("Enter CSV separator", ",")
                df = pd.read_csv(uploaded_file, sep=sep)
            else:  # Excel
                sheet_name = st.text_input("Enter sheet name (leave blank for first sheet)", "")
                if sheet_name:
                    df = pd.read_excel(uploaded_file, sheet_name=sheet_name, engine='openpyxl')
                else:
                    df = pd.read_excel(uploaded_file, engine='openpyxl')

            # Convert SMILES column to string if it exists
            if 'SMILES' in df.columns:
                df['SMILES'] = df['SMILES'].astype(str)
                df['molecule_img'] = df['SMILES'].apply(get_molecule_image_src)

            st.write("Data Preview:")
            st.dataframe(df.head())

            # Select numerical columns
            numeric_cols = [col for col in df.columns if is_numeric(df[col])]
            selected_cols = st.multiselect("Select numerical columns for Pareto analysis", numeric_cols)

            if len(selected_cols) >= 2:
                # Optimization direction
                directions = {}
                for col in selected_cols:
                    direction = st.radio(f"Optimize {col}", ["Minimize", "Maximize"])
                    directions[col] = direction

                # Perform Pareto analysis
                if st.button("Perform Pareto Analysis"):
                    costs = df[selected_cols].values
                    for i, col in enumerate(selected_cols):
                        if directions[col] == "Maximize":
                            costs[:, i] = -costs[:, i]  # Invert for maximization
                    
                    is_efficient = pareto_efficient(costs)
                    
                    # Add Pareto efficiency column to the DataFrame
                    df['Pareto_Efficient'] = is_efficient
                    
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
                    p.circle(x=selected_cols[0], y=selected_cols[1], size=10, color='gray', alpha=0.5, source=source_inefficient)
                    
                    # Plot efficient points
                    p.circle(x=selected_cols[0], y=selected_cols[1], size=15, color='red', alpha=0.8, source=source_efficient)

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

                    # Create Excel file
                    def create_excel():
                        wb = Workbook()
                        ws = wb.active
                        ws.title = "Pareto Analysis Results"
                        
                        for r in dataframe_to_rows(df, index=False, header=True):
                            ws.append(r)
                        
                        excel_file = io.BytesIO()
                        wb.save(excel_file)
                        excel_file.seek(0)
                        return excel_file

                    # Download button
                    excel_file = create_excel()
                    st.download_button(
                        label="Download Excel file",
                        data=excel_file,
                        file_name="pareto_analysis_results.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                    )

        except Exception as e:
            st.error(f"An error occurred: {str(e)}")
            st.error("Please check your file format and try again.")

if __name__ == "__main__":
    main()
    