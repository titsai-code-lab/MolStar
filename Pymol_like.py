import streamlit as st
import py3Dmol
from stmol import showmol
import json
import os

st.set_page_config(layout="wide")
st.title("Boltz-2 scFv vs Human/Cyno Antigen Comparison")

def create_local_view(structure_data, file_format="cif", surface_opacity=0.4):
    view = py3Dmol.view(width=500, height=500)
    # Load the local string data into py3Dmol (supports 'cif' or 'pdb')
    view.addModel(structure_data, file_format)
    
    # Set PyMOL-like styling and transparency
    view.setStyle({'cartoon': {'color': 'spectrum'}})
    view.setBackgroundColor('white', 0)
    view.addSurface(py3Dmol.VDW, {'opacity': surface_opacity, 'color': 'white'})
    view.zoomTo()
    return view

def extract_iptm(directory, base_name):
    # Scan for the Boltz-2 confidence JSON file associated with this model
    for f in os.listdir(directory):
        if base_name in f and f.endswith('.json') and 'confidence' in f:
            with open(os.path.join(directory, f), 'r') as json_file:
                data = json.load(json_file)
                return data.get('iptm', None)
    return None

# User inputs the directory path where Boltz-2 results are stored
target_dir = st.text_input("Enter the path to your Boltz-2 output directory:", value="./")

if os.path.isdir(target_dir):
    # Find all CIF/PDB files in the folder and separate them by prefix
    files = [f for f in os.listdir(target_dir) if f.endswith(('.cif', '.pdb'))]
    hu_files = [f for f in files if f.startswith("hu_")]
    cy_files = [f for f in files if f.startswith("cy_")]

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Human Antigen Complex")
        selected_hu = st.selectbox("Select 'hu_' Model", hu_files) if hu_files else None
        
        if selected_hu:
            # Extract and display ipTM
            base_name = selected_hu.replace('.cif', '').replace('.pdb', '')
            iptm_hu = extract_iptm(target_dir, base_name)
            
            if iptm_hu is not None:
                st.metric(label="ipTM Score (Human)", value=f"{iptm_hu:.4f}")
            else:
                st.warning("ipTM score not found in JSON.")
                
            opacity_1 = st.slider("Surface Opacity", 0.0, 1.0, 0.4, key="op1")
            
            # Read structure string and render
            file_format = "pdb" if selected_hu.endswith(".pdb") else "cif"
            with open(os.path.join(target_dir, selected_hu), 'r') as f:
                view1 = create_local_view(f.read(), file_format, opacity_1)
                showmol(view1, height=500, width=500)
        elif not hu_files:
            st.info("No files starting with 'hu_' found in this directory.")

    with col2:
        st.subheader("Cyno Antigen Complex")
        selected_cy = st.selectbox("Select 'cy_' Model", cy_files) if cy_files else None
        
        if selected_cy:
            # Extract and display ipTM
            base_name = selected_cy.replace('.cif', '').replace('.pdb', '')
            iptm_cy = extract_iptm(target_dir, base_name)
            
            if iptm_cy is not None:
                st.metric(label="ipTM Score (Cyno)", value=f"{iptm_cy:.4f}")
            else:
                st.warning("ipTM score not found in JSON.")
                
            opacity_2 = st.slider("Surface Opacity", 0.0, 1.0, 0.4, key="op2")
            
            # Read structure string and render
            file_format = "pdb" if selected_cy.endswith(".pdb") else "cif"
            with open(os.path.join(target_dir, selected_cy), 'r') as f:
                view2 = create_local_view(f.read(), file_format, opacity_2)
                showmol(view2, height=500, width=500)
        elif not cy_files:
            st.info("No files starting with 'cy_' found in this directory.")
else:
    st.error("Directory not found. Please provide a valid path containing your Boltz-2 files.")
