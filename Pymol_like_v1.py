import streamlit as st
import py3Dmol
from stmol import showmol

st.set_page_config(layout="wide")
st.title("Interactive PyMOL-like Web Viewer")

def create_molecular_view(pdb_id, surface_opacity):
    if not pdb_id:
        return None
        
    # Initialize WebGL viewer fetching data from PDB
    view = py3Dmol.view(query=f"pdb:{pdb_id}", width=500, height=500)
    
    # Set rendering style (similar to PyMOL's cartoon representation)
    view.setStyle({'cartoon': {'color': 'spectrum'}})
    
    # Background Transparency (0 sets alpha to fully transparent)
    view.setBackgroundColor('white', 0)
    
    # Surface Transparency
    view.addSurface(py3Dmol.VDW, {'opacity': surface_opacity, 'color': 'white'})
    
    view.zoomTo()
    return view

# Create side-by-side layout using Streamlit columns
col1, col2 = st.columns(2)

with col1:
    st.subheader("Viewer 1")
    # Option A: Use a selectbox with a predefined list of interesting PDBs
    pdb_options = ["1IGT", "1KB5", "1YY9", "7BWJ", "Other"]
    selected_pdb_1 = st.selectbox("Select an antibody structure:", pdb_options, key="box1")
    
    # Allow custom typing if "Other" is selected
    if selected_pdb_1 == "Other":
        selected_pdb_1 = st.text_input("Enter a 4-letter PDB ID:", value="1FVK", key="text1").upper()
    
    opacity_1 = st.slider("Surface Opacity", 0.0, 1.0, 0.4, key="slider1")
    
    if selected_pdb_1:
        view1 = create_molecular_view(selected_pdb_1, surface_opacity=opacity_1)
        showmol(view1, height=500, width=500)

with col2:
    st.subheader("Viewer 2")
    # Option B: Directly use a text input allowing the user to type any PDB ID
    selected_pdb_2 = st.text_input("Enter any 4-letter PDB ID:", value="1KB5", key="text2").upper()
    
    opacity_2 = st.slider("Surface Opacity", 0.0, 1.0, 0.6, key="slider2")
    
    if selected_pdb_2:
        view2 = create_molecular_view(selected_pdb_2, surface_opacity=opacity_2)
        showmol(view2, height=500, width=500)
