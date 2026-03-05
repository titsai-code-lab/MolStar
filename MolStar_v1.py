import base64
import os
import re
import tempfile
import pandas as pd
from dash import Dash, html, dcc, Input, Output, State, dash_table, no_update, ctx, ALL
import dash_molstar
from dash_molstar.utils import molstar_helper

# Biopython imports for structural processing
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select, NeighborSearch
from Bio.PDB.mmcifio import MMCIFIO
from Bio.SeqUtils import seq1

app = Dash(__name__)

# --- BIOPYTHON HELPER ---
class NonSolventSelect(Select):
    """Filters out water and typical solvents during Biopython PDB save"""
    def accept_residue(self, residue):
        if residue.id[0] in ('W', 'H_HOH', 'H_WAT') or residue.resname.strip() in ('HOH', 'WAT'):
            return 0
        return 1

def analyze_structure(filepath, ext, side, sort_val='A'):
    """Parses structure, builds interaction table, and generates interactive sequence UI."""
    if ext == '.cif':
        parser = MMCIFParser(QUIET=True)
        io = MMCIFIO()
    else:
        parser = PDBParser(QUIET=True)
        io = PDBIO()

    structure = parser.get_structure('protein', filepath)

    # 1. Clean the structure
    clean_path = filepath + "_clean" + ext
    io.set_structure(structure)
    io.save(clean_path, NonSolventSelect())
    clean_structure = parser.get_structure('clean', clean_path)

    # 2. Calculate inter-chain interactions and track interacting residues
    atoms = list(clean_structure.get_atoms())
    ns = NeighborSearch(atoms)
    contact_dict = {}
    interacting_set = set() # To store (chain_id, res_num) of interacting residues

    for a1, a2 in ns.search_all(5.0):
        c1 = a1.get_parent().get_parent()
        c2 = a2.get_parent().get_parent()

        if c1.id != c2.id: # Different chains only
            r1 = a1.get_parent()
            r2 = a2.get_parent()
            
            # Save for UI highlighting
            interacting_set.add((c1.id, r1.id[1]))
            interacting_set.add((c2.id, r2.id[1]))

            # Maintain alphabetical consistency
            if c1.id > c2.id:
                c1, c2 = c2, c1
                r1, r2 = r2, r1
                a1, a2 = a2, a1

            r1_aa = seq1(r1.resname, custom_map={'UNK': 'X'})
            r1_name = f"{r1_aa}{r1.id[1]}"
            r2_aa = seq1(r2.resname, custom_map={'UNK': 'X'})
            r2_name = f"{r2_aa}{r2.id[1]}"

            key = (c1.id, r1_name, r1.id[1], c2.id, r2_name, r2.id[1])
            dist = a1 - a2

            if key not in contact_dict or dist < contact_dict[key]:
                contact_dict[key] = dist

    # 3. Build Interactive FASTA Sequence UI
    seq_components = []
    for model in clean_structure:
        for chain in model:
            # Add FASTA Header
            seq_components.append(html.Div(f">Chain {chain.id}", style={'fontWeight': 'bold', 'color': '#0074d9', 'marginTop': '10px'}))
            
            chain_spans = []
            residues = [res for res in chain if res.id[0] == ' ']
            for res in residues:
                res_num = res.id[1]
                aa = seq1(res.resname, custom_map={'UNK': 'X'})
                is_interacting = (chain.id, res_num) in interacting_set
                
                # Style logic: Red & Bold if interacting, standard if not
                style = {'cursor': 'pointer', 'display': 'inline-block'}
                if is_interacting:
                    style.update({'color': '#d9534f', 'fontWeight': 'bold', 'textDecoration': 'underline'})
                else:
                    style.update({'color': '#333'})
                    
                chain_spans.append(html.Span(
                    aa,
                    id={'type': f'seq-aa-{side}', 'chain': chain.id, 'res': res_num},
                    style=style,
                    title=f"{aa}{res_num} (Chain {chain.id})" + (" - INTERACTING" if is_interacting else "")
                ))
            
            # Wrap chain sequence in a div block
            seq_components.append(html.Div(chain_spans, style={'lineHeight': '1.5'}))
        break  # First model only

    # 4. Structure exact requested columns
    df_data = []
    for (c1, r1, r1_num, c2, r2, r2_num), d in contact_dict.items():
        df_data.append({
            'Column A (Chain A)': c1,
            'Column B (Res A)': r1,
            '_sort_res_A': r1_num,
            'Column C (Chain B)': c2,
            'Column D (Res B)': r2,
            '_sort_res_B': r2_num,
            'Column E (Distance Å)': f"{d:.1f}"
        })

    df = pd.DataFrame(df_data)
    if not df.empty:
        if sort_val == 'B': # Primary sort by Chain B residue position
            df = df.sort_values(by=[
                'Column C (Chain B)',
                '_sort_res_B',
                'Column A (Chain A)',
                '_sort_res_A',
                'Column E (Distance Å)'
            ])
        else: # Primary sort by Chain A residue position (default)
            df = df.sort_values(by=[
                'Column A (Chain A)',
                '_sort_res_A',
                'Column C (Chain B)',
                '_sort_res_B',
                'Column E (Distance Å)'
            ])
        df = df.drop(columns=['_sort_res_A', '_sort_res_B'])
    
    return clean_path, seq_components, df.to_dict('records')


def compare_tables(left_data, right_data, compare_col):
    """Compare left and right tables to find matching residues in selected column"""
    if not left_data or not right_data:
        return []
    
    # Determine which column to compare
    if compare_col == 'B':
        col_name = 'Column B (Res A)'
        chain_col = 'Column A (Chain A)'
    else:  # 'D'
        col_name = 'Column D (Res B)'
        chain_col = 'Column C (Chain B)'
    
    # Build comparison results
    comparison_results = []
    
    for left_row in left_data:
        left_residue = left_row[col_name]
        left_chain = left_row[chain_col]
        
        for right_row in right_data:
            right_residue = right_row[col_name]
            right_chain = right_row[chain_col]
            
            # Match if residue identity is the same (including chain and position)
            if left_residue == right_residue and left_chain == right_chain:
                comparison_results.append({
                    'Matched Residue': f"{left_chain}:{left_residue}",
                    'Left - Chain A': left_row['Column A (Chain A)'],
                    'Left - Res A': left_row['Column B (Res A)'],
                    'Left - Chain B': left_row['Column C (Chain B)'],
                    'Left - Res B': left_row['Column D (Res B)'],
                    'Left - Distance': left_row['Column E (Distance Å)'],
                    'Right - Chain A': right_row['Column A (Chain A)'],
                    'Right - Res A': right_row['Column B (Res A)'],
                    'Right - Chain B': right_row['Column C (Chain B)'],
                    'Right - Res B': right_row['Column D (Res B)'],
                    'Right - Distance': right_row['Column E (Distance Å)']
                })
    
    return comparison_results


# --- DASH APP LAYOUT ---
upload_style = {
    'width': '100%', 'height': '60px', 'lineHeight': '60px', 'borderWidth': '1px', 
    'borderStyle': 'dashed', 'borderRadius': '5px', 'textAlign': 'center', 
    'marginBottom': '10px', 'cursor': 'pointer', 'backgroundColor': '#f8f9fa'
}

def create_panel(side):
    return html.Div([
        # Upload
        dcc.Upload(id=f'upload-{side}', children=html.Div(['Drop/Click to Select File']), style=upload_style, accept='.pdb,.cif'),
        
        # 3D Viewer
        dash_molstar.MolstarViewer(id=f'viewer-{side}', style={'width': '100%', 'height': '500px'}),
        
        # Interactive Sequence Panel
        html.H4("Protein Sequence (Click to view)", style={'marginTop': '20px', 'fontFamily': 'Arial'}),
        html.Div(id=f'seq-{side}', style={
            'backgroundColor': '#f8f9fa', 'padding': '15px', 'maxHeight': '150px', 
            'overflowY': 'auto', 'fontFamily': 'monospace', 'border': '1px solid #dee2e6',
            'borderRadius': '5px', 'wordBreak': 'break-all'
        }),

        # Interaction Table Header & Sorting Form
        html.Div([
            html.H4("Interactions", style={'fontFamily': 'Arial', 'display': 'inline-block', 'marginRight': '20px'}),
            
            # --- SORTING FORM ---
            html.Div([
                html.Strong("Rank by position in: ", style={'fontSize': '14px'}),
                dcc.RadioItems(
                    id=f'sort-radio-{side}',
                    options=[
                        {'label': ' Column B (Res A)  ', 'value': 'A'},
                        {'label': ' Column D (Res B)  ', 'value': 'B'}
                    ],
                    value='A',
                    inline=True,
                    style={'display': 'inline-block', 'marginLeft': '10px', 'fontSize': '14px'}
                )
            ], style={'display': 'inline-block', 'backgroundColor': '#e9ecef', 'padding': '5px 10px', 'borderRadius': '5px'}),
            
            html.Button("📥 Excel", id=f'btn-excel-{side}', style={'float': 'right', 'marginTop': '5px'}),
        ], style={'width': '100%', 'marginTop': '20px'}),
        
        dcc.Download(id=f'download-excel-{side}'),
        
        dash_table.DataTable(
            id=f'table-{side}', 
            page_size=10, 
            style_table={'overflowX': 'auto'},
            style_cell={'cursor': 'pointer'},
            style_data_conditional=[{'if': {'state': 'active'}, 'backgroundColor': 'rgba(0, 116, 217, 0.3)', 'border': '1px solid blue'}]
        )
    ], style={'width': '48%', 'display': 'flex', 'flexDirection': 'column'})

app.layout = html.Div([
    html.H2("Parallel Molstar Analysis Workspace", style={'textAlign': 'center', 'fontFamily': 'Arial'}),
    html.Div([create_panel('left'), create_panel('right')], style={'display': 'flex', 'justifyContent': 'space-between', 'padding': '20px'}),
    
    # --- COMPARISON TABLE SECTION ---
    html.Div([
        html.Hr(style={'margin': '30px 0'}),
        html.Div([
            html.H3("Shared Residue Interactions (Left ⇄ Right)", style={'fontFamily': 'Arial', 'display': 'inline-block', 'marginRight': '20px'}),
            html.Div([
                html.Strong("Compare by: ", style={'fontSize': '14px'}),
                dcc.Dropdown(
                    id='compare-dropdown',
                    options=[
                        {'label': 'Column B (Res A)', 'value': 'B'},
                        {'label': 'Column D (Res B)', 'value': 'D'}
                    ],
                    value='B',
                    clearable=False,
                    style={'width': '200px', 'display': 'inline-block', 'marginLeft': '10px'}
                )
            ], style={'display': 'inline-block'}),
            html.Button("📥 Export Comparison", id='btn-export-comparison', style={'float': 'right', 'marginTop': '5px'}),
        ], style={'marginBottom': '15px'}),
        
        dcc.Download(id='download-comparison'),
        
        dash_table.DataTable(
            id='comparison-table',
            page_size=15,
            style_table={'overflowX': 'auto'},
            style_cell={'textAlign': 'left', 'padding': '8px', 'fontSize': '13px'},
            style_header={'backgroundColor': '#0074d9', 'color': 'white', 'fontWeight': 'bold'},
            style_data_conditional=[
                {'if': {'column_id': 'Matched Residue'}, 'backgroundColor': '#fff3cd', 'fontWeight': 'bold'}
            ]
        )
    ], style={'padding': '20px', 'backgroundColor': '#f8f9fa', 'borderRadius': '10px', 'margin': '20px'})
])


# --- CORE LOGIC CALLBACKS ---
def process_file_event(contents, filename, side, sort_val):
    if not contents: return None, "", []
    
    content_string = contents.split(',')[1]
    decoded = base64.b64decode(content_string)
    ext = os.path.splitext(filename)[1].lower() 
    
    with tempfile.NamedTemporaryFile(delete=False, suffix=ext) as tmp:
        tmp.write(decoded)
        tmp_path = tmp.name
        
    try:
        clean_path, sequence_ui, table_data = analyze_structure(tmp_path, ext, side, sort_val)
        molstar_data = molstar_helper.parse_molecule(clean_path)
    finally:
        if os.path.exists(tmp_path): os.remove(tmp_path)
        if 'clean_path' in locals() and os.path.exists(clean_path): os.remove(clean_path)
        
    return molstar_data, sequence_ui, table_data

def handle_click_events(active_cell, data, seq_type):
    """Handles both Sequence character clicks AND Data Table row clicks"""
    if not ctx.triggered:
        return no_update, no_update
        
    trigger_id = ctx.triggered_id
    
    if isinstance(trigger_id, dict) and trigger_id.get('type') == seq_type:
        chain = trigger_id['chain']
        res_num = trigger_id['res']
        
    elif isinstance(trigger_id, str) and trigger_id.startswith('table'):
        if not active_cell or not data: return no_update, no_update
        row = data[active_cell['row']]
        col_id = active_cell['column_id']
        
        if col_id in ['Column C (Chain B)', 'Column D (Res B)']:
            chain = row['Column C (Chain B)']
            res_str = row['Column D (Res B)']
        else:
            chain = row['Column A (Chain A)']
            res_str = row['Column B (Res A)']
            
        match = re.search(r'-?\d+', res_str)
        if not match: return no_update, no_update
        res_num = int(match.group())
    else:
        return no_update, no_update
        
    target = molstar_helper.get_targets(chain=chain, residue=res_num)
    return molstar_helper.get_focus(target), molstar_helper.get_selection(target)

def resort_table_data(sort_val, table_data):
    """Fast backend re-sorting avoiding full PDB re-parse"""
    if not table_data: return no_update
    df = pd.DataFrame(table_data)
    
    # Extract numeric position (e.g. A123 -> 123) for proper numeric sorting
    df['_sort_A'] = df['Column B (Res A)'].str.extract(r'(-?\d+)').astype(float)
    df['_sort_B'] = df['Column D (Res B)'].str.extract(r'(-?\d+)').astype(float)
    
    if sort_val == 'B':
        df = df.sort_values(by=[
            'Column C (Chain B)',
            '_sort_B',
            'Column A (Chain A)',
            '_sort_A',
            'Column E (Distance Å)'
        ])
    else:
        df = df.sort_values(by=[
            'Column A (Chain A)',
            '_sort_A',
            'Column C (Chain B)',
            '_sort_B',
            'Column E (Distance Å)'
        ])
        
    df = df.drop(columns=['_sort_A', '_sort_B'])
    return df.to_dict('records')


# --- LEFT PANEL CALLBACKS ---
@app.callback(
    Output('viewer-left', 'data'), Output('seq-left', 'children'), Output('table-left', 'data'),
    Input('upload-left', 'contents'), State('upload-left', 'filename'), State('sort-radio-left', 'value'), prevent_initial_call=True
)
def update_left(c, f, sort_val): return process_file_event(c, f, 'left', sort_val)

@app.callback(
    Output('table-left', 'data', allow_duplicate=True),
    Input('sort-radio-left', 'value'), State('table-left', 'data'), prevent_initial_call=True
)
def left_resort_click(sort_val, table_data): return resort_table_data(sort_val, table_data)

@app.callback(
    Output('viewer-left', 'focus'), Output('viewer-left', 'selection'),
    Input('table-left', 'active_cell'), 
    Input({'type': 'seq-aa-left', 'chain': ALL, 'res': ALL}, 'n_clicks'),
    State('table-left', 'data'), prevent_initial_call=True
)
def focus_left_viewer(cell, seq_clicks, data): return handle_click_events(cell, data, 'seq-aa-left')

@app.callback(
    Output('download-excel-left', 'data'),
    Input('btn-excel-left', 'n_clicks'), State('table-left', 'data'), prevent_initial_call=True
)
def download_left(n_clicks, table_data):
    if table_data: return dcc.send_data_frame(pd.DataFrame(table_data).to_excel, "Left_Interactions.xlsx", index=False)


# --- RIGHT PANEL CALLBACKS ---
@app.callback(
    Output('viewer-right', 'data'), Output('seq-right', 'children'), Output('table-right', 'data'),
    Input('upload-right', 'contents'), State('upload-right', 'filename'), State('sort-radio-right', 'value'), prevent_initial_call=True
)
def update_right(c, f, sort_val): return process_file_event(c, f, 'right', sort_val)

@app.callback(
    Output('table-right', 'data', allow_duplicate=True),
    Input('sort-radio-right', 'value'), State('table-right', 'data'), prevent_initial_call=True
)
def right_resort_click(sort_val, table_data): return resort_table_data(sort_val, table_data)

@app.callback(
    Output('viewer-right', 'focus'), Output('viewer-right', 'selection'),
    Input('table-right', 'active_cell'), 
    Input({'type': 'seq-aa-right', 'chain': ALL, 'res': ALL}, 'n_clicks'),
    State('table-right', 'data'), prevent_initial_call=True
)
def focus_right_viewer(cell, seq_clicks, data): return handle_click_events(cell, data, 'seq-aa-right')

@app.callback(
    Output('download-excel-right', 'data'),
    Input('btn-excel-right', 'n_clicks'), State('table-right', 'data'), prevent_initial_call=True
)
def download_right(n_clicks, table_data):
    if table_data: return dcc.send_data_frame(pd.DataFrame(table_data).to_excel, "Right_Interactions.xlsx", index=False)


# --- COMPARISON TABLE CALLBACKS ---
@app.callback(
    Output('comparison-table', 'data'),
    Input('table-left', 'data'),
    Input('table-right', 'data'),
    Input('compare-dropdown', 'value')
)
def update_comparison_table(left_data, right_data, compare_col):
    return compare_tables(left_data, right_data, compare_col)

@app.callback(
    Output('download-comparison', 'data'),
    Input('btn-export-comparison', 'n_clicks'),
    State('comparison-table', 'data'),
    prevent_initial_call=True
)
def download_comparison(n_clicks, comp_data):
    if comp_data:
        return dcc.send_data_frame(pd.DataFrame(comp_data).to_excel, "Shared_Residue_Interactions.xlsx", index=False)

if __name__ == '__main__':
    app.run(debug=True)
