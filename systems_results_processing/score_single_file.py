import os
import pandas as pd
import sys
import networkx as nx
from pathlib import Path

def score_grn(input_csv, boolean_csv):
    # Read input CSV
    df = pd.read_csv(input_csv)
    
    # Read boolean data
    # Note: This assumes the boolean CSV has gene names in the first column
    # and boolean values in the second column, regardless of column names
    booldata = pd.read_csv(boolean_csv)
    bool_states = dict(zip(booldata.iloc[:, 0], booldata.iloc[:, 1]))

    # Create graph
    G = nx.DiGraph()
    effect_mapping = {'Activation': 1, 'Inhibition': -1}
    for _, row in df.iterrows():
        G.add_edge(row['TF'], row['target'], weight=effect_mapping[row['effect']])

    def determine_state(node, node_states):
        parent_edges = G.in_edges(node, data=True)
        active_activations = [edge_data['weight'] > 0 for u, _, edge_data in parent_edges if node_states[u] == 'active']
        active_inhibitions = [edge_data['weight'] < 0 for u, _, edge_data in parent_edges if node_states[u] == 'active']
        
        if any(active_inhibitions):
            return 'inactive'
        elif any(active_activations):
            return 'active'
        else:
            return node_states[node]

    def hash_node_states(node_states):
        return hash(frozenset(node_states.items()))

    def dfs_update_states(start_node, node_states, visited=None):
        if visited is None:
            visited = set()
        visited.add(start_node)
        
        for _, child_node, _ in G.out_edges(start_node, data=True):
            if child_node not in visited:
                node_states[child_node] = determine_state(child_node, node_states)
                dfs_update_states(child_node, node_states, visited)

    node_scores = {}
    for node in G.nodes:
        node_states = {n: 'neutral' for n in G.nodes()}
        node_states[node] = 'active'
        
        seen_states = set()
        previous_states = {}
        
        while hash_node_states(previous_states) != hash_node_states(node_states):
            current_state_hash = hash_node_states(node_states)
            if current_state_hash in seen_states:
                break
            
            seen_states.add(current_state_hash)
            previous_states = node_states.copy()
            dfs_update_states(node, node_states)
        
        score = sum(1 for gene, state in node_states.items() 
                    if gene in bool_states and 
                    ((state == 'active' and bool_states[gene] == 1) or 
                     (state == 'inactive' and bool_states[gene] == 0))) / len(node_states)
        
        node_scores[node] = score

    # Add scores to dataframe
    df['Score'] = df['TF'].map(node_scores)

    return df

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python score_single_file.py <path_to_GRNI_csvs> <path_to_bool_csvs> <output_dir>")
        sys.exit(1)
    
    grni_dir = Path(sys.argv[1])
    bool_dir = Path(sys.argv[2])
    output_dir = Path(sys.argv[3])
    
    print("Note: This script assumes that in the boolean CSV file, gene names are in the first column")
    print("and boolean values (0 or 1) are in the second column, regardless of column names.")
    
    #outputdircreated if it doesn't exist
    output_dir.mkdir(exist_ok=True)
    
    #get all .csv in grni_dir
    grni_files = list(grni_dir.glob("*.csv"))
    
    for grni_file in grni_files:
        #generate a boolean file name
        #e.g., if your GRNI is "GRNI5_NeurpOE_vs_NeurpOC.csv"
        #and boolean "bool_NeurpOE_vs_NeurpOC.csv" this will work
        bool_file_name = grni_file.name.replace("GRNI5_", "bool_")
        
        bool_file = bool_dir / bool_file_name
        
        if bool_file.exists():
            print(f"Processing {grni_file.name} with boolean {bool_file_name}")
            result = score_grn(str(grni_file), str(bool_file))
            
            
            output_file = output_dir / f"Scored_{grni_file.name}"
            result.to_csv(output_file, index=False)
            print(f"Saved scored result to {output_file}")
        else:
            print(f"Warning: No matching boolean file found for {grni_file.name} (looked for {bool_file_name})")
    
    print("Processing complete. Scored GRNIs are saved in the output directory.")

