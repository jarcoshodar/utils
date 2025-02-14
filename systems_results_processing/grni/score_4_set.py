import os
import sys
import pandas as pd
import networkx as nx
from pathlib import Path

def score_grn(input_csv, boolean_csv, active_tfs):

    df = pd.read_csv(input_csv)
    

    booldata = pd.read_csv(boolean_csv)
    bool_states = dict(zip(booldata.iloc[:, 0], booldata.iloc[:, 1]))
    
    #use netwrokx for the graph, remember here 1 is activation and -1, not 0 inhibition
    G = nx.DiGraph()
    effect_mapping = {'Activation': 1, 'Inhibition': -1}
    for _, row in df.iterrows():
        G.add_edge(row['TF'], row['target'], weight=effect_mapping[row['effect']])
    
    #subraph only of direct targets
    nodes_to_include = set(active_tfs)
    for tf in active_tfs:
        nodes_to_include.update(G.successors(tf))
    subG = G.subgraph(nodes_to_include).copy()
    
    def determine_state(node, node_states):
        if node in active_tfs:
            return 'active'
            
        parent_edges = subG.in_edges(node, data=True)
        active_activations = [
            edge_data['weight'] > 0
            for u, _, edge_data in parent_edges
            if node_states[u] == 'active'
        ]
        active_inhibitions = [
            edge_data['weight'] < 0
            for u, _, edge_data in parent_edges
            if node_states[u] == 'active'
        ]
    
        if any(active_inhibitions):
            return 'inactive'
        elif any(active_activations):
            return 'active'
        else:
            return node_states[node]
    
    def hash_node_states(node_states):
        return hash(frozenset(node_states.items()))
    
    def dfs_update_states(active_tfs, node_states):
        visited = set()
        stack = list(active_tfs)
        while stack:
            node = stack.pop()
            if node not in visited:
                visited.add(node)
                #only update the state if node is not an initial active TF
                if node not in active_tfs:
                    node_states[node] = determine_state(node, node_states)
                stack.extend([v for _, v in subG.out_edges(node)])
    
    #intilize all as inactive
    node_states = {n: 'inactive' for n in subG.nodes()}
    for tf in active_tfs:
        node_states[tf] = 'active'
    
    seen_states = set()
    previous_states = {}
    
    #prevent cycles with hash, ending early
    while hash_node_states(previous_states) != hash_node_states(node_states):
        current_state_hash = hash_node_states(node_states)
        if current_state_hash in seen_states:
            break  
        seen_states.add(current_state_hash)
        previous_states = node_states.copy()
        dfs_update_states(active_tfs, node_states)
    
    #estimate total boolean state explained
    matching_genes = sum(
        1 for gene, state in node_states.items()
        if gene in bool_states and (
            (state == 'active' and bool_states[gene] == 1) or
            (state == 'inactive' and bool_states[gene] == 0)
        )
    )
    score = matching_genes / len(bool_states) if len(bool_states) > 0 else 0
    

    for node in subG.nodes():
        subG.nodes[node]['state'] = node_states[node]
    
    return subG, score


if __name__ == '__main__':
    if len(sys.argv) < 4:  #changed from 5 to 4 since we need at least 1 TF
        print("Usage: python score_tfs.py <boolean_csv> <grni_csv> <output_dir> <TF1> [TF2] [TF3] [TF4]")
        print("Note: You can provide 1-4 transcription factors")
        sys.exit(1)
    
    boolean_csv = sys.argv[1]
    grni_csv = sys.argv[2]
    output_dir = Path(sys.argv[3])
    active_tfs = sys.argv[4:]  #non file, string seprated by spaces
    
    if len(active_tfs) > 4:
        print("Error: Please provide at most 4 transcription factors.")
        sys.exit(1)
    elif len(active_tfs) < 1:
        print("Error: Please provide at least 1 transcription factor.")
        sys.exit(1)

    output_dir.mkdir(exist_ok=True)

    if os.path.exists(grni_csv) and os.path.exists(boolean_csv):
        print(f"Processing GRNI file: {grni_csv}")
        print(f"Using boolean file: {boolean_csv}")
        print(f"Activating TFs: {', '.join(active_tfs)}")

        #socre for set
        subG, score = score_grn(grni_csv, boolean_csv, active_tfs)

        output_file = output_dir / f"Subnetwork_{Path(grni_csv).stem}_{'_'.join(active_tfs)}.csv"

        output_data = []
        for u, v, data in subG.edges(data=True):
            output_data.append({
                'TF': u,
                'effect': 'Activation' if data['weight'] > 0 else 'Inhibition',
                'target': v,
                'TF_state': subG.nodes[u]['state'],
                'target_state': subG.nodes[v]['state']
            })
        
        output_df = pd.DataFrame(output_data)
        output_df.to_csv(output_file, index=False)
        print(f"Saved subnetwork to {output_file}")
        print(f"Cumulative score for TFs {', '.join(active_tfs)}: {score}")

    else:
        print("Error: GRNI CSV file or boolean CSV file does not exist.")
        sys.exit(1)
