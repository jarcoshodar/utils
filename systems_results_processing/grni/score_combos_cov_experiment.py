import os
import pandas as pd
import networkx as nx
import itertools
from pathlib import Path
import sys

def score_grn(input_csv, boolean_csv, active_tfs):

    df = pd.read_csv(input_csv)
    

    booldata = pd.read_csv(boolean_csv)
    bool_states = dict(zip(booldata.iloc[:, 0], booldata.iloc[:, 1]))

    #use netwrokx for the graph, remember here 1 is activation and -1, not 0 inhibition
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
    
    def dfs_update_states(active_tfs, node_states):
        visited = set()
        stack = list(active_tfs)
        while stack:
            node = stack.pop()
            if node not in visited:
                visited.add(node)
                node_states[node] = determine_state(node, node_states)
                stack.extend([v for _, v in G.out_edges(node)])
    
    #intilize all as inactive
    node_states = {n: 'inactive' for n in G.nodes()}
    for tf in active_tfs:
        node_states[tf] = 'active'
    
    seen_states = set()
    previous_states = {}
    
    #prevent cycles with hash, ending early
    while hash_node_states(previous_states) != hash_node_states(node_states):
        current_state_hash = hash_node_states(node_states)
        if current_state_hash in seen_states:
            break  # Cycle detected
        seen_states.add(current_state_hash)
        previous_states = node_states.copy()
        dfs_update_states(active_tfs, node_states)
    
    #estimate total boolean state explained
    matching_genes = sum(1 for gene, state in node_states.items()
                         if gene in bool_states and 
                         ((state == 'active' and bool_states[gene] == 1) or 
                          (state == 'inactive' and bool_states[gene] == 0)))
    score = matching_genes / len(bool_states)
    
    return score


if __name__ == '__main__':
    if len(sys.argv) != 5:  #parameter with combination size, to check all
        print("Usage: python score_combos.py <boolean_csv> <grni_csv> <output_dir> <combo_size>")
        print("combo_size should be 2 or 3")
        sys.exit(1)

    boolean_csv = sys.argv[1]
    grni_csv = sys.argv[2]
    output_dir = Path(sys.argv[3])
    combo_size = int(sys.argv[4])  
    
    if combo_size not in [2, 3]:
        print("Error: combo_size must be 2 or 3")
        sys.exit(1)

    output_dir.mkdir(exist_ok=True)
    
    if os.path.exists(grni_csv) and os.path.exists(boolean_csv):
        print(f"Processing GRNI file: {grni_csv}")
        print(f"Using boolean file: {boolean_csv}")
        
        #HARCODED WARNING: this will be versatile later, just quick fix
        tf_list = ["PBX1", "FEZF2", "NR3C2", "L3MBTL4", "ETV5"]
        
        #all allcombinations of that size
        tf_combinations = list(itertools.combinations(tf_list, combo_size))
        total_combinations = len(tf_combinations)
        print(f"Total combinations to process: {total_combinations}")

        results = []
        update_interval = 10  
        counter = 0

        for combination in tf_combinations:
            score = score_grn(grni_csv, boolean_csv, combination)
            results.append((combination, score))
            counter += 1

            if counter % update_interval == 0 or counter == total_combinations:
                print(f"Processed {counter}/{total_combinations} combinations "
                      f"({(counter/total_combinations)*10:.2f}%)")

        
        results.sort(key=lambda x: x[1], reverse=True)
        output_file = output_dir / f"Ranked_Combinations_{Path(grni_csv).stem}.csv"

        df_results = pd.DataFrame({
            'Combination': [';'.join(combo) for combo, _ in results],
            'Score': [score for _, score in results]
        })
        df_results.to_csv(output_file, index=False)
        print(f"Saved ranked combinations to {output_file}")

    else:
        print("Error: GRNI CSV file or boolean CSV file does not exist.")
        sys.exit(1)
