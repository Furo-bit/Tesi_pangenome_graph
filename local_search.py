# Dependencies
from Bio import SeqIO
import networkx as nx
import matplotlib.pyplot as plt
import random
import numpy as np
import math
from typing import List, Dict, Tuple
import argparse
import configparser

import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'Utils')))
import utils

#----------------------------------------------------------------------------------------
# Local search

def local_search(block_dict: Dict, block_id_matrix: np.ndarray, params: Dict) -> Tuple[Dict, np.ndarray]:
     # Ask user input
    seed = int(params['seed'])

    rows = len(block_id_matrix)
    cols = len(block_id_matrix[0])
    cell_total = rows * cols
    # Generate a random numbers using the user-provided seed
    random_numbers = utils.generate_random_numbers(seed, 1, cell_total, int(cell_total*2))
    
    for x in random_numbers:
        # Transform the number in matrix indeces 
        row_index = (x - 1) // cols
        col_index = (x - 1) % cols
        list_mergeable = []
        block_1 = block_dict[block_id_matrix[row_index, col_index]]

        # Check left block
        if col_index-1 >= 0:
            block_2 = block_dict[block_id_matrix[row_index, col_index-1]] 
            if block_1["id"] != block_2["id"]:
                if block_1["sequence_ids"] == block_2["sequence_ids"]:
                    list_mergeable += [[row_index, col_index-1]] 

        # Check right block
        if col_index+1 < cols:
            block_2 = block_dict[block_id_matrix[row_index, col_index+1]] 
            if block_1["id"] != block_2["id"]:
                if block_1["sequence_ids"] == block_2["sequence_ids"]:
                    list_mergeable += [[row_index, col_index+1]] 

        # Check the column
        for row in range(rows):
            block_2 = block_dict[block_id_matrix[row, col_index]] 
            if block_1["id"] != block_2["id"]:
                if (block_1["begin_column"] == block_2["begin_column"] and
                    block_1["end_column"] == block_2["end_column"] and
                    block_1["label"] == block_2["label"]):
                    list_mergeable += [[row, col_index]] 
        
        # Check column split
        if block_1["end_column"] - block_1["begin_column"] > 0:
            list_mergeable += [["column_split"]]
        
        # Check row split
        if len(block_1["sequence_ids"]) > 1:
            list_mergeable += [["row_split"]]

        # Select a random mergeable block
        if list_mergeable != []:
            random.seed(seed)
            n_rand = random.randint(0, len(list_mergeable)-1)

            if list_mergeable[n_rand][0] == "column_split" or list_mergeable[n_rand][0] == "row_split":
                
                if list_mergeable[n_rand][0] == "column_split":                        
                    # column split
                    new_block_1, new_block_2 = utils.split_block_by_column(block_1)
                
                elif list_mergeable[n_rand][0] == "row_split":
                    # row split
                    new_block_1, new_block_2 = utils.split_block_by_row(block_1)
                
                # update new_block_1
                for sequence in new_block_1["sequence_ids"]:
                    for column in range(new_block_1["begin_column"], new_block_1["end_column"] + 1):
                        block_id_matrix[sequence, column] = new_block_1["id"]
                # update new_block_2
                for sequence in new_block_2["sequence_ids"]:
                    for column in range(new_block_2["begin_column"], new_block_2["end_column"] + 1):
                        block_id_matrix[sequence, column] = new_block_2["id"]
                # update dict data
                del block_dict[block_1["id"]]
                block_dict[new_block_1["id"]] = new_block_1
                block_dict[new_block_2["id"]] = new_block_2
                               
            else:
                row_to_merge = list_mergeable[n_rand][0]
                column_to_merge = list_mergeable[n_rand][1]
                if row_to_merge == row_index:
                    how = "column_union"
                else:
                    how = "row_union"
                block_2 = block_dict[block_id_matrix[row_to_merge][column_to_merge]]
                # Update
                if row_to_merge < row_index or column_to_merge < col_index:
                    new_block = utils.merge_two_blocks(block_2, block_1, how)
                    block_dict = utils.update_block_dict_with_same_id(block_dict, new_block, block_2["id"], block_1["id"])
                    del block_dict[block_1["id"]]
                else: 
                    new_block = utils.merge_two_blocks(block_1, block_2, how)
                    block_dict = utils.update_block_dict_with_same_id(block_dict, new_block, block_1["id"], block_2["id"])
                    del block_dict[block_2["id"]]
                for sequence in new_block["sequence_ids"]:
                    for column in range(new_block["begin_column"], new_block["end_column"] + 1):
                        block_id_matrix[sequence, column] = new_block["id"]
            
    return block_dict, block_id_matrix

#----------------------------------------------------------------------------------------
def main(params_file: str, alignment_file: str, output_file: str, quality_file: str) -> None:

    # Reading the MSA
    sequences = utils.read_fasta(alignment_file) 

    # Convert sequences to matrix
    sequence_matrix = utils.sequences_to_matrix(sequences)    
    
    # Convert each label in a block
    block_dict, block_id_matrix = utils.build_blocks_from_sequence_matrix(sequence_matrix)
    
    # Local search
    config = configparser.ConfigParser()
    config.read(params_file)
    params = config['parameters']
    block_dict, block_id_matrix = local_search(block_dict, block_id_matrix, params)

    # Graph
    graph = utils.build_graph(block_dict, block_id_matrix)
    graph = utils.remove_indle_from_graph(graph)
      
    # Computing objective function
    objective_value = 0
    for node, data in graph.nodes(data=True):
        label = data.get('label', '')
        objective_value += utils.of_pangeblocks(int(params['threshold']), int(params['penalization']), len(label))
    
    # Save graph
    utils.graph_to_gfa(graph, output_file)

    # Save the objective function value
    with open(quality_file, 'w') as f:
        f.write(str(objective_value))


    '''
    # Draw the graph 
    nx.draw(graph, pos=pos, labels=nx.get_node_attributes(graph, 'label'), with_labels=True)
    nx.draw_networkx_edge_labels(graph, pos=pos, edge_labels=nx.get_edge_attributes(graph, 'label'), font_color='red')

    # Show the graph
    plt.show()
    '''

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--params', required=True, help='Parameters file (.txt)')
    parser.add_argument('--input', required=True, help='Alignment file, input file (.fa)')
    parser.add_argument('--output', required=True, help='Output file (.gfa)')
    parser.add_argument('--quality', required=True, help='Objective function value output file (.txt)')

    args = parser.parse_args()

    main(args.params, args.input, args.output, args.quality)