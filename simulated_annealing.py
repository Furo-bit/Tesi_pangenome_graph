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
# Simulated annealing

def simulated_annealing(block_dict: Dict, block_id_matrix: np.ndarray, params: Dict) -> Tuple[Dict, np.ndarray]:
    # Ask user inputs
    seed = int(params['seed'])
    random.seed(seed)
    temperature = int(params['temperature']) 
    termination_temperature = int(params['termination_temperature']) 
    cooling_factor = float(params['cooling_factor']) 
    threshold = int(params['threshold'])
    penalization = int(params['penalization'])


    rows = len(block_id_matrix)
    cols = len(block_id_matrix[0])
    cell_total = rows * cols

    while True:
        random_number = random.randint(1, cell_total)
        row_index = (random_number - 1) // cols
        col_index = (random_number - 1) % cols
        list_possible_moves = []
        block_1 = block_dict[block_id_matrix[row_index, col_index]]

        # Check left block
        if col_index-1 >= 0:
            block_2 = block_dict[block_id_matrix[row_index, col_index-1]] 
            if block_1["id"] != block_2["id"]:
                if block_1["sequence_ids"] == block_2["sequence_ids"]:
                    list_possible_moves += [[row_index, col_index-1]] 

        # Check right block
        if col_index+1 < cols:
            block_2 = block_dict[block_id_matrix[row_index, col_index+1]] 
            if block_1["id"] != block_2["id"]:
                if block_1["sequence_ids"] == block_2["sequence_ids"]:
                    list_possible_moves += [[row_index, col_index+1]] 

        # Check the column
        for row in range(rows):
            block_2 = block_dict[block_id_matrix[row, col_index]] 
            if block_1["id"] != block_2["id"]:
                if (block_1["begin_column"] == block_2["begin_column"] and
                    block_1["end_column"] == block_2["end_column"] and
                    block_1["label"] == block_2["label"]):
                    list_possible_moves += [[row, col_index]] 
        
        # Check column split
        if block_1["end_column"] - block_1["begin_column"] > 0:
            list_possible_moves += [["column_split"]]
        
        # Check row split
        if len(block_1["sequence_ids"]) > 1:
            list_possible_moves += [["row_split"]]
        
        if list_possible_moves != []:
            #print(list_possible_moves)
            n_rand = random.randint(0, len(list_possible_moves)-1)

            if list_possible_moves[n_rand][0] == "column_split" or list_possible_moves[n_rand][0] == "row_split":
                

                if list_possible_moves[n_rand][0] == "column_split":                        
                    # column split
                    new_block_1, new_block_2 = utils.split_block_by_column(block_1)

                elif list_possible_moves[n_rand][0] == "row_split":
                    # row split
                    new_block_1, new_block_2 = utils.split_block_by_row(block_1)

                label_b1 = block_1["label"]  
                label_nb1 = new_block_1["label"]
                label_nb2 = new_block_2["label"]  
                
                label_b1 = label_b1.translate(str.maketrans("", "", "-"))
                label_nb1 = label_nb1.translate(str.maketrans("", "", "-"))
                label_nb2 = label_nb2.translate(str.maketrans("", "", "-"))

                delta = (utils.of_min_label_length_threshold(threshold, penalization, len(label_nb1)) +
                        utils.of_min_label_length_threshold(threshold, penalization, len(label_nb2)) -
                        utils.of_min_label_length_threshold(threshold, penalization, len(label_b1))
                        )                  

                probability = utils.acceptance_probability(delta, temperature)

                if random.random() <= probability:

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
                row_to_merge = list_possible_moves[n_rand][0]
                column_to_merge = list_possible_moves[n_rand][1]
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
          
            # Cooling condition
            if list_possible_moves[n_rand][0] == "column_split" or list_possible_moves[n_rand][0] == "row_split":
                temperature = cooling_factor * temperature 
        #print("Temperature:", temperature)
        if temperature <= termination_temperature:
            break

    return block_dict, block_id_matrix

#----------------------------------------------------------------------------------------
def main(params_file: str, alignment_file: str, output_file: str, quality_file: str) -> None:


    # Reading the MSA
    sequences = utils.read_fasta(alignment_file) 

    # Convert sequences to matrix
    sequence_matrix = utils.sequences_to_matrix(sequences)    
    
    config = configparser.ConfigParser()
    config.read(params_file)
    params = config['parameters']
    start = str(params['start'])

    if start == "unitary_blocks":
        block_dict, block_id_matrix = utils.build_blocks_from_sequence_matrix(sequence_matrix)
    elif start == "greedy_rows":
        block_dict, block_id_matrix = utils.build_blocks_from_sequence_matrix(sequence_matrix)
        block_dict, block_id_matrix = utils.greedy_row_merge(block_dict, block_id_matrix)
    elif start == "sequence_blocks":
        block_dict, block_id_matrix = utils.build_long_blocks_from_sequence_matrix(sequence_matrix)
    else:
        raise ValueError(f"Unexpected value for start: {start}")

    # Simulated annealing
    block_dict, block_id_matrix = simulated_annealing(block_dict, block_id_matrix, params)
    

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