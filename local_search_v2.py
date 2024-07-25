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

    threshold = int(params['threshold'])
    penalization = int(params['penalization'])
    seed = int(params['seed'])
    rows = len(block_id_matrix)
    cols = len(block_id_matrix[0])
    cell_total = rows * cols
    limit = cols
    # Generate a random numbers using the user-provided seed
    random_numbers = utils.generate_random_numbers(seed, 1, cell_total, int(cell_total/3))
    
    for x in random_numbers:
        # Transform the number in matrix indeces 
        row_index = (x - 1) // cols
        col_index = (x - 1) % cols
        block_1 = block_dict[block_id_matrix[row_index, col_index]]

        list_mergeable_horizontal = []
        list_mergeable_horizontal_left = []
        list_mergeable_horizontal_right = []
        list_mergeable_vertical = []
        operation_gains = {}
        block_1_label = block_1["label"]
        block_1_label = block_1_label.translate(str.maketrans("", "", "-"))
        block_1_label_value = utils.of_pangeblocks(threshold, penalization, len(block_1_label))    

        # Check the column for merge
        for row in range(rows):
            block_2 = block_dict[block_id_matrix[row, col_index]] 
            if block_1["id"] != block_2["id"] and block_2["id"] not in list_mergeable_vertical:
                if (block_1["begin_column"] == block_2["begin_column"] and
                    block_1["end_column"] == block_2["end_column"] and
                    block_1["label"] == block_2["label"]):
                    list_mergeable_vertical += [block_2["id"]]
        
        # Calculate gain
        vertical_merge_gain = block_1_label_value
        new_block_vertical = block_1
        if list_mergeable_vertical != [] :
            for id in list_mergeable_vertical:
                block_2 = block_dict[id]
                block_2_label = block_2["label"]
                block_2_label = block_2_label.translate(str.maketrans("", "", "-"))
                block_2_label_value = utils.of_pangeblocks(threshold, penalization, len(block_2_label))    
                vertical_merge_gain += block_2_label_value
                new_block_vertical = utils.merge_two_blocks(new_block_vertical, block_2, "row_union")
            
            new_block_vertical_label = new_block_vertical["label"]
            new_block_vertical_label = new_block_vertical_label.translate(str.maketrans("", "", "-"))
            new_block_vertical_label_value = utils.of_pangeblocks(threshold, penalization, len(new_block_vertical_label))    
            vertical_merge_gain = new_block_vertical_label_value - vertical_merge_gain
            operation_gains["row_merge"] = vertical_merge_gain


        # Check the row for merge
        # Check left
        if block_1["begin_column"] <= limit :
            limit_left = 0
        else:
            limit_left = block_1["begin_column"] - limit
        for column in reversed(range(limit_left, block_1["begin_column"])):
            block_2 = block_dict[block_id_matrix[row_index, column]] 
            if block_1["id"] != block_2["id"] and block_2["id"] not in list_mergeable_horizontal_left:
                if block_1["sequence_ids"] == block_2["sequence_ids"]:
                    list_mergeable_horizontal_left += [block_2["id"]]
                else:
                    break
        list_mergeable_horizontal_left.reverse()
        
        # Check right
        if block_1["end_column"] + limit >= cols :
            limit_right = cols
        else:
            limit_right = block_1["end_column"] + limit
        for column in range(block_1["end_column"], limit_right):
            block_2 = block_dict[block_id_matrix[row_index, column]] 
            if block_1["id"] != block_2["id"] and block_2["id"] not in list_mergeable_horizontal_right:
                if block_1["sequence_ids"] == block_2["sequence_ids"]:
                    list_mergeable_horizontal_right += [block_2["id"]]
                else:
                    break
        
        # Calculate gain
        list_mergeable_horizontal = list_mergeable_horizontal_left + [block_1["id"]] + list_mergeable_horizontal_right
        horizontal_merge_gain = block_1_label_value
        if list_mergeable_horizontal_left != [] or list_mergeable_horizontal_right != []:
            new_block_horizontal = block_dict[list_mergeable_horizontal[0]]
            #del list_mergeable_horizontal[0]
            for id in list_mergeable_horizontal:
                if new_block_horizontal["id"] != id :
                    block_2 = block_dict[id]
                    block_2_label = block_2["label"]
                    block_2_label = block_2_label.translate(str.maketrans("", "", "-"))
                    block_2_label_value = utils.of_pangeblocks(threshold, penalization, len(block_2_label))    
                    horizontal_merge_gain += block_2_label_value
                    new_block_horizontal = utils.merge_two_blocks(new_block_horizontal, block_2, "column_union")

            new_block_horizontal_label = new_block_horizontal["label"]
            new_block_horizontal_label = new_block_horizontal_label.translate(str.maketrans("", "", "-"))
            new_block_horizontal_label_value = utils.of_pangeblocks(threshold, penalization, len(new_block_horizontal_label))    
            horizontal_merge_gain = new_block_horizontal_label_value - horizontal_merge_gain
            operation_gains["column_merge"] = horizontal_merge_gain



        # Check row split
        if len(block_1["sequence_ids"]) > 1:

            new_block_1r, new_block_2r = utils.split_block_by_row(block_1)

            label_nb1r = new_block_1r["label"]
            label_nb2r = new_block_2r["label"]  
            
            label_nb1r = label_nb1r.translate(str.maketrans("", "", "-"))
            label_nb2r = label_nb2r.translate(str.maketrans("", "", "-"))

            delta = (utils.of_pangeblocks(threshold, penalization, len(label_nb1r)) +
                    utils.of_pangeblocks(threshold, penalization, len(label_nb2r)) -
                    block_1_label_value
                    ) 
            operation_gains["row_split"] = delta
        
        # Check column split
        if block_1["end_column"] - block_1["begin_column"] > 0:

            new_block_1c, new_block_2c = utils.split_block_by_column(block_1)

            label_nb1c = new_block_1c["label"]
            label_nb2c = new_block_2c["label"]  
            
            label_nb1c = label_nb1c.translate(str.maketrans("", "", "-"))
            label_nb2c = label_nb2c.translate(str.maketrans("", "", "-"))

            delta = (utils.of_pangeblocks(threshold, penalization, len(label_nb1c)) +
                    utils.of_pangeblocks(threshold, penalization, len(label_nb2c)) -
                    block_1_label_value
                    ) 
            operation_gains["column_split"] = delta
        
        # Select operation
        
        if operation_gains != {}:
            operation_gains = dict(sorted(operation_gains.items(), key=lambda item: item[1]))
            operation, gain = next(iter(operation_gains.items()))

            if operation == "row_merge":
                # row merge
                for id in list_mergeable_vertical:
                    del block_dict[id]
                block_dict[new_block_vertical["id"]] = new_block_vertical  

                for sequence in new_block_vertical["sequence_ids"]:
                    for column in range(new_block_vertical["begin_column"], new_block_vertical["end_column"] + 1):
                        block_id_matrix[sequence, column] = new_block_vertical["id"]   

            elif operation == "column_merge":
                # column merge
                for id in list_mergeable_horizontal:
                    del block_dict[id]
                block_dict[new_block_horizontal["id"]] = new_block_horizontal 

                for sequence in new_block_horizontal["sequence_ids"]:
                    for column in range(new_block_horizontal["begin_column"], new_block_horizontal["end_column"] + 1):
                        block_id_matrix[sequence, column] = new_block_horizontal["id"] 

            elif operation == "row_split":
                # row split
                # update new_block_1
                for sequence in new_block_1r["sequence_ids"]:
                    for column in range(new_block_1r["begin_column"], new_block_1r["end_column"] + 1):
                        block_id_matrix[sequence, column] = new_block_1r["id"]
                # update new_block_2
                for sequence in new_block_2r["sequence_ids"]:
                    for column in range(new_block_2r["begin_column"], new_block_2r["end_column"] + 1):
                        block_id_matrix[sequence, column] = new_block_2r["id"]   
                # update dict data
                del block_dict[block_1["id"]]
                block_dict[new_block_1r["id"]] = new_block_1r
                block_dict[new_block_2r["id"]] = new_block_2r  

            elif operation == "column_split":
                # column split
                # update new_block_1
                for sequence in new_block_1c["sequence_ids"]:
                    for column in range(new_block_1c["begin_column"], new_block_1c["end_column"] + 1):
                        block_id_matrix[sequence, column] = new_block_1c["id"]
                # update new_block_2
                for sequence in new_block_2c["sequence_ids"]:
                    for column in range(new_block_2c["begin_column"], new_block_2c["end_column"] + 1):
                        block_id_matrix[sequence, column] = new_block_2c["id"]   
                # update dict data
                del block_dict[block_1["id"]]
                block_dict[new_block_1c["id"]] = new_block_1c
                block_dict[new_block_2c["id"]] = new_block_2c  
  
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

    
    # Local search
    block_dict, block_id_matrix = local_search(block_dict, block_id_matrix, params)

    # Graph
    graph, pos = utils.build_graph(block_dict, block_id_matrix)
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


    
    # Draw the graph 
    nx.draw(graph, pos=pos, labels=nx.get_node_attributes(graph, 'label'), with_labels=True)
    nx.draw_networkx_edge_labels(graph, pos=pos, edge_labels=nx.get_edge_attributes(graph, 'label'), font_color='red')

    # Show the graph
    plt.show()
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--params', required=True, help='Parameters file (.txt)')
    parser.add_argument('--input', required=True, help='Alignment file, input file (.fa)')
    parser.add_argument('--output', required=True, help='Output file (.gfa)')
    parser.add_argument('--quality', required=True, help='Objective function value output file (.txt)')

    args = parser.parse_args()

    main(args.params, args.input, args.output, args.quality)