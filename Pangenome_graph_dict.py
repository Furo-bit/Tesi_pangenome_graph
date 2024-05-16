# Pangenome graph

#----------------------------------------------------------------------------------------
# Dependencies
from Bio import SeqIO
import networkx as nx
import matplotlib.pyplot as plt
import random
import copy
import numpy as np

#----------------------------------------------------------------------------------------
# Block
def create_block(id, label, sequence_id, begin_column, end_column):
    block_dict = {
        "id": id,
        "label": label,
        "sequence_ids": sequence_id if isinstance(sequence_id, list) else [sequence_id],
        "begin_column": begin_column,
        "end_column": end_column
    }
    return block_dict
#----------------------------------------------------------------------------------------
# Utility functions

# Read msa file
def read_fasta(filename):
    sequences = {}
    for record in SeqIO.parse(filename, "fasta"):
        sequences[record.id] = list(str(record.seq))  # Convert sequence to a list of characters
    return sequences

# Put the sequences in a matrix
def sequences_to_matrix(sequences):
    max_length = max(len(seq) for seq in sequences.values())
    matrix = []
    for sequence_id, sequence in sequences.items():
        padded_sequence = sequence + [''] * (max_length - len(sequence))  # Pad shorter sequences with empty strings
        matrix.append(padded_sequence)
    return matrix

# From labels to block dict and block id matrix
def build_blocks_from_sequence_matrix(sequence_matrix):
    block_dict = {}
    num_seq = len(sequence_matrix)
    num_bases = len(sequence_matrix[0])    
    block_id_matrix = np.empty((num_seq, num_bases), dtype=object)

    for sequence_index in range(num_seq):
        for base_index in range(num_bases):
            block_id = "B" + str(sequence_index) + "." + str(base_index)
            block_dict[block_id] = create_block(
                block_id,
                sequence_matrix[sequence_index][base_index],
                sequence_index,
                base_index,
                base_index
            )
            block_id_matrix[sequence_index, base_index] = block_id
    
    return block_dict, block_id_matrix

# Merge two blocks in one
def merge_two_blocks(block_1, block_2, how):
    if how == "column_union":
        new_labels = block_1["label"] + block_2["label"]
        new_block = {
            "id": block_1["id"],
            "label": new_labels,
            "sequence_ids": block_1["sequence_ids"],
            "begin_column": block_1["begin_column"],
            "end_column": block_2["end_column"]
        }
        print(block_1["id"], "with", block_2["id"], "by columns, new label:", new_labels)
    elif how == "row_union":
        new_sequence_ids = block_1["sequence_ids"] + block_2["sequence_ids"]
        new_block = {
            "id": block_1["id"],
            "label": block_1["label"],
            "sequence_ids": new_sequence_ids,
            "begin_column": block_1["begin_column"],
            "end_column": block_1["end_column"]
        }
        print(block_1["id"], "with", block_2["id"], "by columns, new sequences:", new_sequence_ids)
    return new_block


# Update block dict data
def update_block_dict_with_same_id(block_dict, new_block, id1, id2):
    for block_id, block_data in block_dict.items():
        if block_id == id1 or block_id == id2:
            block_dict[block_id] = new_block
    return block_dict

# Update block matrix data
def update_block_matrix_with_same_id(block_id_matrix, new_id, id1, id2):
    rows = len(block_id_matrix)
    cols = len(block_id_matrix[0])
    for i in range(rows):
        for j in range(cols):  
            if block_id_matrix[i, j] == id1 or block_id_matrix[i, j] == id2:
                block_id_matrix[i, j] = new_id
 
    return block_id_matrix

'''
# From dictionary to graph
def build_graph(block_dict, block_id_matrix):
    G = nx.DiGraph()
    rows = len(block_id_matrix)
    cols = len(block_id_matrix[0])

    # Add nodes
    G.add_node("source", label = "source")
    G.add_node("sink", label = "sink")
    for i in range(rows):
        for j in range(cols):
            G.add_node(block_id_matrix[i][j].id, label=matrix[i][j].base, row=i, col=j)
            
    # Calculate positions for nodes based on "row" and "col" attributes
    pos = {}
    for node, data in G.nodes(data=True):
        if node == "source":
            pos[node] = (-1, 0)
        elif node == "sink":
            pos[node] = (cols + 1, 0)
        else:
            pos[node] = (data['col'], -data['row'])  # Invert row value for upward display
       
    # Add edges
    for i in range(rows):
        G.add_edge("source", matrix[i][0].id, label = matrix[i][0].sequences_ids)
        G.add_edge(matrix[i][cols-1].id, "sink", label = matrix[i][cols-1].sequences_ids)
        for j in range(cols):
            # Connect horizontally adjacent nodes
            if j < cols - 1 and matrix[i][j].id != matrix[i][j+1].id:
                G.add_edge(matrix[i][j].id,
                            matrix[i][j+1].id,
                             label = list(set(matrix[i][j].sequences_ids).intersection(set(matrix[i][j+1].sequences_ids))) 
                        )
    
    return G, pos
'''
#----------------------------------------------------------------------------------------
# Euristichs
def local_search(block_dict, block_id_matrix):
    rows, cols = block_id_matrix.shape
    
    for i in range(rows-1):
        for j in range(cols):  
            # Try to merge two blocks by rows
            if i+1 != rows:
                block_1 = block_dict[block_id_matrix[i, j]]
                block_2 = block_dict[block_id_matrix[i+1, j]] 
                if (block_1["begin_column"] == block_2["begin_column"] and
                    block_1["end_column"] == block_2["end_column"] and
                    block_1["label"] == block_2["label"]):
                    new_block = merge_two_blocks(block_1, block_2, "row_union")
                    block_dict = update_block_dict_with_same_id(block_dict, new_block, block_1["id"], block_2["id"])
                    # Devo aggiornare tutta la matrice di id, non solo le due posizioni mergiate
                    block_id_matrix = update_block_matrix_with_same_id(block_id_matrix, new_block["id"], block_1["id"], block_2["id"])
                    #block_id_matrix[i, j] = new_block["id"]
                    #block_id_matrix[i+1, j] = new_block["id"]
                    del block_dict[block_2["id"]]
                    
            # Try to merge two blocks by columns
            if j+1 != cols:
                block_1 = block_dict[block_id_matrix[i, j]]
                block_2 = block_dict[block_id_matrix[i, j+1]]
                if block_1["sequence_ids"] == block_2["sequence_ids"]:
                    new_block = merge_two_blocks(block_1, block_2, "column_union")
                    block_dict = update_block_dict_with_same_id(block_dict, new_block, block_1["id"], block_2["id"])
                    block_id_matrix = update_block_matrix_with_same_id(block_id_matrix, new_block["id"], block_1["id"], block_2["id"])
                    #block_id_matrix[i, j] = new_block["id"]
                    #block_id_matrix[i, j+1] = new_block["id"]
                    del block_dict[block_2["id"]]

    return block_dict, block_id_matrix, "ls"




#----------------------------------------------------------------------------------------

def main():
    # Reading the MSA
    filename = "test"  # Change this to your .fa file
    sequences = read_fasta("Test_allignments/"+filename+".fa")

    # Convert sequences to matrix
    sequence_matrix = sequences_to_matrix(sequences) 

    # Convert each label in a block
    block_dict, block_id_matrix = build_blocks_from_sequence_matrix(sequence_matrix)
    
    # Euristich
    block_dict, block_id_matrix, euristich = local_search(block_dict, block_id_matrix)
    print(block_dict, block_id_matrix)


if __name__ == "__main__":
    main()