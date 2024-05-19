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
        #print(block_1["id"], "with", block_2["id"], "by columns, new label:", new_labels)
    elif how == "row_union":
        new_sequence_ids = block_1["sequence_ids"] + block_2["sequence_ids"]
        new_block = {
            "id": block_1["id"],
            "label": block_1["label"],
            "sequence_ids": new_sequence_ids,
            "begin_column": block_1["begin_column"],
            "end_column": block_1["end_column"]
        }
        #print(block_1["id"], "with", block_2["id"], "by rows, new sequences:", new_sequence_ids)
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

def graph_to_gfa(graph, filename):
    with open(filename, 'w') as f:
        f.write("H\tVN:Z:1.0\n")  # GFA header line
        for node, data in graph.nodes(data=True):
            f.write(f"S\t{node}\t*\tLN:i:{data['label']}\n")  # S-line for each node
        for u, v, data in graph.edges(data=True):
            edge_label = data.get('label', '*')  # Get the edge label, or use '*' if not present
            f.write(f"L\t{u}\t+\t{v}\t+\t{edge_label}\n")  # L-line for each edge with label

# Generate count numbers from seed
def generate_random_numbers(seed, start, end, count):
    random.seed(seed)
    return [random.randint(start, end) for _ in range(count)]

def build_graph(block_dict, block_id_matrix):
    G = nx.DiGraph()
    rows, cols = block_id_matrix.shape

    # Add nodes
    G.add_node("source", label="source")
    G.add_node("sink", label="sink")
    for i in range(rows):
        for j in range(cols):
            block_id = block_id_matrix[i, j]
            block = block_dict[block_id]
            G.add_node(block_id, label=block["label"], row=i, col=j)
            
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
        G.add_edge("source", block_id_matrix[i, 0], label=block_dict[block_id_matrix[i, 0]]["sequence_ids"])
        G.add_edge(block_id_matrix[i, cols-1], "sink", label=block_dict[block_id_matrix[i, cols-1]]["sequence_ids"])
        for j in range(cols - 1):
            current_block_id = block_id_matrix[i, j]
            next_block_id = block_id_matrix[i, j+1]
            if current_block_id != next_block_id:
                current_block = block_dict[current_block_id]
                next_block = block_dict[next_block_id]
                common_sequences = list(set(current_block["sequence_ids"]).intersection(set(next_block["sequence_ids"])))
                G.add_edge(current_block_id, next_block_id, label=common_sequences)

    return G, pos
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
                    block_id_matrix = update_block_matrix_with_same_id(block_id_matrix, new_block["id"], block_1["id"], block_2["id"])
                    del block_dict[block_2["id"]]
                    
            # Try to merge two blocks by columns
            if j+1 != cols:
                block_1 = block_dict[block_id_matrix[i, j]]
                block_2 = block_dict[block_id_matrix[i, j+1]]
                if block_1["sequence_ids"] == block_2["sequence_ids"]:
                    new_block = merge_two_blocks(block_1, block_2, "column_union")
                    block_dict = update_block_dict_with_same_id(block_dict, new_block, block_1["id"], block_2["id"])
                    block_id_matrix = update_block_matrix_with_same_id(block_id_matrix, new_block["id"], block_1["id"], block_2["id"])
                    del block_dict[block_2["id"]]

    return block_dict, block_id_matrix, "ls"


def local_search_random_2(block_dict, block_id_matrix):
     # Ask user input
    seed = int(input("Insert the seed for the random number: "))

    rows = len(block_id_matrix)
    cols = len(block_id_matrix[0])
    cell_total = rows * cols
    # Generate a random numbers using the user-provided seed
    random_numbers = generate_random_numbers(seed, 1, cell_total, cell_total*2)
    
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

        # Select a random mergeable block
        if list_mergeable != []:
            random.seed(seed)
            n_rand = random.randint(0, len(list_mergeable)-1)
            row_to_merge = list_mergeable[n_rand][0]
            column_to_merge = list_mergeable[n_rand][1]
            if row_to_merge == row_index:
                how = "column_union"
            else:
                how = "row_union"
            block_2 = block_dict[block_id_matrix[row_to_merge][column_to_merge]]
            # Update 
            if row_to_merge < row_index or column_to_merge < col_index:
                new_block = merge_two_blocks(block_2, block_1, how)
                block_dict = update_block_dict_with_same_id(block_dict, new_block, block_2["id"], block_1["id"])
                block_id_matrix = update_block_matrix_with_same_id(block_id_matrix, new_block["id"], block_2["id"], block_1["id"])
                del block_dict[block_1["id"]]
            else: 
                new_block = merge_two_blocks(block_1, block_2, how)
                block_dict = update_block_dict_with_same_id(block_dict, new_block, block_1["id"], block_2["id"])
                block_id_matrix = update_block_matrix_with_same_id(block_id_matrix, new_block["id"], block_1["id"], block_2["id"])
                del block_dict[block_2["id"]]
            
    return block_dict, block_id_matrix, "lsr2", str(seed)
 


#----------------------------------------------------------------------------------------

def main():
    # Reading the MSA
    filename = "test3"  # Change this to your .fa file
    sequences = read_fasta("Test_allignments/"+filename+".fa")

    # Convert sequences to matrix
    sequence_matrix = sequences_to_matrix(sequences) 

    # Convert each label in a block
    block_dict, block_id_matrix = build_blocks_from_sequence_matrix(sequence_matrix)
    
    # Euristich
    block_dict, block_id_matrix, euristich, seed = local_search_random_2(block_dict, block_id_matrix)
    #print(block_dict, block_id_matrix)

    #Graph
    graph, pos = build_graph(block_dict, block_id_matrix)

    # Draw the graph 
    nx.draw(graph, pos=pos, labels=nx.get_node_attributes(graph, 'label'), with_labels=True)
    nx.draw_networkx_edge_labels(graph, pos=pos, edge_labels=nx.get_edge_attributes(graph, 'label'), font_color='red')

    # Save the image to a file
    plt.savefig('Graphs_from_test/png_files/graph_image_'+filename+'_'+euristich+'_seed_'+seed+'.png')

    # Show the graph
    plt.show()

    # Save graph
    file_path = 'Graphs_from_test/gfa_files/graph_'+filename+'_'+euristich+'_seed_'+seed+'.gfa'
    graph_to_gfa(graph, file_path)


if __name__ == "__main__":
    main()