# Dependencies
from Bio import SeqIO
import networkx as nx
import matplotlib.pyplot as plt
import random
import numpy as np
import math
from typing import List, Dict, Tuple

#----------------------------------------------------------------------------------------
# Block
def create_block(id: str, label: str, sequence_id: List[str], begin_column: int, end_column: int):
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
def read_fasta(filename: str) -> Dict[str, List[str]]:
    sequences = {}
    for record in SeqIO.parse(filename, "fasta"):
        sequences[record.id] = list(str(record.seq))  # Convert sequence to a list of characters
    return sequences

# Put the sequences in a matrix
def sequences_to_matrix(sequences: Dict[str, List[str]]) -> List[List[str]]:
    max_length = max(len(seq) for seq in sequences.values())
    matrix = []
    for sequence_id, sequence in sequences.items():
        padded_sequence = sequence + [''] * (max_length - len(sequence))  # Pad shorter sequences with empty strings
        matrix.append(padded_sequence)
    return matrix

# From labels to block dict and block id matrix, every single label is a block
def build_blocks_from_sequence_matrix(sequence_matrix: List[List[str]]) -> Tuple[Dict[str, Dict], np.ndarray]:
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

# From labels to block dict and block id matrix, every sequence is a block
def build_long_blocks_from_sequence_matrix(sequence_matrix: List[List[str]]) -> Tuple[Dict[str, Dict], np.ndarray]:
    block_dict = {}
    num_seq = len(sequence_matrix)
    num_bases = len(sequence_matrix[0])    
    block_id_matrix = np.empty((num_seq, num_bases), dtype=object)

    for sequence_index in range(num_seq):
        block_id = "B" + str(sequence_index)
        sequence = ''.join(sequence_matrix[sequence_index])
        for base_index in range(num_bases):  
            block_dict[block_id] = create_block(
                block_id,
                sequence,
                sequence_index,
                0,
                num_bases-1
            )
            block_id_matrix[sequence_index, base_index] = block_id
    
    return block_dict, block_id_matrix

# From dict block and matrix, do greedy row block merge
def greedy_row_merge(block_dict: Dict, block_id_matrix: np.ndarray) -> Tuple[Dict[str, Dict], np.ndarray]:
    num_seq, num_bases = block_id_matrix.shape
    for i in range(num_seq):
        for j in range(num_bases):
            for z in range(num_seq):
                block_1 = block_dict[block_id_matrix[i, j]] 
                block_2 = block_dict[block_id_matrix[z, j]]
                if block_1["id"] != block_2["id"]:
                    if block_1["label"] == block_2["label"]:
                        if z < i:
                            new_block = merge_two_blocks(block_2, block_1, "row_union")
                            block_dict = update_block_dict_with_same_id(block_dict, new_block, block_2["id"], block_1["id"])
                            del block_dict[block_1["id"]]
                        else: 
                            new_block = merge_two_blocks(block_1, block_2, "row_union")
                            block_dict = update_block_dict_with_same_id(block_dict, new_block, block_1["id"], block_2["id"])
                            del block_dict[block_2["id"]]
                        for sequence in new_block["sequence_ids"]:
                            for column in range(new_block["begin_column"], new_block["end_column"] + 1):
                                block_id_matrix[sequence, column] = new_block["id"]

    return block_dict, block_id_matrix

# Merge two blocks in one
def merge_two_blocks(block_1: Dict, block_2: Dict, how: str) -> Dict:
    valid_values = {"column_union", "row_union"}
    if how not in valid_values:
        raise ValueError(f"Invalid value for param: {how}, the accepted values are column_union and row_union")
    
    if how == "column_union":
        new_labels = block_1["label"] + block_2["label"]
        new_block = {
            "id": block_1["id"],
            "label": new_labels,
            "sequence_ids": block_1["sequence_ids"],
            "begin_column": block_1["begin_column"],
            "end_column": block_2["end_column"]
        }
    elif how == "row_union":
        new_sequence_ids = block_1["sequence_ids"] + block_2["sequence_ids"]
        new_block = {
            "id": block_1["id"],
            "label": block_1["label"],
            "sequence_ids": new_sequence_ids,
            "begin_column": block_1["begin_column"],
            "end_column": block_1["end_column"]
        }
    return new_block

# Split one block in two by column
def split_block_by_column(block: Dict) -> Tuple[Dict, Dict]:
    label = block["label"]
    begin_column = block["begin_column"]
    end_column = block["end_column"]
    
    if end_column - begin_column <= 0:
        raise ValueError("Block must span at least two columns to be split.")
    
    # Choose a random split point ensuring both parts are non-empty
    split_point = random.randint(begin_column + 1, end_column)
    
    part1 = {
        "id": block["id"] + "_1",
        "label": label[:split_point - begin_column],
        "sequence_ids": block["sequence_ids"],
        "begin_column": begin_column,
        "end_column": split_point - 1
    }
    
    part2 = {
        "id": block["id"] + "_2",
        "label": label[split_point - begin_column:],
        "sequence_ids": block["sequence_ids"],
        "begin_column": split_point,
        "end_column": end_column
    }
    
    return part1, part2

# Split one block in two by row
def split_block_by_row(block: Dict) -> Tuple[Dict, Dict]:
    sequence_ids = block["sequence_ids"]
    
    if len(sequence_ids) <= 1:
        raise ValueError("Block must contain at least two sequence IDs to be split.")
    
    # Choose a random split point ensuring both parts are non-empty
    split_point = random.randint(1, len(sequence_ids) - 1)
    
    part1 = {
        "id": block["id"] + "_1",
        "label": block["label"],
        "sequence_ids": sequence_ids[:split_point],
        "begin_column": block["begin_column"],
        "end_column": block["end_column"]
    }
    
    part2 = {
        "id": block["id"] + "_2",
        "label": block["label"],
        "sequence_ids": sequence_ids[split_point:],
        "begin_column": block["begin_column"],
        "end_column": block["end_column"]
    }
    
    return part1, part2

# Update block dict data
def update_block_dict_with_same_id(block_dict: Dict, new_block: Dict, id1: str, id2: str) -> Dict:
    for block_id, block_data in block_dict.items():
        if block_id == id1 or block_id == id2:
            block_dict[block_id] = new_block
    return block_dict

# Update block matrix data
def update_block_matrix_with_same_id(block_id_matrix: np.ndarray, new_id: str, id1: str, id2: str) -> np.ndarray:
    rows = len(block_id_matrix)
    cols = len(block_id_matrix[0])
    for i in range(rows):
        for j in range(cols):  
            if block_id_matrix[i, j] == id1 or block_id_matrix[i, j] == id2:
                block_id_matrix[i, j] = new_id
 
    return block_id_matrix

def update_block_submatrix_with_same_id(block_id_matrix: np.ndarray, new_id: str, sequences: List[int], first_column: int, last_column: int) -> np.ndarray:
    for sequence in sequences:
        for column in range(first_column, last_column + 1):
            block_id_matrix[sequence, column] = new_id

# From graph to gfa file
def graph_to_gfa(graph: nx.DiGraph, filename: str) -> None:
    with open(filename, 'w') as f:
        f.write("H\tVN:Z:1.0\n")  # GFA header line
        for node, data in graph.nodes(data=True):
            f.write(f"S\t{node}\t*\tLN:i:{data['label']}\n")  # S-line for each node
        for u, v, data in graph.edges(data=True):
            edge_label = data.get('label', '*')  # Get the edge label, or use '*' if not present
            f.write(f"L\t{u}\t+\t{v}\t+\t{edge_label}\n")  # L-line for each edge with label

# Generate count numbers from seed
def generate_random_numbers(seed: int, start: int, end: int, count: int) -> List[int]:
    random.seed(seed)
    return [random.randint(start, end) for _ in range(count)]

# Generate count number from seed
def generate_random_number(seed: int, start: int, end: int) -> int:
    random.seed(seed)
    return random.randint(start, end)

# Build graph from block dict and block id matrix
def build_graph(block_dict: Dict, block_id_matrix: np.ndarray) -> nx.DiGraph:
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

def acceptance_probability(delta: float, temperature: float) -> float:

    boltzmann_constant = 1.380649e-23  # Boltzmann constant (J/K)

    return math.exp(-delta / (boltzmann_constant * temperature))

def remove_indle_from_graph(graph: nx.DiGraph) -> nx.DiGraph:
    nodes_to_remove = [node for node, data in graph.nodes(data=True) if set(data.get('label', '')) == {"-"}]
    new_edges = []
    for node in nodes_to_remove:
        predecessors = list(graph.predecessors(node))
        successors = list(graph.successors(node))
        for pred in predecessors:
            for succ in successors:
                if graph.has_edge(node, succ):
                    edge_data = graph.get_edge_data(pred, node)
                    graph.add_edge(pred, succ, **edge_data)
     
    graph.remove_nodes_from(nodes_to_remove)
 
    for node, data in graph.nodes(data=True):
        if 'label' in data: 
            cleaned_label = data['label'].replace("-", "")
            if cleaned_label != data['label']:
                data['label'] = cleaned_label
    
    nodes_without_data = [node for node in graph.nodes() if 'label' not in graph.nodes[node]]
    graph.remove_nodes_from(nodes_without_data)

    return graph
#----------------------------------------------------------------------------------------
# Objective functions

def of_min_label_length_threshold(threshold: int, penalization: int, label_length: int) -> int:
    if label_length < threshold :
        return penalization * (threshold - label_length)
    else :
        return 1

def of_pangeblocks(threshold: int, penalization: int, label_length: int) -> int:
    if label_length < threshold :
        return penalization
    else :
        return 1    
