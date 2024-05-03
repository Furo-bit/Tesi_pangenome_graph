# Pangenome graph 

# Dependencies
from Bio import SeqIO
import networkx as nx
import matplotlib.pyplot as plt
import random

#----------------------------------------------------------------------------------------
# Block class
class block:
    def __init__(self, id, base, sequence_id, begin_column, end_column):
        self.id = id
        self.base = base
        if isinstance(sequence_id, list):
            self.sequences_ids = sequence_id
        else:
            self.sequences_ids = [sequence_id]
        self.begin_column = begin_column
        self.end_column = end_column
#----------------------------------------------------------------------------------------
# Functions

# read msa file
def read_fasta(filename):
    sequences = {}
    for record in SeqIO.parse(filename, "fasta"):
        sequences[record.id] = list(str(record.seq))  # Convert sequence to a list of characters
    return sequences
     
# put the sequences in a matrix
def sequences_to_matrix(sequences):
    max_length = max(len(seq) for seq in sequences.values())
    matrix = []
    for sequence_id, sequence in sequences.items():
        padded_sequence = sequence + [''] * (max_length - len(sequence))  # Pad shorter sequences with empty strings
        matrix.append(padded_sequence)
    return matrix

# convert each base in a block
def build_blocks_from_sequence_matrix(sequence_matrix):
    block_matrix = [row[:] for row in sequence_matrix]
    num_seq = len(sequence_matrix)
    num_bases = len(sequence_matrix[0])
    for sequence_index in range(num_seq):
        for base_index in range(num_bases):
            block_matrix[sequence_index][base_index] = block(
                "B"+str(sequence_index)+"."+str(base_index),
                sequence_matrix[sequence_index][base_index],
                sequence_index,
                base_index,
                base_index 
            )
    return block_matrix

# build pangenome graph from block matrix
def build_graph(matrix):
    G = nx.DiGraph()
    rows = len(matrix)
    cols = len(matrix[0])

    # Add nodes
    G.add_node("source", label = "source")
    G.add_node("sink", label = "sink")
    for i in range(rows):
        for j in range(cols):
            G.add_node(matrix[i][j].id, label=matrix[i][j].base, row=i, col=j)
            
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

# Convert the nx graph into a gfa file
# TO DO: aggiungere percorso per ogni sequenza
def graph_to_gfa(graph, filename):
    with open(filename, 'w') as f:
        f.write("H\tVN:Z:1.0\n")  # GFA header line
        for node, data in graph.nodes(data=True):
            f.write(f"S\t{node}\t*\tLN:i:{data['label']}\n")  # S-line for each node
        for u, v, data in graph.edges(data=True):
            edge_label = data.get('label', '*')  # Get the edge label, or use '*' if not present
            f.write(f"L\t{u}\t+\t{v}\t+\t{edge_label}\n")  # L-line for each edge with label

#Merge two block in one
def merge_two_blocks(block_1, block_2, how):
    if how == "column_union":
        new_bases = block_1.base + block_2.base 
        new_block = block(
            block_1.id,
            new_bases,
            block_1.sequences_ids,
            block_1.begin_column,
            block_2.end_column
        )
        print("Merged block", block_1.id, "with block", block_2.id, "by columns")
    if how == "row_union":
        new_sequences_ids = block_1.sequences_ids + block_2.sequences_ids
        new_block = block(
            block_1.id,
            block_1.base,
            new_sequences_ids,
            block_1.begin_column,
            block_1.end_column
        )
        print("Merged block", block_1.id, "with block", block_2.id, "by rows")
    return new_block

#Update blocks data
def update_blocks_with_same_id(block_matrix, new_block, id1, id2):
    rows = len(block_matrix)
    cols = len(block_matrix[0])
    for i in range(rows):
        for j in range(cols):  
            if block_matrix[i][j].id == id1 or block_matrix[i][j].id == id2:
                block_matrix[i][j] = new_block
    return block_matrix  

# Random number from range and seed
def generate_random_number(seed, start, end):
    random.seed(seed)
    return random.randint(start, end)


def generate_random_numbers(seed, start, end, count):
    random.seed(seed)
    return [random.randint(start, end) for _ in range(count)]


# Euristichs
def local_search(block_matrix):
    rows = len(block_matrix)
    cols = len(block_matrix[0])
    for i in range(rows):
        for j in range(cols):
            # Try to merge two blocks by rows
            if i+1 != rows:
                if (block_matrix[i][j].begin_column == block_matrix[i+1][j].begin_column and
                     block_matrix[i][j].end_column == block_matrix[i+1][j].end_column and
                       block_matrix[i][j].base == block_matrix[i+1][j].base):
                    new_block = merge_two_blocks(block_matrix[i][j], block_matrix[i+1][j], "row_union")
                    block_matrix = update_blocks_with_same_id(block_matrix, new_block, block_matrix[i][j].id, block_matrix[i+1][j].id)
            # Try to merge two blocks by columns
            if j+1 != cols:
                if block_matrix[i][j].sequences_ids == block_matrix[i][j+1].sequences_ids:
                    new_block = merge_two_blocks(block_matrix[i][j], block_matrix[i][j+1], "column_union")
                    block_matrix = update_blocks_with_same_id(block_matrix, new_block, block_matrix[i][j].id, block_matrix[i][j+1].id)
    return block_matrix, "ls"

def local_search_random(block_matrix):
    
    # Ask the user for a seed
    seed = int(input("Insert the seed for the random number: "))

    rows = len(block_matrix)
    cols = len(block_matrix[0])
    cell_total = rows * cols
    # Generate a random numbers using the user-provided seed
    random_numbers = generate_random_numbers(seed, 1, cell_total, cell_total*2)
    # Create a boolean matrix 
    boolean_matrix = [[False] * cols for _ in range(rows)]

    for x in random_numbers:
        # Transform the number in matrix indeces 
        row_index = (x - 1) // cols
        col_index = (x - 1) % cols
        # Check if the block is visited, if not then visit and try to merge
        if boolean_matrix[row_index][col_index] == False :
            boolean_matrix[row_index][col_index] = True
            # Check up block
            if row_index-1 >= 0:
                if block_matrix[row_index][col_index].id != block_matrix[row_index-1][col_index].id:
                    if (block_matrix[row_index][col_index].begin_column == block_matrix[row_index-1][col_index].begin_column and
                            block_matrix[row_index][col_index].end_column == block_matrix[row_index-1][col_index].end_column and
                            block_matrix[row_index][col_index].base == block_matrix[row_index-1][col_index].base):
                        new_block = merge_two_blocks(block_matrix[row_index-1][col_index], block_matrix[row_index][col_index], "row_union")
                        block_matrix = update_blocks_with_same_id(block_matrix, new_block, block_matrix[row_index][col_index].id, block_matrix[row_index-1][col_index].id)
            # Check down block       
            if row_index+1 < rows:
                if block_matrix[row_index][col_index].id != block_matrix[row_index+1][col_index].id:
                    if (block_matrix[row_index][col_index].begin_column == block_matrix[row_index+1][col_index].begin_column and
                            block_matrix[row_index][col_index].end_column == block_matrix[row_index+1][col_index].end_column and
                            block_matrix[row_index][col_index].base == block_matrix[row_index+1][col_index].base):
                        new_block = merge_two_blocks(block_matrix[row_index][col_index], block_matrix[row_index+1][col_index], "row_union")
                        block_matrix = update_blocks_with_same_id(block_matrix, new_block, block_matrix[row_index][col_index].id, block_matrix[row_index+1][col_index].id)
            # Check left block
            if col_index-1 >= 0:
                if block_matrix[row_index][col_index].id != block_matrix[row_index][col_index-1].id:
                    if block_matrix[row_index][col_index].sequences_ids == block_matrix[row_index][col_index-1].sequences_ids:
                        new_block = merge_two_blocks(block_matrix[row_index][col_index-1], block_matrix[row_index][col_index], "column_union")
                        block_matrix = update_blocks_with_same_id(block_matrix, new_block, block_matrix[row_index][col_index].id, block_matrix[row_index][col_index-1].id)
            # Check right block
            if col_index+1 < cols:
                if block_matrix[row_index][col_index].id != block_matrix[row_index][col_index+1].id:
                    if block_matrix[row_index][col_index].sequences_ids == block_matrix[row_index][col_index+1].sequences_ids:
                        new_block = merge_two_blocks(block_matrix[row_index][col_index], block_matrix[row_index][col_index+1], "column_union")
                        block_matrix = update_blocks_with_same_id(block_matrix, new_block, block_matrix[row_index][col_index].id, block_matrix[row_index][col_index+1].id)

    return block_matrix, "lsr", str(seed)



#----------------------------------------------------------------------------------------
# Main

# Read msa from file
filename = "test4"  # Change this to your .fa file
sequences = read_fasta("Test_allignments/"+filename+".fa")

# Convert sequences to matrix
sequence_matrix = sequences_to_matrix(sequences) 

# Convert each base in a block
block_matrix = build_blocks_from_sequence_matrix(sequence_matrix)

# Euristich
block_matrix, euristich, seed = local_search_random(block_matrix)

# Create the pangenome graph using the blocks
# Each block is a node
graph, pos = build_graph(block_matrix)

# Draw the graph 
nx.draw(graph, pos=pos, labels=nx.get_node_attributes(graph, 'label'), with_labels=True)
nx.draw_networkx_edge_labels(graph, pos=pos, edge_labels=nx.get_edge_attributes(graph, 'label'), font_color='red')

# Save the image to a file
plt.savefig('Graphs_from_test/png_files/graph_image_'+filename+'_'+euristich+'_'+seed+'.png')

# Show the graph
plt.show()

#Save graph
file_path = 'Graphs_from_test/gfa_files/graph_'+filename+'_'+euristich+'_'+seed+'.gfa'
graph_to_gfa(graph, file_path)

