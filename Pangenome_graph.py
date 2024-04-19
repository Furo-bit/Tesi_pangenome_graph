# Pangenome graph 

# Dependencies
from Bio import SeqIO
import networkx as nx
import matplotlib.pyplot as plt

#----------------------------------------------------------------------------------------
# Block class
class block:
    def __init__(self, base, sequence_id, begin_column, end_column):
        self.base = base
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
            ij_string = str(i) + str(j)
            G.add_node(ij_string, label=matrix[i][j].base, row=i, col=j)

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
        G.add_edge("source", str(i)+str(0), label = str(i))
        G.add_edge(str(i)+str(cols-1), "sink", label = str(i))
        for j in range(cols):
            ij_string = str(i) + str(j)
            ij1_string = str(i) + str(j+1)
            # Connect horizontally adjacent nodes
            if j < cols - 1:
                G.add_edge(ij_string, ij1_string, label = str(i))
    
    return G, pos

# Convert the nx graph into a gfa file
# Da rivedere
def graph_to_gfa(graph, filename):
    with open(filename, 'w') as f:
        f.write("H\tVN:Z:1.0\n")  # GFA header line
        for node, data in graph.nodes(data=True):
            f.write(f"S\t{node}\t*\tLN:i:{data['label']}\n")  # S-line for each node
        for u, v, data in graph.edges(data=True):
            f.write(f"L\t{u}\t+\t{v}\t+\t*\n")  # L-line for each edge

#----------------------------------------------------------------------------------------
# Main
# Read msa from file
filename = "Test_allignments/test.fa"  # Change this to your .fa file path
sequences = read_fasta(filename)

# Convert sequences to matrix
sequence_matrix = sequences_to_matrix(sequences) 
# Convert each base in a block
block_matrix = build_blocks_from_sequence_matrix(sequence_matrix)

# Create the pangenome graph using the blocks
# Each node is a block
graph, pos = build_graph(block_matrix)

# Draw the graph 
nx.draw(graph, pos=pos, labels=nx.get_node_attributes(graph, 'label'), with_labels=True)
nx.draw_networkx_edge_labels(graph, pos=pos, edge_labels=nx.get_edge_attributes(graph, 'label'), font_color='red')

# Save the image to a file
plt.savefig('Graphs_from_test/graph_image.png')

# Show the graph
plt.show()

#Save graph
file_path = "Graphs_from_test/graph.gfa"
graph_to_gfa(graph, file_path)

