import os
import sys
import matplotlib.pyplot as plt

def count_vertices_in_gfa(file_path):
    """Count the number of vertices (segments) in a GFA file."""
    vertex_count = 0
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('S'):
                vertex_count += 1
    return vertex_count

def plot_vertices_count(directories, output_path):
    """Plot the number of vertices for each GFA file in the given directories and save the plot."""
    vertex_counts = []

    # Process each directory
    for directory in directories:
        if not os.path.isdir(directory):
            print(f"{directory} is not a valid directory.")
            continue

        # List and sort only files ending with .gfa
        files = sorted(f for f in os.listdir(directory) if f.endswith('.gfa'))

        for filename in files:
            file_path = os.path.join(directory, filename)
            try:
                num_vertices = count_vertices_in_gfa(file_path)
                # Append the result with directory name as prefix
                vertex_counts.append((f"{os.path.basename(directory)}/{filename}", num_vertices))
            except Exception as e:
                print(f"Error processing {file_path}: {e}")

    if not vertex_counts:
        print("No GFA files found or all are empty in the provided directories.")
        return

    # Separate file names and vertex counts
    file_names, vertex_counts = zip(*vertex_counts)

    # Plotting
    plt.figure(figsize=(12, 6))
    plt.bar(file_names, vertex_counts, color='skyblue')
    plt.xlabel('File Name (Directory/File)')
    plt.ylabel('Number of Vertices')
    plt.title('Number of Vertices in GFA Files from Multiple Directories')
    plt.xticks(rotation=90)
    plt.tight_layout()

    # Save the plot instead of showing it
    #plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_path, 'vertex_plot.png'))
    plt.close()

if __name__ == "__main__":
    # Read directories and output file path from command line arguments
    directories = sys.argv[1:-1]  # All but the last argument
    output_file = sys.argv[-1]    # The last argument
    plot_vertices_count(directories, output_file)
