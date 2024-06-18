import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the CSV file
df = pd.read_csv('output/all_quality.csv')

# Convert the 'Quality' column to numeric
df['Quality'] = pd.to_numeric(df['Quality'])

# Set the style of the plot
plt.style.use('ggplot')

# Create the figure and the axis
fig, ax = plt.subplots(figsize=(10, 6))

# Define colors for each program
colors = {'local_search.py': 'blue', 'simulated_annealing.py': 'green'}

# Define the width of a single bar
bar_width = 0.4

# Get unique alignments and the number of programs
alignments = df['Alignment'].unique()
num_programs = len(df['Program'].unique())

# Create a dictionary to store bar positions for each program
bar_positions = {}
for i, program in enumerate(df['Program'].unique()):
    bar_positions[program] = np.arange(len(alignments)) + i * bar_width

# Plot each program's bars
for program in df['Program'].unique():
    subset = df[df['Program'] == program]
    ax.bar(bar_positions[program], subset['Quality'], label=program, width=bar_width, color=colors[program])

# Set the position of x ticks and their labels
ax.set_xticks(np.arange(len(alignments)) + bar_width / (num_programs - 1))
ax.set_xticklabels(alignments)

# Add labels and title
ax.set_xlabel('Alignment')
ax.set_ylabel('Weight')
ax.set_title('Weight scores by Alignment and Heuristic')

# Add the legend
ax.legend(title='Heuristic')

# Rotate x-axis labels for better readability
plt.xticks(rotation=90)

# Adjust layout to make room for the rotated x-axis labels
plt.tight_layout()

# Save the plot as an image
plt.savefig('output/weight_scores_plot.png')

# Show the plot
#plt.show()
