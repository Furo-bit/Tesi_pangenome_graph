import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file
df = pd.read_csv('output/all_quality.csv')

# Convert the 'Quality' column to numeric
df['Quality'] = pd.to_numeric(df['Quality'])

# Set the style of the plot
plt.style.use('ggplot')

# Create the figure and the axis
fig, ax = plt.subplots(figsize=(10, 6))

# Create the bar plot
colors = {'local_search.py': 'blue', 'simulated_annealing.py': 'green'}
df['Color'] = df['Program'].map(colors)
for program in df['Program'].unique():
    subset = df[df['Program'] == program]
    ax.bar(subset['Alignment'], subset['Quality'], label=program, color=subset['Color'])

# Add labels and title
ax.set_xlabel('Alignment')
ax.set_ylabel('Quality')
ax.set_title('Quality Scores by Alignment and Heuristic')

# Add the legend
ax.legend(title='Heuristic')

# Rotate x-axis labels for better readability
plt.xticks(rotation=90)

# Adjust layout to make room for the rotated x-axis labels
plt.tight_layout()

# Save the plot as an image
plt.savefig('output/quality_scores_plot.png')

# Show the plot
plt.show()
