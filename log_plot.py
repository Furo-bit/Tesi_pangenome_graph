import os
import sys
import pandas as pd
import matplotlib.pyplot as plt

def parse_log_file(filepath):
    data = {}
    with open(filepath, 'r') as file:
        for line in file:
            if "	User time (seconds):" in line:
                # Convert seconds to minutes
                user_time = float(line.split(":")[1].strip()) / 60
                data['User Time'] = user_time
            elif "	System time (seconds):" in line:
                # Convert seconds to minutes
                system_time = float(line.split(":")[1].strip()) / 60
                data['System Time'] = system_time
            elif "	Maximum resident set size (kbytes):" in line:
                # Convert kilobytes to gigabytes
                max_resident_size = int(line.split(":")[1].strip()) / (1024 * 1024)
                data['Max Resident Size'] = max_resident_size

    if 'User Time' in data and 'System Time' in data:
        data['Elapsed Time'] = data['User Time'] + data['System Time']
    return data

log_dir = sys.argv[1]
output_dir = sys.argv[2]

os.makedirs(output_dir, exist_ok=True)

data_list = []
file_names = []

for filename in sorted(os.listdir(log_dir)):
    if filename.endswith('.log'):
        file_path = os.path.join(log_dir, filename)
        parsed_data = parse_log_file(file_path)
        if parsed_data:
            data_list.append(parsed_data)
            file_names.append(filename)
        else:
            print(f"No data parsed from {file_path}")

df = pd.DataFrame(data_list, index=file_names)

if 'Elapsed Time' in df.columns and 'Max Resident Size' in df.columns:
    # Plotting Elapsed Time
    plt.figure(figsize=(14, 8))
    plt.bar(df.index, df['Elapsed Time'], color='tab:blue', label='Elapsed Time (minutes)')
    plt.xlabel('Log File')
    plt.ylabel('Elapsed Time (minutes)')
    plt.title('Elapsed Time from Log Files')
    plt.xticks(rotation=90)  
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'elapsed_time.png'))
    plt.close()

    # Plotting Max Resident Size
    plt.figure(figsize=(14, 8))
    plt.bar(df.index, df['Max Resident Size'], color='tab:red', label='Max Resident Size (GB)')
    plt.xlabel('Log File')
    plt.ylabel('Max Resident Size (GB)')
    plt.title('Max Resident Size from Log Files')
    plt.xticks(rotation=90)  
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'max_resident_size.png'))
    plt.close()
else:
    print("Errore: Le colonne 'Elapsed Time' e/o 'Max Resident Size' non sono presenti nel DataFrame")
