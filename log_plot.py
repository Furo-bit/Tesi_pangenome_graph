import os
import sys
import pandas as pd
import matplotlib.pyplot as plt

# Funzione per estrarre i dati da un file di log
def parse_log_file(filepath):
    data = {}
    with open(filepath, 'r') as file:
        for line in file:
            if "	User time (seconds):" in line:
                user_time = float(line.split(":")[1].strip())
                data['User Time'] = user_time
            elif "	System time (seconds):" in line:
                system_time = float(line.split(":")[1].strip())
                data['System Time'] = system_time
            elif "	Maximum resident set size (kbytes):" in line:
                max_resident_size = int(line.split(":")[1].strip())
                data['Max Resident Size'] = max_resident_size

    if 'User Time' in data and 'System Time' in data:
        data['Elapsed Time'] = data['User Time'] + data['System Time']
    return data

# Leggi i percorsi delle directory dagli argomenti della riga di comando
log_dir = sys.argv[1]
output_dir = sys.argv[2]

# Crea la cartella di output se non esiste
os.makedirs(output_dir, exist_ok=True)

# Lista per conservare i dati estratti
data_list = []
file_names = []

# Lettura dei file di log in ordine alfabetico
for filename in sorted(os.listdir(log_dir)):
    if filename.endswith('.log'):
        file_path = os.path.join(log_dir, filename)
        parsed_data = parse_log_file(file_path)
        if parsed_data:
            data_list.append(parsed_data)
            file_names.append(filename)
        else:
            print(f"No data parsed from {file_path}")

# Creazione di un DataFrame
df = pd.DataFrame(data_list, index=file_names)

# Verifica che le colonne esistano nel DataFrame
if 'Elapsed Time' in df.columns and 'Max Resident Size' in df.columns:
    # Creazione del grafico a barre per il tempo
    plt.figure(figsize=(14, 8))
    plt.bar(df.index, df['Elapsed Time'], color='tab:blue', label='Elapsed Time (s)')
    plt.xlabel('Log File')
    plt.ylabel('Elapsed Time (s)')
    plt.title('Elapsed Time from Log Files')
    plt.xticks(rotation=90)  # Ruota le etichette dell'asse x per leggibilità
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig('output/elapsed_time.png')
    plt.close()

    # Creazione del grafico a barre per la memoria
    plt.figure(figsize=(14, 8))
    plt.bar(df.index, df['Max Resident Size'], color='tab:red', label='Max Resident Size (KB)')
    plt.xlabel('Log File')
    plt.ylabel('Max Resident Size (KB)')
    plt.title('Max Resident Size from Log Files')
    plt.xticks(rotation=90)  # Ruota le etichette dell'asse x per leggibilità
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig('output/max_resident_size.png')
    plt.close()
else:
    print("Errore: Le colonne 'Elapsed Time' e/o 'Max Resident Size' non sono presenti nel DataFrame")
