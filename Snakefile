import os
import csv

# Specifica quale script utilizzare
configfile: "config.yaml"

rule all:
    input:
        expand("output/{alignment}.gfa", alignment=glob_wildcards("data/alignments/{alignment}.fa").alignment),
        expand("output/{alignment}_quality.txt", alignment=glob_wildcards("data/alignments/{alignment}.fa").alignment),
        "output/all_quality.csv"  # Cambiato in .csv

rule run_alignment:
    input:
        alignment="data/alignments/{alignment}.fa",
        params=lambda wildcards: config["parameters"][config["script"]]  # File di parametri basato sullo script selezionato
    output:
        gfa="output/{alignment}.gfa",
        quality="output/{alignment}_quality.txt"
    params:
        alignment_name=glob_wildcards("data/alignments/{alignment}.fa").alignment,
        script=config["script"]  # Utilizza lo script specificato nel file di configurazione
    shell:
        """
        python {params.script} --params {input.params} --input {input.alignment} --output {output.gfa} --quality {output.quality}
        """

rule concatenate_quality:
    input:
        quality_files=expand("output/{alignment}_quality.txt", alignment=glob_wildcards("data/alignments/{alignment}.fa").alignment)
    output:
        "output/all_quality.csv"  # Cambiato in .csv
    params:
        script=config["script"]  # Nome dello script dal file di configurazione
    run:
        with open(output[0], 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(['Alignment', 'Quality', 'Program'])  # Scrive l'intestazione
            for quality_file in input.quality_files:
                alignment_name = os.path.basename(quality_file).split("_quality.txt")[0]
                with open(quality_file, 'r') as infile:
                    for line in infile:
                        csvwriter.writerow([alignment_name, line.strip(), params.script])
