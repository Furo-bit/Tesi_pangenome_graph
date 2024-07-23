import os
import csv

configfile: "config.yaml"


rule all:
    input:
        expand("output/local_search/{alignment}.gfa", alignment=glob_wildcards("data/alignments/{alignment}.fa").alignment),
        expand("output/local_search/{alignment}_quality.txt", alignment=glob_wildcards("data/alignments/{alignment}.fa").alignment),
        expand("output/simulated_annealing/{alignment}.gfa", alignment=glob_wildcards("data/alignments/{alignment}.fa").alignment),
        expand("output/simulated_annealing/{alignment}_quality.txt", alignment=glob_wildcards("data/alignments/{alignment}.fa").alignment),
        expand("output/tabu_search/{alignment}.gfa", alignment=glob_wildcards("data/alignments/{alignment}.fa").alignment),
        expand("output/tabu_search/{alignment}_quality.txt", alignment=glob_wildcards("data/alignments/{alignment}.fa").alignment),
        "output/all_quality.csv",
        "output/weight_scores_plot.png",
        "output/elapsed_time.png",
        "output/max_resident_size.png",
        "output/vertex_plot.png"


rule run_alignment_local_search:
    input:
        alignment="data/alignments/{alignment}.fa",
        params=config["parameters"]["local_search.py"]
    output:
        gfa="output/local_search/{alignment}.gfa",
        quality="output/local_search/{alignment}_quality.txt",
        monitor_log="logs/local_search_{alignment}.log"
    params:
        script="local_search.py"
    shell:
        """
        /usr/bin/time --verbose python {params.script} --params {input.params} --input {input.alignment} --output {output.gfa} --quality {output.quality} > {output.monitor_log} 2>&1
        """

rule run_alignment_simulated_annealing:
    input:
        alignment="data/alignments/{alignment}.fa",
        params=config["parameters"]["simulated_annealing.py"]
    output:
        gfa="output/simulated_annealing/{alignment}.gfa",
        quality="output/simulated_annealing/{alignment}_quality.txt",
        monitor_log="logs/simulated_annealing_{alignment}.log"
    params:
        script="simulated_annealing.py"
    shell:
        """
        /usr/bin/time --verbose python {params.script} --params {input.params} --input {input.alignment} --output {output.gfa} --quality {output.quality} > {output.monitor_log} 2>&1
        """
rule run_alignment_tabu_search:
    input:
        alignment="data/alignments/{alignment}.fa",
        params=config["parameters"]["tabu_search.py"]
    output:
        gfa="output/tabu_search/{alignment}.gfa",
        quality="output/tabu_search/{alignment}_quality.txt",
        monitor_log="logs/tabu_search_{alignment}.log"
    params:
        script="tabu_search.py"
    shell:
        """
        /usr/bin/time --verbose python {params.script} --params {input.params} --input {input.alignment} --output {output.gfa} --quality {output.quality} > {output.monitor_log} 2>&1
        """


rule concatenate_quality:
    input:
        local_search_quality_files=expand("output/local_search/{alignment}_quality.txt", alignment=glob_wildcards("data/alignments/{alignment}.fa").alignment),
        simulated_annealing_quality_files=expand("output/simulated_annealing/{alignment}_quality.txt", alignment=glob_wildcards("data/alignments/{alignment}.fa").alignment),
        tabu_search_quality_files=expand("output/tabu_search/{alignment}_quality.txt", alignment=glob_wildcards("data/alignments/{alignment}.fa").alignment)
    output:
        "output/all_quality.csv"
    run:
        with open(output[0], 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(['Alignment', 'Quality', 'Program'])
            for quality_file in input.local_search_quality_files:
                alignment_name = os.path.basename(quality_file).split("_quality.txt")[0]
                with open(quality_file, 'r') as infile:
                    for line in infile:
                        csvwriter.writerow([alignment_name, line.strip(), "local_search.py"])
            for quality_file in input.simulated_annealing_quality_files:
                alignment_name = os.path.basename(quality_file).split("_quality.txt")[0]
                with open(quality_file, 'r') as infile:
                    for line in infile:
                        csvwriter.writerow([alignment_name, line.strip(), "simulated_annealing.py"])
            for quality_file in input.tabu_search_quality_files:
                alignment_name = os.path.basename(quality_file).split("_quality.txt")[0]
                with open(quality_file, 'r') as infile:
                    for line in infile:
                        csvwriter.writerow([alignment_name, line.strip(), "tabu_search.py"])

rule plot_quality:
    input:
        "output/all_quality.csv"
    output:
        "output/weight_scores_plot.png"
    shell:
        """
        python plot_quality.py
        """

rule generate_plots:
    input:
        "output/weight_scores_plot.png",
        logs="logs"
    output:
        elapsed_time="output/elapsed_time.png",
        max_resident_size="output/max_resident_size.png"
    shell:
        "python log_plot.py {input.logs} output"

rule plot_vertices:
    input:
        "output/elapsed_time.png",
        directories=["output/local_search", "output/simulated_annealing", "output/tabu_search"]
    output:
        vertex_plot="output/vertex_plot.png"
    shell:
        "python vertex_plot.py {input.directories} output"