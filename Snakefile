import os

rule all:
    input:
        expand("output/{alignment}.gfa", alignment=glob_wildcards("data/alignments/{alignment}.fa").alignment),
        expand("output/{alignment}_quality.txt", alignment=glob_wildcards("data/alignments/{alignment}.fa").alignment),
        "output/all_quality.txt"

rule run_alignment:
    input:
        params="data/params.txt",
        alignment="data/alignments/{alignment}.fa"
    output:
        gfa="output/{alignment}.gfa",
        quality="output/{alignment}_quality.txt"
    params:
        alignment_name=glob_wildcards("data/alignments/{alignment}.fa").alignment
    shell:
        """
        python Pangenome_graph_dict.py --params {input.params} --input {input.alignment} --output {output.gfa} --quality {output.quality}
        """

rule concatenate_quality:
    input:
        quality_files=expand("output/{alignment}_quality.txt", alignment=glob_wildcards("data/alignments/{alignment}.fa").alignment)
    output:
        "output/all_quality.txt"
    run:
        with open(output[0], 'w') as outfile:
            for quality_file in input.quality_files:
                alignment_name = os.path.basename(quality_file).split("_quality.txt")[0]
                with open(quality_file, 'r') as infile:
                    for line in infile:
                        outfile.write(f"{alignment_name}\t{line.strip()}\n")
