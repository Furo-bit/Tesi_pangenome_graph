#include "utils.h"

void simulated_annealing(const unordered_map<string, Block>& block_dict, const vector<vector<string>>& block_id_matrix){


    }



int main() {

    // Leggi le sequenze dal file FASTA
    const string filename = "test.fa";
    auto sequences = read_fasta(filename);

    // Converti le sequenze in una matrice
    auto sequence_matrix = sequences_to_matrix(sequences);
    

    // Costruisci i blocchi e la matrice degli ID dei blocchi
    auto [block_dict, block_id_matrix] = build_blocks_from_sequence_matrix(sequence_matrix);

    // Inserire qua euristica
    simulated_annealing(block_dict, block_id_matrix);

    // Costruisce il grafo e lo salva 
    Agraph_t* g = build_graph(block_dict, block_id_matrix);

    return 0;
}