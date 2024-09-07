#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <iomanip>
#include <algorithm>
#include <stdexcept>

using namespace std;

// Definizione della struttura Block
struct Block {
    string id;
    string label;
    vector<int> sequence_ids;
    int begin_column;
    int end_column;
};

// Funzione per creare un Block
Block create_block(const string& id, const string& label, 
                   const vector<int>& sequence_ids, int begin_column, int end_column) {
    Block block;
    block.id = id;
    block.label = label;
    block.sequence_ids = sequence_ids;
    block.begin_column = begin_column;
    block.end_column = end_column;
    
    return block;
}

// Funzione per stampare i dettagli di un Block
void print_block(const Block& block) {
    cout << "ID: " << block.id << ", "
         << "Label: " << block.label << ", "
         << "Sequence IDs: ";
    for (int id : block.sequence_ids) {
        cout << id << " ";
    }
    cout << ", Begin Column: " << block.begin_column
         << ", End Column: " << block.end_column
         << endl;
}

// Funzione per leggere un file FASTA e restituire le sequenze
unordered_map<string, vector<char>> read_fasta(const string& filename) {
    unordered_map<string, vector<char>> sequences;
    ifstream file(filename);

    if (!file.is_open()) {
        cerr << "Errore nell'apertura del file: " << filename << endl;
        return sequences;
    }

    string line;
    string current_id;

    while (getline(file, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '>') {
            current_id = line.substr(1);
            sequences[current_id] = vector<char>();
        } else {
            sequences[current_id].insert(sequences[current_id].end(), line.begin(), line.end());
        }
    }

    file.close();
    return sequences;
}

// Funzione per convertire le sequenze in una matrice
vector<vector<char>> sequences_to_matrix(const unordered_map<string, vector<char>>& sequences) {
    vector<vector<char>> matrix;
    
    if (sequences.empty()) {
        return matrix;
    }
    
    size_t length = sequences.begin()->second.size();
    
    for (const auto& pair : sequences) {
        if (pair.second.size() != length) {
            cerr << "Le sequenze non hanno la stessa lunghezza!" << endl;
            return {};
        }
        matrix.push_back(pair.second);
    }
    
    return matrix;
}

// Funzione per costruire blocchi da una matrice di sequenze
pair<unordered_map<string, Block>, vector<vector<string>>> build_blocks_from_sequence_matrix(
    const vector<vector<char>>& sequence_matrix) {
    
    unordered_map<string, Block> block_dict;
    size_t num_seq = sequence_matrix.size();
    size_t num_bases = sequence_matrix.empty() ? 0 : sequence_matrix[0].size() - 1;
    vector<vector<string>> block_id_matrix(num_seq, vector<string>(num_bases, ""));
    
    for (size_t sequence_index = 0; sequence_index < num_seq; ++sequence_index) {
        for (size_t base_index = 0; base_index < num_bases; ++base_index) {
            string block_id = "B" + to_string(sequence_index) + "." + to_string(base_index);
            Block block = create_block(
                block_id,
                string(1, sequence_matrix[sequence_index][base_index]),
                vector<int>{static_cast<int>(sequence_index)},
                static_cast<int>(base_index),
                static_cast<int>(base_index)
            );
            block_dict[block_id] = block;
            block_id_matrix[sequence_index][base_index] = block_id;
        }
    }
    
    return {block_dict, block_id_matrix};
}

// Funzione per ordinare i blocchi in base al sequence_id
vector<Block> sort_blocks_by_sequence_id(const unordered_map<string, Block>& block_dict) {
    vector<Block> blocks;
    
    for (const auto& pair : block_dict) {
        blocks.push_back(pair.second);
    }
    
    sort(blocks.begin(), blocks.end(), [](const Block& a, const Block& b) {
        return a.sequence_ids[0] < b.sequence_ids[0];
    });
    
    return blocks;
}

// Funzione per unire due blocchi
Block merge_two_blocks(const Block& block_1, const Block& block_2, const string& how) {
    const vector<string> valid_values = {"column_union", "row_union"};
    if (find(valid_values.begin(), valid_values.end(), how) == valid_values.end()) {
        throw invalid_argument("Invalid value for param: " + how + ", the accepted values are column_union and row_union");
    }

    Block new_block;

    if (how == "column_union") {
        new_block.id = block_1.id;
        new_block.label = block_1.label + block_2.label;
        new_block.sequence_ids = block_1.sequence_ids; // Assumendo che sequence_ids siano rappresentati come un vettore di int
        new_block.begin_column = block_1.begin_column;
        new_block.end_column = block_2.end_column;
    } else if (how == "row_union") {
        new_block.id = block_1.id;
        new_block.label = block_1.label;
        new_block.sequence_ids = block_1.sequence_ids; // Unione dei sequence_ids
        new_block.sequence_ids.insert(new_block.sequence_ids.end(), block_2.sequence_ids.begin(), block_2.sequence_ids.end());
        new_block.begin_column = block_1.begin_column;
        new_block.end_column = block_1.end_column;
    }

    return new_block;
}

int main() {
    // Leggi le sequenze dal file FASTA
    const string filename = "test.fa";
    auto sequences = read_fasta(filename);

    // Converti le sequenze in una matrice
    auto sequence_matrix = sequences_to_matrix(sequences);

    // Costruisci i blocchi e la matrice degli ID dei blocchi
    auto [block_dict, block_id_matrix] = build_blocks_from_sequence_matrix(sequence_matrix);

    // Ordina i blocchi per sequence_id
    auto sorted_blocks = sort_blocks_by_sequence_id(block_dict);

    // Stampa dei blocchi ordinati
    cout << "Blocks (sorted by Sequence ID):" << endl;
    for (const auto& block : sorted_blocks) {
        print_block(block);
    }
    
    // Stampa della matrice degli ID dei blocchi
    cout << "Block ID Matrix:" << endl;
    for (const auto& row : block_id_matrix) {
        for (const auto& id : row) {
            cout << setw(10) << id; // Format per l'allineamento
        }
        cout << endl;
    }

    // Test della funzione di unione dei blocchi
    if (sorted_blocks.size() > 1) {
        Block block_1 = sorted_blocks[0];
        Block block_2 = sorted_blocks[1];
        string how = "column_union"; // o "row_union"
        
        try {
            Block merged_block = merge_two_blocks(block_1, block_2, how);
            cout << "Merged Block:" << endl;
            print_block(merged_block);
        } catch (const invalid_argument& e) {
            cerr << "Errore: " << e.what() << endl;
        }
    } else {
        cout << "Non ci sono abbastanza blocchi per eseguire l'unione." << endl;
    }

    return 0;
}
