#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <iomanip>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <random>

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

// Funzione aggiornata per eliminare i blocchi id1 e id2 e aggiungere il nuovo blocco
unordered_map<string, Block> update_block_dict_with_same_id(
    unordered_map<string, Block>& block_dict, const Block& new_block, const string& id1, const string& id2) {
    
    // Elimina i blocchi con ID id1 e id2
    block_dict.erase(id1);
    block_dict.erase(id2);
    
    // Aggiungi il nuovo blocco al dizionario
    block_dict[new_block.id] = new_block;
    
    return block_dict;
}

void update_block_submatrix_with_same_id(vector<vector<string>>& block_id_matrix, 
                                         const string& new_id, 
                                         const vector<int>& sequences, 
                                         int first_column, 
                                         int last_column) {
    // Aggiorna gli ID dei blocchi nelle righe specificate e nelle colonne comprese tra first_column e last_column
    for (int sequence : sequences) {
        for (int column = first_column; column <= last_column; ++column) {
            block_id_matrix[sequence][column] = new_id;
        }
    }
}

double acceptance_probability(double delta, double temperature, double beta) {
    return std::exp(-beta * delta / temperature);
}

// Genera una lista di numeri casuali
vector<int> generate_random_numbers(int seed, int start, int end, int count) {
    vector<int> random_numbers;
    mt19937 generator(seed);  // Generatore di numeri casuali con seme
    uniform_int_distribution<int> distribution(start, end);  // Distribuzione uniforme tra start e end

    for (int i = 0; i < count; ++i) {
        random_numbers.push_back(distribution(generator));
    }

    return random_numbers;
}

// Genera un singolo numero casuale, se rifatto con gli stessi parametri dà lo stesso risultato
int generate_random_number(int seed, int start, int end) {
    mt19937 generator(seed);  // Generatore di numeri casuali con seme
    uniform_int_distribution<int> distribution(start, end);  // Distribuzione uniforme tra start e end
    return distribution(generator);
}

// Funzione per dividere un blocco per riga
tuple<Block, Block> split_block_by_row(const Block& block, int split_number) {
    const vector<int>& sequence_ids = block.sequence_ids;

    if (sequence_ids.size() <= 1) {
        throw invalid_argument("Block must contain at least two sequence IDs to be split.");
    }

    // Generatore di numeri casuali per il punto di divisione
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> dist(1, sequence_ids.size() - 1);
    int split_point = dist(gen);

    // Creazione dei nuovi ID dei blocchi
    string id1, id2;
    size_t underscore_index = block.id.find('_');

    if (underscore_index != string::npos) {
        id1 = block.id.substr(0, underscore_index + 1) + to_string(split_number);
        id2 = block.id.substr(0, underscore_index + 1) + to_string(split_number + 1);
    } else {
        id1 = block.id + '_' + to_string(split_number);
        id2 = block.id + '_' + to_string(split_number + 1);
    }

    // Creazione dei due nuovi blocchi
    Block part1 = {
        id1,
        block.label,
        vector<int>(sequence_ids.begin(), sequence_ids.begin() + split_point),
        block.begin_column,
        block.end_column
    };

    Block part2 = {
        id2,
        block.label,
        vector<int>(sequence_ids.begin() + split_point, sequence_ids.end()),
        block.begin_column,
        block.end_column
    };

    return make_tuple(part1, part2);
}


// Funzione per dividere un blocco per colonna
tuple<Block, Block> split_block_by_column(const Block& block, int split_number) {
    const string& label = block.label;
    int begin_column = block.begin_column;
    int end_column = block.end_column;

    if (end_column - begin_column <= 0) {
        throw invalid_argument("Block must span at least two columns to be split.");
    }

    // Generatore di numeri casuali per il punto di divisione
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> dist(begin_column + 1, end_column);
    int split_point = dist(gen);

    // Creazione dei nuovi ID dei blocchi
    string id1, id2;
    size_t underscore_index = block.id.find('_');

    if (underscore_index != string::npos) {
        id1 = block.id.substr(0, underscore_index + 1) + to_string(split_number);
        id2 = block.id.substr(0, underscore_index + 1) + to_string(split_number + 1);
    } else {
        id1 = block.id + '_' + to_string(split_number);
        id2 = block.id + '_' + to_string(split_number + 1);
    }

    // Creazione dei due nuovi blocchi
    Block part1 = {
        id1,
        label.substr(0, split_point - begin_column),
        block.sequence_ids,
        begin_column,
        split_point - 1
    };

    Block part2 = {
        id2,
        label.substr(split_point - begin_column),
        block.sequence_ids,
        split_point,
        end_column
    };

    return make_tuple(part1, part2);
}


unordered_map<string, Block> greedy_row_merge(unordered_map<string, Block>& block_dict, vector<vector<string>>& block_id_matrix) {
    size_t num_seq = block_id_matrix.size();
    size_t num_bases = block_id_matrix[0].size();

    for (size_t i = 0; i < num_seq; ++i) {
        for (size_t j = 0; j < num_bases; ++j) {
            for (size_t z = 0; z < num_seq; ++z) {
                string id_1 = block_id_matrix[i][j];
                string id_2 = block_id_matrix[z][j];

                if (id_1 != id_2 && block_dict[id_1].label == block_dict[id_2].label) {
                    Block& block_1 = block_dict[id_1];
                    Block& block_2 = block_dict[id_2];
                    
                    Block new_block;
                    if (z < i) {
                        new_block = merge_two_blocks(block_2, block_1, "row_union");
                        block_dict = update_block_dict_with_same_id(block_dict, new_block, block_2.id, block_1.id);
                    } else {
                        new_block = merge_two_blocks(block_1, block_2, "row_union");
                        block_dict = update_block_dict_with_same_id(block_dict, new_block, block_1.id, block_2.id);
                    }

                    for (int sequence : new_block.sequence_ids) {
                        for (int column = new_block.begin_column; column <= new_block.end_column; ++column) {
                            block_id_matrix[sequence][column] = new_block.id;
                        }
                    }
                }
            }
        }
    }

    return block_dict;
}


int main() {
    // Creazione di un dizionario di blocchi e una matrice di test
    unordered_map<string, Block> block_dict = {
        {"B0.0", {"B0.0", "A", {0}, 0, 0}},
        {"B1.0", {"B1.0", "A", {1}, 0, 0}},
        {"B2.0", {"B2.0", "A", {2}, 0, 0}},
        {"B0.1", {"B0.1", "T", {0}, 1, 1}},
        {"B1.1", {"B1.1", "T", {1}, 1, 1}},
        {"B2.1", {"B2.1", "T", {2}, 1, 1}}
    };

    // Matrice di ID blocco
    vector<vector<string>> block_id_matrix = {
        {"B0.0", "B0.1"},
        {"B1.0", "B1.1"},
        {"B2.0", "B2.1"}
    };

    // Test della funzione greedy_row_merge
    block_dict = greedy_row_merge(block_dict, block_id_matrix);

    // Stampa dei risultati
    cout << "Block ID Matrix after row merge:\n";
    for (const auto& row : block_id_matrix) {
        for (const auto& id : row) {
            cout << id << " ";
        }
        cout << endl;
    }

    // Stampa dei blocchi dopo il merge
    cout << "\nBlocks after row merge:\n";
    for (const auto& pair : block_dict) {
        print_block(pair.second);
    }


    return 0;
}