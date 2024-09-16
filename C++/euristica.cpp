#include "utils.h"

// Funzione per convertire stringa in un tipo numerico
template<typename T>
T convert_to(const string& value) {
    T result;
    stringstream(value) >> result; // Converte la stringa in un valore del tipo T
    return result;
}

void simulated_annealing(const unordered_map<string, Block>& block_dict, const vector<vector<string>>& block_id_matrix){

    string params = "params_simulated_annealing_c.txt";
    unordered_map<string, string> config = read_config(params);

 // Variabili da configurare
    int seed = 0;
    float temperature = 0.0f;
    float termination_temperature = 0.0f;
    int limit_horizontal = 0;
    int limit_vertical = 0;
    float beta = 0.0f;
    float cooling_factor = 0.0f;
    int threshold = 0;
    int penalization = 0;
    int splits_percent_for_cooling = 1000;

    // Assegna i valori letti dalla mappa alle variabili
    if (config.find("seed") != config.end()) {
        seed = convert_to<int>(config["seed"]);
    }
    if (config.find("temperature") != config.end()) {
        temperature = convert_to<float>(config["temperature"]);
    }
    if (config.find("termination_temperature") != config.end()) {
        termination_temperature = convert_to<float>(config["termination_temperature"]);
    }
    if (config.find("limit_horizontal") != config.end()) {
        limit_horizontal = convert_to<int>(config["limit_horizontal"]);
    }
    if (config.find("limit_vertical") != config.end()) {
        limit_vertical = convert_to<int>(config["limit_vertical"]);
    }
    if (config.find("beta") != config.end()) {
        beta = convert_to<float>(config["beta"]);
    }
    if (config.find("cooling_factor") != config.end()) {
        cooling_factor = convert_to<float>(config["cooling_factor"]);
    }
    if (config.find("threshold") != config.end()) {
        threshold = convert_to<int>(config["threshold"]);
    }
    if (config.find("penalization") != config.end()) {
        penalization = convert_to<int>(config["penalization"]);
    }
    if (config.find("splits_percent") != config.end()) {
        splits_percent_for_cooling = convert_to<int>(config["splits_percent"]);
    }

    int rows = block_id_matrix.size();
    int cols = block_id_matrix[0].size();
    int cell_total = rows * cols;

    // Creazione di un vettore di dimensione rows * 15, con valori di default null (o nullptr per puntatori)
    vector<string> operations_for_cooling(rows * 15, "new_start");

    if (limit_horizontal == -1) {
        limit_horizontal = cols;
    } else if (limit_horizontal == -2) {
        limit_horizontal = rows / 2;
    }

    if (limit_vertical == -1) {
        limit_vertical = cols;
    } else if (limit_vertical == -2) {
        limit_vertical = rows / 2;
    }

    // Inizializza il generatore di numeri casuali con un seed
    mt19937 generator(seed);  // Usa il seed passato come parametro
    uniform_int_distribution<int> distribution(1, cell_total);


    while (true)
    { 
        // Genera un numero casuale tra 1 e cell_total
        int random_number = distribution(generator);

        // Trasforma il numero casuale in indici della matrice
        int row_index = (random_number - 1) / cols;
        int col_index = (random_number - 1) % cols;

        // Accedi al blocco corrispondente usando gli indici
        Block block_1 = block_dict.at(block_id_matrix[row_index][col_index]);

        vector<string> list_mergeable_horizontal;
        vector<string> list_mergeable_horizontal_left;
        vector<string> list_mergeable_horizontal_right;
        vector<string> list_mergeable_vertical;
        map<string, int> operation_gains;
       
        string block_1_label = block_1.label;
        block_1_label = remove_chars(block_1_label, "-");
        int block_1_label_value = of_pangeblocks(threshold, penalization, block_1_label.size());
        
        // greedy vertical merge
        int limit_vert_top = limit_vertical;
        int limit_vert_bot = limit_vertical;

        if (limit_vertical > row_index) {
            limit_vert_top = 0;
        } else {
            limit_vert_top = row_index - limit_vertical;
        }

        if (limit_vertical + row_index > rows) {
            limit_vert_bot = rows;
        } else {
            limit_vert_bot = row_index + limit_vertical;
        }
        
        vector<string> list_mergeable_vertical_greedy;
        vector<string> list_used_ids;
        for (int i = limit_vert_top; i < limit_vert_bot; ++i){
            vector<string> list_i_mergeable;
            Block block_a = block_dict.at(block_id_matrix[i][col_index]);

            // Verifica se block_a["id"] non è in list_used_ids
            if (find(list_used_ids.begin(), list_used_ids.end(), block_a_id) == list_used_ids.end()) {

                // Codice da eseguire se block_a["id"] non è trovato
                for (int z = i; z < limit_vert_bot; ++z){

                    Block block_b = block_dict.at(block_id_matrix[z][col_index]);
                    if (block_a.id != block_b.id && 
                        find(list_used_ids.begin(), list_used_ids.end(), block_b.id) == list_used_ids.end()) {

                        if (block_a.begin_column == block_b.begin_column &&
                            block_a.end_column == block_b.end_column &&
                            block_a.label == block_b.label){
                                
                                list_i_mergeable.push_back(block_b.id);
                                list_used_ids.push_back(block_b.id);

                                if (find(list_used_ids.begin(), list_used_ids.end(), block_a_id) == list_used_ids.end()) {
                                    list_used_ids.push_back(block_a.id);
                            }
                        }
                }
            }

            // Verifica se list_i_mergeable non è vuoto
            if (!list_i_mergeable.empty()) {
                // Aggiungi l'id di block_a all'inizio di list_i_mergeable
                list_i_mergeable.insert(list_i_mergeable.begin(), block_a.id);

                // Aggiungi list_i_mergeable a list_mergeable_vertical_greedy
                list_mergeable_vertical_greedy.push_back(list_i_mergeable);
            }

        }
        // Calculate greedy vertical gain
        if (!list_mergeable_vertical_greedy.empty()) {
            int vertical_greedy_merge_gain = 0;
            // continua qua
        }

        break;
    }
    

}



int main() {

    // Leggi le sequenze dal file FASTA
    const string alignment = "test.fa";
    auto sequences = read_fasta(alignment);

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