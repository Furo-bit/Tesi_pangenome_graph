#include "utils.h"
#include <chrono>
#include <sys/resource.h> 

using namespace std::chrono;

// Funzione per convertire stringa in un tipo numerico
template<typename T>
T convert_to(const string& value) {
    T result;
    stringstream(value) >> result; // Converte la stringa in un valore del tipo T
    return result;
}



void simulated_annealing(unordered_map<string, Block>& block_dict, vector<vector<string>>& block_id_matrix){

    string params = "params_simulated_annealing_c.txt";
    unordered_map<string, string> config = read_config(params);
    int split_number = 0;
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
        //print_block(block_1);

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
        
        vector<vector<string>> list_mergeable_vertical_greedy;
        vector<string> list_used_ids;
        for (int i = limit_vert_top; i < limit_vert_bot; ++i){
            vector<string> list_i_mergeable;
            Block block_a = block_dict.at(block_id_matrix[i][col_index]);

            // Verifica se block_a["id"] non è in list_used_ids
            if (find(list_used_ids.begin(), list_used_ids.end(), block_a.id) == list_used_ids.end()) {

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

                                if (find(list_used_ids.begin(), list_used_ids.end(), block_a.id) == list_used_ids.end()) {
                                    list_used_ids.push_back(block_a.id);
                                }
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
            for (vector<string> comp_list : list_mergeable_vertical_greedy) {
                string& block_id = comp_list[0];
                Block& block_greedy = block_dict.at(block_id);
                string block_greedy_label = block_greedy.label;
                
                // Rimuovi i caratteri '-' dalla label
                block_greedy_label = remove_chars(block_greedy_label, "-");

                // Calcola il guadagno usando la funzione obiettivo di pangeblocks
                int block_greedy_label_value = of_pangeblocks(threshold, penalization, block_greedy_label.length());
                vertical_greedy_merge_gain += - (block_greedy_label_value * (comp_list.size()- 1));
            }
            operation_gains["row_merge_greedy"] = vertical_greedy_merge_gain;
        }

        // Column merge
        // Check left
        int limit_left;
        if (block_1.begin_column <= limit_horizontal) {
            limit_left = 0;
        }
        else {
            limit_left = block_1.begin_column - limit_horizontal;
        }

        for (int column = block_1.begin_column - 1; column >= limit_left; --column){
            Block block_2 = block_dict.at(block_id_matrix[row_index][column]);
            if (block_1.id != block_2.id && 
                find(list_mergeable_horizontal_left.begin(), list_mergeable_horizontal_left.end(), block_2.id) == list_mergeable_horizontal_left.end()){
                    if (block_1.sequence_ids == block_2.sequence_ids){
                        list_mergeable_horizontal_left.insert(list_mergeable_horizontal_left.begin(), block_2.id);
                    }
                    else{
                        break;
                    }
                }
        }

        // Check right
        int limit_right;
        if (block_1.end_column + limit_horizontal >= cols){
            limit_right = cols;
        }
        else{
            limit_right = block_1.end_column + limit_horizontal;
        }
        for (int column = block_1.end_column; column < limit_right; ++column){ 
            Block block_2 = block_dict.at(block_id_matrix[row_index][column]);
            if (block_1.id != block_2.id && 
                find(list_mergeable_horizontal_right.begin(), list_mergeable_horizontal_right.end(), block_2.id) == list_mergeable_horizontal_right.end()){
                    if (block_1.sequence_ids == block_2.sequence_ids){
                        list_mergeable_horizontal_right.push_back(block_2.id);
                    }
                    else{
                        break;
                    }
                }
        }

        // Calculate gain
        // Unisce il vettore sinistro
        list_mergeable_horizontal.insert(list_mergeable_horizontal.end(),
                                        list_mergeable_horizontal_left.begin(),
                                        list_mergeable_horizontal_left.end());

        // Inserisce block_1["id"]
        list_mergeable_horizontal.push_back(block_1.id);

        // Unisce il vettore destro
        list_mergeable_horizontal.insert(list_mergeable_horizontal.end(),
                                        list_mergeable_horizontal_right.begin(),
                                        list_mergeable_horizontal_right.end());


        int horizontal_merge_gain;
        Block new_block_horizontal;
        if (!list_mergeable_horizontal_left.empty() || !list_mergeable_horizontal_right.empty()) {
            new_block_horizontal = block_dict.at(list_mergeable_horizontal[0]); 
            string new_block_horizontal_label = new_block_horizontal.label;
            new_block_horizontal_label  = remove_chars(new_block_horizontal_label, "-");
            horizontal_merge_gain = of_pangeblocks(threshold, penalization, new_block_horizontal_label.size());
            for (string& id : list_mergeable_horizontal) {
                if (new_block_horizontal.id != id) {
                    Block block_2 = block_dict.at(id); 
                    string block_2_label = block_2.label;
                    block_2_label  = remove_chars(block_2_label, "-");
                    int block_2_label_value = of_pangeblocks(threshold, penalization, block_2_label.size());    
                    horizontal_merge_gain += block_2_label_value;
                    new_block_horizontal = merge_two_blocks(new_block_horizontal, block_2, "column_union");
                }
            }
        
            new_block_horizontal_label = new_block_horizontal.label;
            new_block_horizontal_label  = remove_chars(new_block_horizontal_label, "-");
            int new_block_horizontal_label_value = of_pangeblocks(threshold, penalization, new_block_horizontal_label.size());
            horizontal_merge_gain = new_block_horizontal_label_value - horizontal_merge_gain;
            operation_gains["column_merge"] = horizontal_merge_gain;
        }



        Block new_block_1r;
        Block new_block_2r;


        // Check row split
        if (block_1.sequence_ids.size() > 1) {
            
            uniform_int_distribution<int> distribution_row_split(1, block_1.sequence_ids.size() - 1);
            int split_point_row = distribution_row_split(generator);

            tie(new_block_1r, new_block_2r) = split_block_by_row(block_1, split_number, split_point_row);
            split_number += 2;

            string label_nb1r = new_block_1r.label;
            string label_nb2r = new_block_2r.label;

            label_nb1r = (label_nb1r, "-");
            label_nb2r = (label_nb2r, "-");

            int delta = (of_pangeblocks(threshold, penalization, label_nb1r.size()) +
                        of_pangeblocks(threshold, penalization, label_nb2r.size()) -
                        block_1_label_value);
            
            operation_gains["row_split"] = delta;
        }
        

        Block new_block_1c;
        Block new_block_2c;

        // Check column split
        if (block_1.end_column - block_1.begin_column > 0){
            
            uniform_int_distribution<int> distribution_column_split(block_1.begin_column + 1, block_1.end_column);
            int split_point_column = distribution_column_split(generator);

            tie(new_block_1c, new_block_2c) = split_block_by_column(block_1, split_number, split_point_column);
            split_number += 2;

            string label_nb1c = new_block_1c.label;
            string label_nb2c = new_block_2c.label;

            label_nb1c = (label_nb1c, "-");
            label_nb2c = (label_nb2c, "-");


            int delta = (of_pangeblocks(threshold, penalization, label_nb1c.size()) +
                        of_pangeblocks(threshold, penalization, label_nb2c.size()) -
                        block_1_label_value);
            
            operation_gains["column_split"] = delta;

        }

        // Select operation
        if (!operation_gains.empty()) {
            // Converti unordered_map in un vettore di coppie per poter ordinare
            vector<pair<string, int>> sorted_gains(operation_gains.begin(), operation_gains.end());

            // Ordina il vettore per il secondo elemento (gain)
            sort(sorted_gains.begin(), sorted_gains.end(),
                    [](pair<string, int>& a, pair<string, int>& b) {
                        return a.second < b.second;
                    });

            // Prendi la prima operazione e il suo guadagno dopo aver ordinato
            string operation = sorted_gains.front().first;
            float gain = sorted_gains.front().second;
 

            // Aggiorna il vettore operations_for_cooling
            operations_for_cooling.insert(operations_for_cooling.begin(), operation);
            operations_for_cooling.pop_back(); // Rimuove l'ultimo elemento

            float probability = 1.0f;
            
            if (operation == "column_split" || operation == "row_split"){
                probability = acceptance_probability(gain, temperature, beta);
           
            }

            if (operation == "column_merge")
            {

                // Column merge
                for (auto& id : list_mergeable_horizontal) {
                    if (id != new_block_horizontal.id)
                    {
                        block_dict.erase(id);  // Rimuove l'elemento dal dizionario
                    }                 
                }

                // Assegna il nuovo blocco orizzontale al dizionario
                block_dict[new_block_horizontal.id] = new_block_horizontal; 

                // Aggiorna il block_id_matrix per ciascuna sequenza e colonna
                for (auto& sequence : new_block_horizontal.sequence_ids) {
                    for (int column = new_block_horizontal.begin_column; column <= new_block_horizontal.end_column; ++column) {
                        block_id_matrix[sequence][column] = new_block_horizontal.id;
                    }
                }
            }
            else if (operation == "row_split"){
                float x = distribution(generator) / static_cast<float>(cell_total);
                if (x < probability) {
                    
                    // Row split - Aggiorna new_block_1r
                    for (auto& sequence : new_block_1r.sequence_ids) {
                        for (int column = new_block_1r.begin_column; column <= new_block_1r.end_column; ++column) {
                            block_id_matrix[sequence][column] = new_block_1r.id;
                        }
                    }

                    // Row split - Aggiorna new_block_2r
                    for (auto& sequence : new_block_2r.sequence_ids) {
                        for (int column = new_block_2r.begin_column; column <= new_block_2r.end_column; ++column) {
                            block_id_matrix[sequence][column] = new_block_2r.id;
                        }
                    }

                    // Aggiorna il dizionario
                    block_dict.erase(block_1.id);  // Elimina il vecchio blocco
                    block_dict[new_block_1r.id] = new_block_1r;  // Aggiungi new_block_1r
                    block_dict[new_block_2r.id] = new_block_2r;  // Aggiungi new_block_2r
                }

            }
            else if (operation == "column_split"){
                float x = distribution(generator) / static_cast<float>(cell_total);
                if (x / static_cast<float>(cell_total) < probability) {
                    
                    // Column split - Aggiorna new_block_1c
                    for (auto& sequence : new_block_1c.sequence_ids) {
                        for (int column = new_block_1c.begin_column; column <= new_block_1c.end_column; ++column) {
                            block_id_matrix[sequence][column] = new_block_1c.id;
                        }
                    }

                    // Column split - Aggiorna new_block_2c
                    for (auto& sequence : new_block_2c.sequence_ids) {
                        for (int column = new_block_2c.begin_column; column <= new_block_2c.end_column; ++column) {
                            block_id_matrix[sequence][column] = new_block_2c.id;
                        }
                    }

                    // Aggiorna il dizionario
                    block_dict.erase(block_1.id);  // Elimina il vecchio blocco
                    block_dict[new_block_1c.id] = new_block_1c;  // Aggiungi new_block_1c
                    block_dict[new_block_2c.id] = new_block_2c;  // Aggiungi new_block_2c
                }
            }
            else if (operation == "row_merge_greedy"){
                for (auto& comp_list : list_mergeable_vertical_greedy) {
                    string block_1_id = comp_list[0];
                    Block& block_1 = block_dict[block_1_id];

                    for (string& block_2_id : comp_list) {
                        if (block_1_id != block_2_id) {
                            Block& block_2 = block_dict[block_2_id];
                            block_1 = merge_two_blocks(block_1, block_2, "row_union");
                            block_dict.erase(block_2_id);
                        }
                    }

                    block_dict[block_1_id] = block_1;

                    for (auto& sequence : block_1.sequence_ids) {
                        for (int column = block_1.begin_column; column <= block_1.end_column; ++column) {
                            block_id_matrix[sequence][column] = block_1.id;
                        }
                    }
                }
            }



        }


        /*
        for (auto& entry : operation_gains) {
            cout << entry.first << ": " << entry.second << endl;
        }
        */
        // Visualizza la matrice di ID
        //print_block_id_matrix(block_id_matrix);
        
        // Visualizza il dizionario di blocchi
        //print_block_dict(block_dict);

        //cout << "\n====================================\n\n";

        // Cooling condition
        // Inizializza i contatori per gli split di riga e colonna
        int column_split_number = 0;
        int row_split_number = 0;

        // Conta il numero di "column_split" e "row_split" nel vettore operations_for_cooling
        column_split_number = count(operations_for_cooling.begin(), operations_for_cooling.end(), "column_split");
        row_split_number = count(operations_for_cooling.begin(), operations_for_cooling.end(), "row_split");


        // Verifica se la condizione di cooling è soddisfatta
        if (column_split_number + row_split_number >= static_cast<int>((operations_for_cooling.size() / 100.0) * splits_percent_for_cooling)) {
            temperature *= cooling_factor;

            // Reinizializza il vettore operations_for_cooling
            operations_for_cooling = vector<string>(rows * 15, "new_start");
        }

        // Controlla se la temperatura ha raggiunto il livello di terminazione
        if (temperature <= termination_temperature) {
            break;
        }
    }
    
}

// Funzione per ottenere l'uso della memoria su Linux/macOS
long getMemoryUsage() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss;  // Restituisce l'uso massimo della memoria in kilobyte
}



int main() {

    auto start = high_resolution_clock::now();
    long initial_memory = getMemoryUsage();


    // Leggi le sequenze dal file FASTA
    const string alignment = "MSAs/30-sars-cov-2-ena.fa"; //10-sars-cov-2-ena
    auto sequences = read_fasta(alignment);

    // Converti le sequenze in una matrice
    auto sequence_matrix = sequences_to_matrix(sequences);
    
    // Costruisci i blocchi e la matrice degli ID dei blocchi
    auto [block_dict, block_id_matrix] = build_blocks_from_sequence_matrix(sequence_matrix);

    // Inserire qua euristica
    simulated_annealing(block_dict, block_id_matrix);

    //check_id_matrix_consistency(block_id_matrix, block_dict, sequence_matrix);
    // Costruisce il grafo e lo salva 
    Agraph_t* g = build_graph(block_dict, block_id_matrix);

    long final_memory = getMemoryUsage();
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    int total_seconds = duration.count();
    int minutes = total_seconds / 60;
    int seconds = total_seconds % 60;
    double memory_usage_mb = static_cast<double>(final_memory - initial_memory) / 1024.0;

    // Stampa dei risultati
    cout << "Tempo di esecuzione: " << minutes << " minuti e " << seconds << " secondi" << endl;
    cout << "Memoria utilizzata: " << memory_usage_mb << " MB" << endl;

    return 0;
}