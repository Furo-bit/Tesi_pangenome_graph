#ifndef UTILS_H
#define UTILS_H

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
#include <graphviz/gvc.h>
#include <graphviz/cgraph.h> 
#include <set>
#include <map>

using namespace std;

struct Block {
    string id;
    string label;
    vector<int> sequence_ids;
    int begin_column;
    int end_column;
};
Block create_block(const string& id, const string& label, const vector<int>& sequence_ids, int begin_column, int end_column);
void print_block(const Block& block);
unordered_map<string, vector<char>> read_fasta(const string& filename);
vector<vector<char>> sequences_to_matrix(const unordered_map<string, vector<char>>& sequences);
pair<unordered_map<string, Block>, vector<vector<string>>> build_blocks_from_sequence_matrix(const vector<vector<char>>& sequence_matrix);
vector<Block> sort_blocks_by_sequence_id(const unordered_map<string, Block>& block_dict);
Block merge_two_blocks(const Block& block_1, const Block& block_2, const string& how);
unordered_map<string, Block> update_block_dict_with_same_id(unordered_map<string, Block>& block_dict, const Block& new_block, const string& id1, const string& id2);
void update_block_submatrix_with_same_id(vector<vector<string>>& block_id_matrix, const string& new_id, const vector<int>& sequences, int first_column, int last_column);
double acceptance_probability(double delta, double temperature, double beta);
vector<int> generate_random_numbers(int seed, int start, int end, int count);
int generate_random_number(int seed, int start, int end);
tuple<Block, Block> split_block_by_row(const Block& block, int split_number);
tuple<Block, Block> split_block_by_column(const Block& block, int split_number);
unordered_map<string, Block> greedy_row_merge(unordered_map<string, Block>& block_dict, vector<vector<string>>& block_id_matrix);
int of_min_label_length_threshold(int threshold, int penalization, int label_length);
int of_pangeblocks(int threshold, int penalization, int label_length);
bool check_number_in_string(int number, const string& str);
bool check_id_matrix_consistency(const vector<vector<string>>& block_id_matrix, const unordered_map<string, Block>& block_dict, const vector<vector<char>>& sequence_matrix);
void graph_to_gfa(Agraph_t* g, const string& filename);
void remove_indel_nodes(Agraph_t* g);
Agraph_t* build_graph(const unordered_map<string, Block>& block_dict, const vector<vector<string>>& block_id_matrix);
char* toCharPointer(const string& str);
void print_graph(Agraph_t* g);
void remove_nodes_with_dash_labels(Agraph_t *g);
unordered_map<string, string> read_config(const string& filename);
string remove_chars(const std::string& str, const std::string& chars_to_remove);

#endif