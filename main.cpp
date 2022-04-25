#include <iostream>
#include <vector>
#include <set>
#include <limits>
#include <cmath>
#include <fstream>
#include <string>
#include <chrono>
#include <numeric>
#include <algorithm>

using namespace std;

vector<string> split_file_line(const string &str) {
    vector<string> tokens;

    string::size_type start = 0;
    string::size_type end;

    while ((end = str.find('\t', start)) != string::npos) {
        tokens.push_back(str.substr(start, end - start));

        start = end + 1;
    }

    tokens.push_back(str.substr(start));

    return tokens;
}

// Base interface
class TSPInterface {
private:
    virtual void run_algorithm(unsigned int curr, unsigned int count, double cost) = 0;

protected:
    const unsigned int nodes;
    unsigned int tried_paths;
    unsigned int tried_edges;
    double min_cost{std::numeric_limits<float>::max()};
    vector<unsigned int> min_path;
    vector<bool> visited_nodes;
    vector<unsigned int> current_path;
    vector<vector<double>> adj_matrix;

    inline explicit TSPInterface(unsigned const int n) : nodes(n), adj_matrix(n, vector<double>(n, 0)) {
        visited_nodes = vector<bool>(n);
        current_path = vector<unsigned int>(n);
        tried_paths = 0;
        tried_edges = 0;
    }

public:

    virtual void update_edge_cost(unsigned int i, unsigned int j, double dist) = 0;

    virtual void start_algorithm() {
        visited_nodes[0] = true;
        run_algorithm(0, 1, 0);
    }

    inline void reset_state() {
        visited_nodes = vector<bool>(nodes);
        current_path = vector<unsigned int>(nodes);
        tried_paths = 0;
        tried_edges = 0;
        min_cost = std::numeric_limits<float>::max();
    }

    string get_min_path() {
        string str = to_string(min_cost) + ',';
        for (unsigned int i = 0; i < nodes; i++) {
            str += to_string(min_path[i]) + "->";
        }
        str += to_string(min_path[0]) + ',' + to_string(tried_paths) + ',' + to_string(tried_edges);
        return str;
    }
};

// Basic TSP, if partial path length + edge cost < current minimal complete path length then don't take it
class TSPBasicHeuristic : public TSPInterface {
private:
    void run_algorithm(const unsigned int curr, const unsigned int count, const double cost) override {
        if (count == nodes) {
            tried_paths += 1;
            const double candidate_cost = cost + adj_matrix[curr][0];
            if (candidate_cost < min_cost) {
                min_cost = candidate_cost;
                min_path = current_path;
            }
            return;
        }
        for (unsigned int i = 0; i < nodes; i++) {
            if (!visited_nodes[i]) {
                const double next_cost = cost + adj_matrix[curr][i];
                if (next_cost < min_cost) {
                    visited_nodes[i] = true;
                    current_path[count] = i;
                    tried_edges += 1;
                    run_algorithm(i, count + 1, next_cost);
                    visited_nodes[i] = false;
                }
            }
        }
    }

public:
    explicit TSPBasicHeuristic(unsigned const int n) : TSPInterface(n) {}

    inline void update_edge_cost(unsigned int i, unsigned int j, double dist) override {
        adj_matrix[i][j] = dist;
        adj_matrix[j][i] = dist;
    }
};

// Auxiliary struct for sorted TSPs
struct AdjacentNode {
public:
    const unsigned int node;
    const double cost;

    inline AdjacentNode(const double cost, const unsigned int node) : node(node), cost(cost) {}

    inline bool operator<(AdjacentNode other) const {
        return cost < other.cost;
    }
};

// Basic TSP but adjacent nodes are checked ordered by edge calculate_distance
class TSPBasicHeuristicSorted : public TSPInterface {
private:
    vector<set<AdjacentNode>> nodes_sorted_by_distance;

    void run_algorithm(const unsigned int curr, const unsigned int count, const double cost) override {
        if (count == nodes) {
            tried_paths += 1;
            const double candidate_cost = cost + adj_matrix[curr][0];
            if (candidate_cost < min_cost) {
                min_cost = candidate_cost;
                min_path = current_path;
            }
            return;
        }
        for (auto adj_node: nodes_sorted_by_distance[curr]) {
            const auto i = adj_node.node;
            if (!visited_nodes[i]) {
                if (cost + adj_node.cost < min_cost) {
                    visited_nodes[i] = true;
                    current_path[count] = i;
                    tried_edges += 1;
                    run_algorithm(i, count + 1, cost + adj_matrix[curr][i]);
                    visited_nodes[i] = false;
                }
            }
        }
    }

public:

    inline explicit TSPBasicHeuristicSorted(unsigned const int n) : TSPInterface(n), nodes_sorted_by_distance(n) {
    }

    inline void update_edge_cost(unsigned int i, unsigned int j, double dist) override {
        adj_matrix[i][j] = dist;
        adj_matrix[j][i] = dist;
        nodes_sorted_by_distance[i].emplace(dist , j);
        nodes_sorted_by_distance[j].emplace(dist , i);
    }
};

// TSP with better heuristic, if partial path length + edge cost + return-to-node-0 cost < current minimal complete path length then don't take it
class TSPReturnHeuristic : public virtual TSPInterface {
private:
    void run_algorithm(const unsigned int curr, const unsigned int count, const double cost) override {
        if (count == nodes) {
            tried_paths += 1;
            const double candidate_cost = cost + adj_matrix[curr][0];
            if (candidate_cost < min_cost) {
                min_cost = candidate_cost;
                min_path = current_path;
            }
            return;
        }
        for (unsigned int i = 0; i < nodes; i++) {
            if (!visited_nodes[i]) {
                const double next_cost = cost + adj_matrix[curr][i];
                if (next_cost + adj_matrix[i][0] < min_cost) {
                    visited_nodes[i] = true;
                    current_path[count] = i;
                    tried_edges += 1;
                    run_algorithm(i, count + 1, next_cost);
                    visited_nodes[i] = false;
                }
            }
        }
    }

public:
    inline explicit TSPReturnHeuristic(unsigned const int n) : TSPInterface(n) {}

    inline void update_edge_cost(unsigned int i, unsigned int j, double dist) override {
        adj_matrix[i][j] = dist;
        adj_matrix[j][i] = dist;
    }

};

// TSP with return-heuristic but adjacent nodes are ordered by return-heuristic. Can stop checking earlier.
class TSPReturnHeuristicSorted : public TSPInterface {
private:
    vector<set<AdjacentNode>> nodes_sorted_by_heuristic;

    void run_algorithm(const unsigned int curr, const unsigned int count, const double cost) override {
        if (count == nodes) {
            tried_paths += 1;
            const double candidate_cost = cost + adj_matrix[curr][0];
            if (candidate_cost < min_cost) {
                min_cost = candidate_cost;
                min_path = current_path;
            }
            return;
        }
        for (auto adj_node: nodes_sorted_by_heuristic[curr]) {
            const auto i = adj_node.node;
            if (!visited_nodes[i]) {
                if (cost + adj_node.cost < min_cost) {
                    visited_nodes[i] = true;
                    current_path[count] = i;
                    tried_edges += 1;
                    run_algorithm(i, count + 1, cost + adj_matrix[curr][i]);
                    visited_nodes[i] = false;
                } else {
                    break;
                }
            }
        }
    }

public:


    inline explicit TSPReturnHeuristicSorted(unsigned const int n) : TSPInterface(n), nodes_sorted_by_heuristic(n) {
    }

    inline void update_edge_cost(unsigned int i, unsigned int j, double dist) override {
        adj_matrix[i][j] = dist;
        adj_matrix[j][i] = dist;
        nodes_sorted_by_heuristic[i].emplace(dist + adj_matrix[0][j], j);
        nodes_sorted_by_heuristic[j].emplace(dist + adj_matrix[0][i], i);
    }
};

// Greedy TSP that applies 2-Opt afterwards
class TSP2Opt : public TSPInterface {
private:
    vector<set<AdjacentNode>> nodes_sorted_by_distance;
    bool path_changed{true};

    void calculate_min_cost(){
        min_cost = adj_matrix[min_path[nodes-1]][min_path[0]];
        for(unsigned int i = 0; i < nodes-1; i++){
            min_cost += adj_matrix[min_path[i]][min_path[i+1]];
        }
    }


    void run_algorithm(const unsigned int curr, const unsigned int count, const double cost) override {
        if (count == nodes) {
            tried_paths += 1;
            return;
        }
        for (auto adj_node: nodes_sorted_by_distance[curr]) {
            const auto i = adj_node.node;
            if (!visited_nodes[i]) {
                visited_nodes[i] = true;
                current_path[count] = i;
                tried_edges += 1;
                run_algorithm(i, count + 1, 0);
                break;
            }
        }
    }

public:
    void start_algorithm() override {
        visited_nodes[0] = true;
        run_algorithm(0, 1, 0);
        while(path_changed){
            path_changed = false;
            for(auto i = 1; i < nodes-1; i++){
                for(auto j = i+1; j < nodes-1; j++){
                    if(adj_matrix[current_path[i-1]][current_path[i]] + adj_matrix[current_path[j]][current_path[j + 1]] > adj_matrix[current_path[i-1]][current_path[j]] + adj_matrix[current_path[i]][current_path[j+1]]){
                        std::reverse(current_path.begin()+i, current_path.begin() + j);
                        path_changed=true;
                    }
                }
            }
        }
        min_path = current_path;
        calculate_min_cost();
    }

    inline explicit TSP2Opt(unsigned const int n) : TSPInterface(n), nodes_sorted_by_distance(n) {
    }

    inline void update_edge_cost(unsigned int i, unsigned int j, double dist) override {
        adj_matrix[i][j] = dist;
        adj_matrix[j][i] = dist;
        nodes_sorted_by_distance[i].emplace(dist , j);
        nodes_sorted_by_distance[j].emplace(dist , i);
    }
};

// Auxiliary struct to store points from input
struct Point {
    const double x;
    const double y;
public:
    inline Point(const double x, const double y) : x(x), y(y) {}

    [[nodiscard]] inline double calculate_distance(Point other) const {
        auto dx = x - other.x;
        auto dy = y - other.y;
        return (dx * dx) + (dy * dy);
    }
};

// Runs function
int main(int argc, char **argv) {
    string usage = "    Usage: ./TSP_Euclid {file} {algorithm} *{nodes} *{times}\n"
                   "         file: contains a first line with nodes, number of cities, and then nodes lines with points \"x y\"\n"
                   "         algorithm: 1 for Simple, 2 for Heuristic, 3 for Heuristic + Sorted\n"
                   "         (optional) nodes: number of cities, overrides the nodes of the file if it is less\n"
                   "         (optional) times: number of executions to time, default 60\n";


    if (argc < 4 || argc > 5) {
        cerr << "Not enough or too many arguments.\n" << usage;
        return -1;
    }

    fstream input;
    input.open(argv[1], ios::in); //open a file to perform read operation using file object
    if (!input) {
        cerr << "Can't open file.\n" << usage;
        return -1;
    }

    if (input.peek() == EOF) {
        cerr << "File empty.\n" << usage;
        return -1;
    }

    string temp_line;
    getline(input, temp_line);
    vector<string> parsed_first_line = split_file_line(temp_line);

    int n = stoi(parsed_first_line[0]);
    if (argc >= 4) {
        n = atoi(argv[3]);
    }

    unsigned int times = 60;
    if (argc == 5) {
        times = atoi(argv[4]);
    }


    // file to vector
    vector<Point> temp_points;
    for (unsigned int i = 0; i < n; i++) {
        getline(input, temp_line);
        vector<string> parsed_line = split_file_line(temp_line);
        auto x = stod(parsed_line[0]);
        auto y = stod(parsed_line[1]);
        temp_points.emplace_back(x, y);
    }
    ifstream output_test("output.csv");
    bool flag = false;
    if (!output_test) {
        flag = true;
    }
    output_test.close();


    ofstream output("output.csv", std::ios_base::app);
    if (flag) {
        output << "algorithm,nodes,init_time,run_time,total_time,min_cost,min_path,complete_paths,tried_edges\n";
    }
    unsigned int algorithm = atoi(argv[2]);
    string output_str;
    if (algorithm == 1) {
        auto tsp = TSPBasicHeuristic(n);
        vector<double> init_times;
        vector<double> run_times;
        for (unsigned int t = 0; t < times; t++) {
            tsp.reset_state();
            auto t1 = chrono::steady_clock::now();
            for (unsigned int i = 0; i < n; i++) {
                for (unsigned int j = 0; j < i; j++) {
                    auto dist = sqrt(temp_points[i].calculate_distance(temp_points[j]));
                    tsp.update_edge_cost(i, j, dist);
                }
            }
            auto t2 = chrono::steady_clock::now();
            tsp.start_algorithm();
            auto t3 = chrono::steady_clock::now();
            init_times.emplace_back(std::chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count());
            run_times.emplace_back(std::chrono::duration_cast<chrono::nanoseconds>(t3 - t2).count());
        }
        double sum1 = accumulate(begin(init_times), end(init_times), 0.0);
        double m1 = sum1 / times;
        double sum2 = accumulate(begin(run_times), end(run_times), 0.0);
        double m2 = sum2 / times;
        output << "TSPBasicHeuristic," + to_string(n) + ',' + to_string(m1) + ',' + to_string(m2)
                  + ',' + to_string(m1 + m2) + ',' + tsp.get_min_path() << '\n';
    } else if (algorithm == 2) {
        auto tsp = TSPBasicHeuristicSorted(n);
        vector<double> init_times;
        vector<double> run_times;
        for (unsigned int t = 0; t < times; t++) {
            tsp.reset_state();
            auto t1 = chrono::steady_clock::now();
            for (unsigned int i = 0; i < n; i++) {
                for (unsigned int j = 0; j < i; j++) {
                    auto dist = sqrt(temp_points[i].calculate_distance(temp_points[j]));
                    tsp.update_edge_cost(i, j, dist);
                }
            }
            auto t2 = chrono::steady_clock::now();
            tsp.start_algorithm();
            auto t3 = chrono::steady_clock::now();
            init_times.emplace_back(std::chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count());
            run_times.emplace_back(std::chrono::duration_cast<chrono::nanoseconds>(t3 - t2).count());
        }
        double sum1 = accumulate(begin(init_times), end(init_times), 0.0);
        double m1 = sum1 / times;
        double sum2 = accumulate(begin(run_times), end(run_times), 0.0);
        double m2 = sum2 / times;
        output << "TSPBasicHeuristicSorted," + to_string(n) + ',' + to_string(m1) + ',' + to_string(m2)
                  + ',' + to_string(m1 + m2) + ',' + tsp.get_min_path() << '\n';
    } else if (algorithm == 3) {
        auto tsp = TSPReturnHeuristic(n);
        vector<double> init_times;
        vector<double> run_times;
        for (unsigned int t = 0; t < times; t++) {
            tsp.reset_state();
            auto t1 = chrono::steady_clock::now();
            for (unsigned int i = 0; i < n; i++) {
                for (unsigned int j = 0; j < i; j++) {
                    auto dist = sqrt(temp_points[i].calculate_distance(temp_points[j]));
                    tsp.update_edge_cost(i, j, dist);
                }
            }
            auto t2 = chrono::steady_clock::now();
            tsp.start_algorithm();
            auto t3 = chrono::steady_clock::now();
            init_times.emplace_back(std::chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count());
            run_times.emplace_back(std::chrono::duration_cast<chrono::nanoseconds>(t3 - t2).count());
        }
        double sum1 = accumulate(begin(init_times), end(init_times), 0.0);
        double m1 = sum1 / times;
        double sum2 = accumulate(begin(run_times), end(run_times), 0.0);
        double m2 = sum2 / times;


        output << "TSPReturnHeuristic," + to_string(n) + ',' + to_string(m1) + ',' + to_string(m2)
                  + ',' + to_string(m1 + m2) + ',' + tsp.get_min_path() << '\n';
    }else if (algorithm == 4) {
        auto tsp = TSPReturnHeuristicSorted(n);
        vector<double> init_times;
        vector<double> run_times;
        for (unsigned int t = 0; t < times; t++) {
            tsp.reset_state();
            auto t1 = chrono::steady_clock::now();
            for (unsigned int i = 0; i < n; i++) {
                for (unsigned int j = 0; j < i; j++) {
                    auto dist = sqrt(temp_points[i].calculate_distance(temp_points[j]));
                    tsp.update_edge_cost(i, j, dist);
                }
            }
            auto t2 = chrono::steady_clock::now();
            tsp.start_algorithm();
            auto t3 = chrono::steady_clock::now();
            init_times.emplace_back(std::chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count());
            run_times.emplace_back(std::chrono::duration_cast<chrono::nanoseconds>(t3 - t2).count());
        }
        double sum1 = accumulate(begin(init_times), end(init_times), 0.0);
        double m1 = sum1 / times;
        double sum2 = accumulate(begin(run_times), end(run_times), 0.0);
        double m2 = sum2 / times;


        output << "TSPReturnHeuristicSorted," + to_string(n) + ',' + to_string(m1) + ',' + to_string(m2)
                  + ',' + to_string(m1 + m2) + ',' + tsp.get_min_path() << '\n';
    }else if (algorithm == 5) {
        auto tsp = TSP2Opt(n);
        vector<double> init_times;
        vector<double> run_times;
        for (unsigned int t = 0; t < times; t++) {
            tsp.reset_state();
            auto t1 = chrono::steady_clock::now();
            for (unsigned int i = 0; i < n; i++) {
                for (unsigned int j = 0; j < i; j++) {
                    auto dist = sqrt(temp_points[i].calculate_distance(temp_points[j]));
                    tsp.update_edge_cost(i, j, dist);
                }
            }
            auto t2 = chrono::steady_clock::now();
            tsp.start_algorithm();
            auto t3 = chrono::steady_clock::now();
            init_times.emplace_back(std::chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count());
            run_times.emplace_back(std::chrono::duration_cast<chrono::nanoseconds>(t3 - t2).count());
        }
        double sum1 = accumulate(begin(init_times), end(init_times), 0.0);
        double m1 = sum1 / times;
        double sum2 = accumulate(begin(run_times), end(run_times), 0.0);
        double m2 = sum2 / times;


        output << "TSP2Opt," + to_string(n) + ',' + to_string(m1) + ',' + to_string(m2)
                  + ',' + to_string(m1 + m2) + ',' + tsp.get_min_path() << '\n';
    } else {
        cerr << "Wrong algorithm.\n" << usage;
        return -1;
    }
    output.close();
    return 0;
}