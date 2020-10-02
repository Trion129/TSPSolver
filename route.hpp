#include "datatypes.hpp"
#include "tour-tree.hpp"

struct Route {
    vector<Point*> cityData;
    long double length;
    SplayTree* tree;
    Graph* graph;

    explicit Route(const vector<Point*>& cities, vector<int>& city_sequence, Graph* graph_given) : length(0), graph(graph_given) {
        cityData = cities;
        Point* previousCity = cityData[city_sequence[city_sequence.size()-1]];
        for (auto city : city_sequence) {
            Point* point = cityData[city];
            length += euclidean(previousCity, point);
            previousCity = point;
        }
        tree = new SplayTree(city_sequence);
    }
    explicit Route(Route* route) : cityData(vector<Point*>(route->cityData.size())), length(0) {
        tree = new SplayTree(route->tree);
        graph = route->graph;
        cityData = route->cityData;
        auto city_sequence = tree->get_sequence();
        Point* previousCity = cityData[(*city_sequence)[city_sequence->size()-1]];
        for (auto city : *city_sequence) {
            Point* point = cityData[city];
            length += euclidean(previousCity, point);
            previousCity = point;
        }
    }

    ~Route() {
        delete tree;
    }

    Point* get_point(int city) {
        return cityData[city];
    }

    void swap_edges(int city1, int city2) const {
        tree->flip(city1, city2);
    }

    void improve_route(int t1, int current_city, int depth, long double gain, long double& max_gain_ptr,
                       vector<pair<int, int>>& swaps_done, vector<pair<int, int>>& best_swaps_done,
                       unordered_set<int>& touched_cities) {
        if (depth == 5) {
            return;
        }
        int next_city = get_next_city(current_city);
        int previous_city = get_previous_city(current_city);
        auto neighbours = boost::adjacent_vertices(current_city, *graph);
        for (auto neighbour : make_iterator_range(neighbours)) {
            if (neighbour == next_city || neighbour == previous_city) {
                continue;
            }
            if (touched_cities.find(neighbour) != touched_cities.end()) {
                continue;
            }
            int neighbour_previous = get_previous_city(neighbour);
            if (touched_cities.find(neighbour_previous) != touched_cities.end()) {
                continue;
            }
            long double profit = calculate_swap_profit(t1, current_city, neighbour_previous, neighbour);
            if (gain + profit > 0) {
                swap_edges(t1, neighbour);
                swaps_done.emplace_back(make_pair(t1, neighbour));
                gain += profit;
                if (gain > max_gain_ptr) {
                    max_gain_ptr = gain;
                    best_swaps_done = swaps_done;
                }
                touched_cities.insert(neighbour);
                touched_cities.insert(neighbour_previous);
                improve_route(t1, neighbour_previous, depth+1, gain, max_gain_ptr, swaps_done, best_swaps_done, touched_cities);
                gain -= profit;
                touched_cities.erase(neighbour_previous);
                touched_cities.erase(neighbour);
                swaps_done.pop_back();
                swap_edges(t1, neighbour);
            }
        }
    }

    void neighbour_profit_swap(int t1) {
        int current_city = get_next_city(t1);
        vector<pair<int, int>> swaps_done;
        vector<pair<int, int>> best_swaps_done;
        long double max_gain = 0;
        unordered_set<int> touched_cities;
        touched_cities.insert(t1);
        improve_route(t1, current_city, 0, 0, max_gain, swaps_done, best_swaps_done, touched_cities);
        for (auto city_swaps : best_swaps_done) {
            swap_edges(city_swaps.first, city_swaps.second);
        }
        length -= max_gain;
    }

    long double calculate_swap_profit(int t1, int t2, int t3, int t4) {
        Point* p1 = get_point(t1);
        Point* p2 = get_point(t2);
        Point* p3 = get_point(t3);
        Point* p4 = get_point(t4);

        return euclidean(p1, p2) + euclidean(p3, p4) - euclidean(p1, p3) - euclidean(p2, p4);
    }

    [[nodiscard]] int get_next_city(int city) const{
        return tree->nextValue(city);
    }

    [[nodiscard]] int get_previous_city(int city) const{
        return tree->previousValue(city);
    }

    [[nodiscard]] long double get_length() const {
        return length;
    }

    long double get_actual_length() {
        long double actual_length = 0;
        auto city_sequence = tree->get_sequence();
        Point* previousCity = cityData[(*city_sequence)[city_sequence->size()-1]];
        for (auto city : *city_sequence) {
            Point* point = cityData[city];
            actual_length += euclidean(previousCity, point);
            previousCity = point;
        }

        return actual_length;
    }

    void print() const {
        tree->print();
    }
    void print_solution() const {
        tree->print_solution();
    }
};