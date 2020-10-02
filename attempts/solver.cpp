#include<iostream>
#include<fstream>
#include<unordered_set>
#include<vector>
#include<unordered_map>
#include<stack>
#include<cmath>
#include<algorithm>
#include<random>
#include<chrono>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Delaunay_triangulation_2.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>

#include <fstream>
#include"lru.hpp"

using namespace std;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> Triangulation;
typedef K::Point_2 Point_2;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                    boost::property< boost::vertex_index_t, size_t>,
                    boost::property< boost::edge_index_t, size_t, boost::property<boost::edge_weight_t,double> >> Graph;
typedef boost::graph_traits<Triangulation>::vertex_descriptor       vertex_descriptor;
typedef boost::graph_traits<Triangulation>::vertex_iterator         vertex_iterator;
typedef boost::graph_traits<Triangulation>::edge_descriptor         edge_descriptor;
typedef boost::graph_traits<Graph>::vertex_descriptor               mst_vertex_descriptor;
typedef boost::graph_traits<Graph>::edge_descriptor                 mst_edge_descriptor;
typedef std::map<vertex_descriptor,int>                             VertexIndexMap;
typedef boost::associative_property_map<VertexIndexMap>             VertexIdPropertyMap;

struct hash_pair {
    size_t operator()(const pair<int, int>& p) { 
        auto hash1 = boost::hash<int>{}(p.first); 
        auto hash2 = boost::hash<int>{}(p.second); 
        return hash1 ^ hash2; 
    } 
};

struct Point {
    int city;
    long double x;
    long double y;

    Point(int city, double x, double y): city(city), x(x), y(y) {}
};

struct Node {
    int city;
    Node* before;
    Node* after;

    explicit Node(int city) : city(city), after(NULL), before(NULL) {}

    explicit Node(int city, Node* before, Node* after) :
        city(city), before(before), after(after) {}
};

long double euclidean(Point* p1, Point* p2) {
    return sqrt(pow(p1->x - p2->x, 2) + pow(p1->y - p2->y, 2));
}

struct Route {
    Node* start;
    Node* end;
    vector<pair<Node*, Point*>> cityData;
    long double length;

    explicit Route(int n_cities) : cityData(vector<pair<Node*, Point*>>(n_cities)), start(NULL), end(NULL), length(0) {}
    explicit Route(Route* route) : cityData(vector<pair<Node*, Point*>>(route->cityData.size())), start(NULL), end(NULL), length(0) {
        Node* p = route->start->after;
        while(p != route->start) {
            this->insert(route->get_point(p->city));
            p = p->after;
        }
        this->insert(route->get_point(p->city));
        this->close();
    }

    ~Route() {
        for (int i = 0; i < cityData.size(); i++) {
            delete cityData[i].first;
        }
    }

    void insert(Point* point) {
        auto node = new Node(point->city);
        cityData[point->city] = make_pair(node, point);

        if (end != NULL) {
            length += euclidean(get_point(end->city), point);
            end->after = node;
            node->before = end;
        }
        if (start == NULL) {
            start = node;
        }
        end = node;
    }

    Node* get_node(int city) {
        return cityData[city].first;
    }

    Point* get_point(int city) {
        return cityData[city].second;
    }

    void close() {
        end->after = start;
        start->before = end;
        length += euclidean(get_point(end->city), get_point(start->city));
    }

    void swap_edges(int city1, int city2) {
        Node* city1_node = get_node(city1);
        Node* city2_node = get_node(city2);
        Node* city1_next = get_node(get_next_city(city1));
        Node* city2_next = get_node(get_next_city(city2));
        Node* city1_previous = get_node(get_previous_city(city1));
        Node* city2_previous = get_node(get_previous_city(city2));

        swap(city1_next->before, city2_next->before);
        swap(city1_previous->after, city2_previous->after);
        swap(city1_node->before, city2_node->before);
        swap(city1_node->after, city2_node->after);
    }

    bool neighbour_profit_swap(int city, int offset, long double prob, long double temperature) {
        int other_city = -1;
        Node* p = get_node(city);
        for (int i = 0; i < offset; i++) {
            p = p->after;
        }
        other_city = p->city;
        long double profit = calculate_swap_profit(city, p->city);

        long double probability = exp(profit/temperature);

        if (prob > probability) {
            return false;
        }

        p = get_node(get_next_city(city));
        while (p->city != other_city) {
            Node* current_city = p;
            p = p->after;
            swap(current_city->after, current_city->before);
        }
        length -= profit;

        Node* city_node = get_node(city);
        Node* current_city_node = p;
        Node* city_next_node = get_node(get_next_city(city));
        Node* current_city_next_node = get_node(get_next_city(other_city));

        city_next_node->after = current_city_next_node;
        current_city_next_node->before = city_next_node;
        city_node->after = current_city_node;
        current_city_node->after = current_city_node->before;
        current_city_node->before = city_node;
        return true;
    }

    long double calculate_swap_profit(int city, int current_city) {
        Point* p1 = get_point(city);
        Point* p2 = get_point(get_next_city(city));
        Point* p3 = get_point(current_city);
        Point* p4 = get_point(get_next_city(current_city));

        return euclidean(p1, p2) + euclidean(p3, p4) - euclidean(p1, p3) - euclidean(p2, p4);
    }

    void confirmed_swap_edges(int city1, int city2) {
        length -= calculate_profit(city1, city2);
        swap_edges(city1, city2);
    }

    int get_next_city(int city){
        return cityData[city].first->after->city;
    }

    int get_previous_city(int city){
        return cityData[city].first->before->city;
    }

    long double calculate_profit(int city1, int city2) {
        Point* city1_node = get_point(city1);
        Point* city2_node = get_point(city2);
        Point* city1_next = get_point(get_next_city(city1));
        Point* city2_next = get_point(get_next_city(city2));
        Point* city1_previous = get_point(get_previous_city(city1));
        Point* city2_previous = get_point(get_previous_city(city2));

        if (city1_next == city2_node) {
            return (euclidean(city1_node, city1_previous) - euclidean(city2_node, city1_previous)) +
                    (euclidean(city2_node, city2_next) - euclidean(city1_node, city2_next));
        } else if (city1_previous == city2_node) {
            return (euclidean(city2_node, city2_previous) - euclidean(city1_node, city2_previous)) +
                    (euclidean(city1_node, city1_next) - euclidean(city2_node, city1_next));
        }

        double edge1 = euclidean(city1_node, city1_next);
        double edge2 = euclidean(city2_node, city2_next);
        double edge1_prev = euclidean(city1_node, city1_previous);
        double edge2_prev = euclidean(city2_node, city2_previous);
        double edge1_new = euclidean(city1_node, city2_next);
        double edge2_new = euclidean(city2_node, city1_next);
        double edge1_prev_new = euclidean(city1_node, city2_previous);
        double edge2_prev_new = euclidean(city2_node, city1_previous);

        return (edge1 + edge2 + edge1_prev + edge2_prev) - (edge1_new + edge2_new + edge1_prev_new + edge2_prev_new);
    }

    long double get_length() {
        return length;
    }

    long double calculate_length() {
        double path_length = 0;
        Node* p = start;
        for(int i = 0; i < cityData.size(); i++) {
            Node* next = p->after;
            path_length += euclidean(get_point(p->city), get_point(next->city));
            p = next;
        }
        return path_length;
    }

    void print(unordered_map<pair<double, double>, int, boost::hash<pair<double, double>>> coordToId) {
        Node* ptr = start;
        for(int i = 0; i < cityData.size(); i++) {
            Point* pt = get_point(ptr->city);
            int actualId = coordToId[make_pair(pt->x, pt->y)];

            cout << actualId;
            if (i == cityData.size()){
                cout << "\n";
            } else {
                cout << " ";
            }
            ptr = ptr->after;
        }
    }
    void print_reverse() {
        Node* ptr = start;
        for(int i = 0; i < cityData.size()-1; i++) {
            cout << ptr->city << " ";
            ptr = ptr->before;
        }
        cout << ptr->city << "\n";
    }
};

Route* solve_it_2_swap(Route* route, const int n_cities) {
    int MAX_NEIGHBOURS = min(1000, n_cities-3);
    int MAX_MINIMUM = 3;
    random_device rd;
    default_random_engine rng(rd());
    uniform_int_distribution<int> rnd(0, n_cities-1);
    uniform_int_distribution<int> rnd_minimum(0, MAX_MINIMUM-1);
    uniform_int_distribution<int> rnd_offset(2, MAX_NEIGHBOURS);
    uniform_real_distribution<long double> prob(0, 1);

    int lru_size = n_cities/10;
    cache::lru_cache<size_t, bool> taboo_lru(lru_size);
    vector<Route*> minimum_solutions(MAX_MINIMUM, NULL);
    int current_fill = 0;

    int unchanged_count = 0;
    long double length = route->get_length();
    long double temperature = 100500.0;
    long double decrement = 0.97;

    Route* minimum_route = new Route(route);
    minimum_solutions[current_fill] = new Route(route);
    current_fill++;

    auto start_time = chrono::steady_clock::now();
    while(true) {
        auto current_time = chrono::steady_clock::now();
        int seconds_passed = chrono::duration_cast<std::chrono::seconds>(current_time - start_time).count();

        if (seconds_passed == 600) {
            break;
        }

        int city = rnd(rng);
        int other_city = rnd_offset(rng);
        // while (taboo_lru.exists(city)){
        //     city = rnd(rng);
        // }

        
        bool selected = route->neighbour_profit_swap(city, other_city, prob(rng), temperature);
        if (selected) {
            taboo_lru.put(city, true);
            temperature *= decrement;
        }


        if (route->get_length() < minimum_route->get_length()) {
            minimum_route = new Route(route);
            minimum_solutions[current_fill] = minimum_route;
            current_fill = (current_fill + 1) % MAX_MINIMUM;
        }

        if (route->get_length() - length < 1) {
            unchanged_count++;
        } else {
            unchanged_count = 0;
        }

        length = route->get_length();

        if(unchanged_count > 500) {
            temperature = 100500.0;
            delete route;
            Route* min_solution = NULL;
            while (min_solution == NULL) {
                int choice = rnd_minimum(rng);
                min_solution = minimum_solutions[choice];
            }
            route = new Route(min_solution);
        }
        if (temperature < 1.0e-6) {
            temperature = 100500.0;
            delete route;
            Route* min_solution = NULL;
            while (min_solution == NULL) {
                int choice = rnd_minimum(rng);
                min_solution = minimum_solutions[choice];
            }
            route = new Route(min_solution);
        }
    }

    return minimum_route;
}

class DFSVisitor : public boost::default_dfs_visitor
{
    public:
    boost::shared_ptr<std::vector<int> > point_sequence;

    DFSVisitor() : point_sequence(new std::vector<int>()) {}

    void discover_vertex(mst_vertex_descriptor v, const Graph& g) {
        this->point_sequence->push_back(int(v));
    }
};

Route* apply_delaunay(vector<Point_2> cities_delaunay, int n_cities) {
    vector<Point*> cities(n_cities, NULL);
    Triangulation tr;
    tr.insert(cities_delaunay.begin(), cities_delaunay.end());

    VertexIndexMap vertex_id_map;
    std::map<int, vertex_descriptor> id_vertex_map;
    VertexIdPropertyMap vertex_index_pmap(vertex_id_map);
    int index = 0;
    for(vertex_descriptor vd : vertices(tr)) {
        vertex_id_map[vd] = index;
        id_vertex_map[index] = vd;
        index++;
    }

    // We use the default edge weight which is the squared length of the edge
    // This property map is defined in graph_traits_Triangulation_2.h
    // In the function call you can see a named parameter: vertex_index_map
    std::vector<edge_descriptor> mst;
    boost::kruskal_minimum_spanning_tree(tr, std::back_inserter(mst), vertex_index_map(vertex_index_pmap));
    
    Graph MST(n_cities);
    for (edge_descriptor ed: mst) {
        vertex_descriptor svd = source(ed, tr);
        vertex_descriptor tvd = target(ed, tr);

        mst_vertex_descriptor svd_mst = boost::vertex(vertex_id_map[svd], MST);
        mst_vertex_descriptor tvd_mst = boost::vertex(vertex_id_map[tvd], MST);
        boost::add_edge(svd_mst, tvd_mst, MST);
    }

    mst_vertex_descriptor start = boost::vertex(0, MST);
    DFSVisitor vis;
    boost::depth_first_search(MST, boost::visitor(vis).root_vertex(start));

    vector<int> final_seq = *vis.point_sequence;
    Route* start_route = new Route(n_cities);
    for (int i = 0; i < final_seq.size(); i++) {
        Point_2 pt = id_vertex_map[final_seq[i]]->point();
        cities[i] = new Point(i, pt.x(), pt.y());
        start_route->insert(cities[i]);
    }
    start_route->close();

    return start_route;
}

int main() {
    int n_cities;
    cin >> n_cities;
    vector<Point_2> cities_delaunay(n_cities);
    unordered_map<pair<double, double>, int, boost::hash<pair<double, double>>> coord_to_id;
    hash_pair hasher;
    for (int i = 0; i < n_cities; i++) {
        double x, y;

        cin >> x >> y;

        coord_to_id[make_pair(x, y)] = i;
        cities_delaunay[i] = Point_2(x, y);
    }

    Route* start_route = apply_delaunay(cities_delaunay, n_cities);

    int RESTARTS = 10;
    Route* best_solution = start_route;
    vector<Route*> solutions(RESTARTS);
    #pragma omp parallel for
    for(int i = 0; i < RESTARTS; i++) {
        // Route* solution = solve_it(n_cities, cities);
        Route* solution = new Route(start_route);
        solutions[i] = solve_it_2_swap(solution, n_cities);
    }

    for(int i = 0; i < RESTARTS; i++) {
        Route* solution = solutions[i];
        if(best_solution == NULL || solution->get_length() < best_solution->get_length()) {
            best_solution = solution;
        }
    }
    long double total_distance = best_solution->get_length();

    cout << total_distance << " " << "0" << "\n";
    best_solution->print(coord_to_id);

    return 0;
}
