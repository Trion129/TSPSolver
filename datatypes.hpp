#include<iostream>
#include<fstream>
#include<unordered_set>
#include<vector>
#include<stack>
#include<cmath>
#include<algorithm>
#include<random>
#include<chrono>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Delaunay_triangulation_2.h>
#include <boost/unordered_map.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <iostream>

using namespace std;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> Triangulation;
typedef K::Point_2 Point_2;
typedef boost::unordered_map<pair<double, double>, int, boost::hash<pair<double, double>>> pair_to_int;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
        boost::property< boost::vertex_index_t, size_t>,
        boost::property< boost::edge_index_t, size_t, boost::property<boost::edge_weight_t,double> >> Graph;
typedef boost::graph_traits<Triangulation>::vertex_descriptor                         triangulation_vertex_descriptor;
typedef boost::graph_traits<Triangulation>::edge_descriptor                           triangulation_edge_descriptor;
typedef boost::graph_traits<Graph>::vertex_descriptor                                 graph_vertex_descriptor;
typedef boost::graph_traits<Graph>::edge_descriptor                                   graph_edge_descriptor;
typedef boost::unordered_map<triangulation_vertex_descriptor,int>                     triangulation_vertex_to_index_map;
typedef boost::unordered_map<int, triangulation_vertex_descriptor>                    triangulation_index_to_vertex_map;
typedef boost::associative_property_map<triangulation_vertex_to_index_map>            vertex_to_id_associative_map;

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