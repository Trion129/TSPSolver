#include "helper.hpp"

using namespace std;

int main() {
    int n_cities;
    cin >> n_cities;
    vector<Point_2> cities_delaunay(n_cities);
    pair_to_int coord_to_id;
    for (int i = 0; i < n_cities; i++) {
        double x, y;

        cin >> x >> y;

        coord_to_id[make_pair(x, y)] = i;
        cities_delaunay[i] = Point_2(x, y);
    }

    auto route = apply_delaunay(cities_delaunay, n_cities, coord_to_id);

    long double currentLength = route->get_length();
    long double prevLength = 0;
    while (currentLength != prevLength) {
        for (int i = 0; i < n_cities; i++) {
            route->neighbour_profit_swap(i);
        }
        prevLength = currentLength;
        currentLength = route->get_length();
    }
    cout << route->get_actual_length() << " 0\n";
    route->print();

    return 0;
}

