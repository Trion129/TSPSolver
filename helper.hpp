#include "route.hpp"

class DFSVisitor : public boost::default_dfs_visitor
{
    public:
    boost::shared_ptr<std::vector<int> > point_sequence;

    DFSVisitor() : point_sequence(new std::vector<int>()) {}

    void discover_vertex(graph_vertex_descriptor v, const Graph& g) const {
        this->point_sequence->push_back(int(v));
    }
};

Route* apply_delaunay(vector<Point_2> cities_delaunay, int n_cities, pair_to_int& coord_to_index) {
    Triangulation tr;
    tr.insert(cities_delaunay.begin(), cities_delaunay.end());

    triangulation_vertex_to_index_map vertex_id_map;
    triangulation_index_to_vertex_map id_vertex_map;
    for(const triangulation_vertex_descriptor& vd : vertices(tr)) {
        auto point_pair = make_pair(vd->point().x(), vd->point().y());
        int index = coord_to_index[point_pair];
        vertex_id_map[vd] = index;
        id_vertex_map[index] = vd;
    }

    auto* delaunay_graph = new Graph();
    boost::graph_traits<Triangulation>::edge_iterator ei, eiend;
    for(boost::tie(ei, eiend) = edges(tr); ei != eiend; ++ei) {
        triangulation_edge_descriptor edge_def = *ei;
        triangulation_vertex_descriptor svd = source(edge_def, tr);
        triangulation_vertex_descriptor tvd = target(edge_def, tr);
        graph_vertex_descriptor svd_del = boost::vertex(vertex_id_map[svd], *delaunay_graph);
        graph_vertex_descriptor tvd_del = boost::vertex(vertex_id_map[tvd], *delaunay_graph);
        boost::add_edge(svd_del, tvd_del, *delaunay_graph);
    }

    std::vector<triangulation_edge_descriptor> mst;
    vertex_to_id_associative_map vertex_index_pmap(vertex_id_map);
    boost::kruskal_minimum_spanning_tree(tr, std::back_inserter(mst), vertex_index_map(vertex_index_pmap));
    
    Graph MST(n_cities);
    for (const triangulation_edge_descriptor& ed: mst) {
        triangulation_vertex_descriptor svd = source(ed, tr);
        triangulation_vertex_descriptor tvd = target(ed, tr);

        graph_vertex_descriptor svd_mst = boost::vertex(vertex_id_map[svd], MST);
        graph_vertex_descriptor tvd_mst = boost::vertex(vertex_id_map[tvd], MST);
        boost::add_edge(svd_mst, tvd_mst, MST);
    }

    graph_vertex_descriptor start = boost::vertex(0, MST);
    DFSVisitor vis;
    boost::depth_first_search(MST, boost::visitor(vis).root_vertex(start));

    vector<int> final_seq = *vis.point_sequence;
    vector<Point*> points(n_cities);
    for (int i = 0; i < n_cities; i++) {
        Point_2 pt = id_vertex_map[i]->point();
        points[i] = new Point(i, pt.x(), pt.y());
    }
    return new Route(points, final_seq, delaunay_graph);
}