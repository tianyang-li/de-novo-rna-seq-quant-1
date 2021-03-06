/*
 * Tianyang Li 2012
 * tmy1018@gmail.com
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * some general stuff
 * for working with sequences, graphs, etc.
 */

#ifndef _GRAPH_SEQ_0_H_
#define _GRAPH_SEQ_0_H_

#include <vector>
#include <string>
#include <cstddef>
#include <iterator>
#include <algorithm>
#include <utility>
#include <iostream>
#include <exception>
#include <boost/functional/hash.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/unordered_set.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/transitive_closure.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/unordered_map.hpp>

#include "_misc_0.h"
#include "_err_0.h"

namespace _graph_seq_0 {

using std::vector;
using std::string;
using boost::num_vertices;
using boost::out_edges;
using std::cerr;
using std::endl;
using boost::print_graph;
using boost::in_degree;
using boost::tie;
using boost::out_edges;
using std::exception;

#ifdef DEBUG
using _err_0::IteratorEndError;
#endif

}

namespace _graph_seq_0 {

/*
 * makes it easy to get graph information from
 * python
 */

/*
 * assumes that we know what what graph it is
 */
class Node {
public:
	Node(ulong node_id_, string seq_) :
			node_id(node_id_), seq(seq_), est_len(seq_.size()) {
	}

	ulong node_id;

	// all the alignment positions are according to
	// @seq, but the length used in estimation is
	// @est_len
	string seq;
	ulong est_len;

	vector<ulong> edges;
};

// graph info from python (act as glue)
class PyGraph {
public:
	PyGraph() {
	}

	ulong graph_id;
	vector<Node> nodes;

	inline ulong get_gene_len() const {
		ulong gene_len = 0;

		for (vector<Node>::const_iterator i = nodes.begin(); i != nodes.end();
				++i) {
			gene_len += i->est_len;
		}

		return gene_len;
	}
};

typedef PyGraph GraphInfo;

/*
 * assumes that we already know which graph it is
 */
class SeqLoc {
public:
	SeqLoc(ulong node_id_, ulong start_, ulong end_) :
			node_id(node_id_), start(start_), end(end_) {
	}

	ulong node_id;
	ulong start;
	ulong end;
};

/*
 * assumes the graph is already known
 *
 * this contains all the nodes that this read is aligned to
 * in the graph
 */
typedef vector<SeqLoc> ReadNodeLoc;

class Fasta {
public:
	string info;
	string seq; // only ACGT allowed
};

// size if the length of the graph
// bit set to 1 if the seq falls into
// the corresponding node
typedef boost::dynamic_bitset<> SeqConstraint;

// TODO: ways to optimize?
class SeqConstraintHash {
public:
	inline size_t operator()(SeqConstraint const &x) const {
		string x_str;
		boost::to_string(x, x_str);
		boost::hash<string> hasher;
		return hasher(x_str);
	}
};

typedef boost::adjacency_list<boost::hash_setS, boost::vecS,
		boost::bidirectionalS> DirectedGraph;

typedef SeqConstraintHash IsoformHash;
typedef boost::dynamic_bitset<> Isoform;
typedef boost::unordered_set<Isoform, IsoformHash> IsoformSet;
typedef boost::graph_traits<DirectedGraph>::vertex_descriptor DGVertex;

typedef boost::graph_traits<DirectedGraph>::vertex_iterator DGVertexIter;
typedef boost::graph_traits<DirectedGraph>::adjacency_iterator DGAdjIter;
typedef boost::graph_traits<DirectedGraph>::in_edge_iterator DGInEdgeIter;
typedef boost::graph_traits<DirectedGraph>::out_edge_iterator DGOutEdgeIter;
typedef boost::graph_traits<DirectedGraph>::edge_descriptor DGEdge;

// stores isoform and its corresponding expression level
typedef boost::unordered_map<Isoform, double, IsoformHash> IsoformMap;

// graph used in actual calculation
class SpliceGraph {
public:

	// for @dist_from_starts
	// the "dist_from_start" for a source
	// vertex
	static double const kSourceDistFromStart = 1.0;

	// for @dist_from_starts
	// distance between 2 vertices
	static double const kVertDist = 1.0;

	ulong graph_id;

	vector<SeqConstraint> read_constraints;

	DirectedGraph graph;
	DirectedGraph tc; // transitive closure
	vector<DGVertex> topo_sort; // topological sort results

	// TODO: get stuff out of MCMC results
	vector<IsoformMap> mcmc_results;

	// "distance" of each node from a
	// source vertex, can be a mixture of
	// distances
	vector<double> dist_from_starts;

	typedef ulong IsofCount; // XXX: is ulong too small?
	// number of isoforms that can start from a node
	vector<IsofCount> vert_isof_counts;

	inline void setup() {
		try {
			boost::transitive_closure(graph, tc);

			boost::topological_sort(graph, std::back_inserter(topo_sort));
			std::reverse(topo_sort.begin(), topo_sort.end());

			ulong vert_num = num_vertices(graph);

			dist_from_starts.assign(vert_num, 0);

			for (vector<DGVertex>::const_iterator i = topo_sort.begin();
					i != topo_sort.end(); ++i) {

				if (in_degree(*i, graph) != 0) {
					break;
				}

				_set_dist_from_starts(*i, kSourceDistFromStart);

			}

			for (ulong i = 0; i != num_vertices(graph); ++i) {

				if (in_degree(i, graph) > 1) {
					dist_from_starts[i] /= double(in_degree(i, graph));
				}

			}

			_count_vert_isof();

#ifdef DEBUG
			cerr << "########\n";
			cerr << "distances from starts\n";
			for (ulong i = 0; i != dist_from_starts.size(); ++i) {
				cerr << i << ":" << dist_from_starts[i] << endl;
			}
			cerr << "######\n";

			cerr << "graph id: " << graph_id << endl;
			print_graph(graph);

			cerr << "########\n";
			cerr << "topological sort\n";
			for (vector<ulong>::const_iterator i = topo_sort.begin();
					i != topo_sort.end(); ++i) {
				cerr << *i << "\t";
			}
			cerr << endl;
			cerr << "########\n";

			cerr << "number of isoforms starting from a vertex\n";
			for (ulong i = 0; i != vert_isof_counts.size(); ++i) {
				cerr << i << ":" << vert_isof_counts[i] << endl;
			}

			cerr << "########\n";
#endif

		} catch (exception &e) {
			cerr << e.what() << endl;
			throw;
		}

	}

private:
	inline void _count_vert_isof() {
		// count vert_isof_counts, similar to Floyd-Warshall

		ulong vert_num = num_vertices(graph);

		vert_isof_counts.assign(vert_num, 0);

		for (vector<ulong>::const_reverse_iterator cur_start_iter =
				topo_sort.rbegin(); cur_start_iter != topo_sort.rend();
				++cur_start_iter) {

			vert_isof_counts[*cur_start_iter] += 1;

			DGOutEdgeIter out_i, out_end;

			for (tie(out_i, out_end) = boost::out_edges(*cur_start_iter, graph);
					out_i != out_end; ++out_i) {

				// edge: *cur_start_iter -> targ
				ulong targ = boost::target(*out_i, graph);

				vert_isof_counts[*cur_start_iter] += vert_isof_counts[targ];
			}

		}

	}

	inline void _set_dist_from_starts(ulong cur_vert,
			double dfs /* dist from start*/) {

		if (dist_from_starts[cur_vert] == 0) {

			dist_from_starts[cur_vert] = dfs;

			DGOutEdgeIter out_i, out_end;
			tie(out_i, out_end) = out_edges(cur_vert, graph);

			while (out_i != out_end) {
				_set_dist_from_starts(target(*out_i, graph), dfs + kVertDist);
				++out_i;
			}

		} else {
			// will be divided by in_degree later
			dist_from_starts[cur_vert] += dfs;
		}

	}

};

/*
 * location of the read
 * in the graph
 *
 * assumes that the read is already known
 */
template<class RNodeLoc> // Read Node Loc
class ReadGraphLoc {
public:
	ReadGraphLoc(ulong graph_id_) :
			graph_id(graph_id_) {
	}

	ulong graph_id;

	// distinct location(s) of the read in the graph
	vector<RNodeLoc> locs;
};

// for a given read, has all the graphs that
// this read is in
template<class RNodeLoc>
class ReadInGraph {
public:
	ReadInGraph() {
	}

	ulong read_id;

	// graphs that contain the read
	vector<ReadGraphLoc<RNodeLoc> > graph_locs;
};

/*
 * given read's @read_id
 * the alignment is at
 * py_reads[read_id].graph_locs[graph_index].locs[align_index]
 */
class ReadIndex {
public:
	ReadIndex(ulong read_id_, ulong graph_index_, ulong align_index_) :
			read_id(read_id_), graph_index(graph_index_), align_index(
					align_index_) {
	}

	ulong read_id;
	ulong graph_index;
	ulong align_index;

	template<class RNodeLoc>
	inline ulong get_align_num(
			vector<ReadInGraph<RNodeLoc> > const &read_in_graph) const {
		ulong align_num = 0;

		ReadInGraph<RNodeLoc> const &cur_read_in_graph = read_in_graph[read_id];

		for (typename vector<ReadGraphLoc<RNodeLoc> >::const_iterator i =
				cur_read_in_graph.graph_locs.begin();
				i != cur_read_in_graph.graph_locs.end(); ++i) {

			for (typename vector<RNodeLoc>::const_iterator j = i->locs.begin();
					j != i->locs.end(); ++j) {
				++align_num;
			}

		}

		return align_num;
	}
};

/*
 * keep the reads aligned to a particular graph
 */
class GraphReads {
public:
	ulong graph_id;
	vector<ReadIndex> reads;
};

void get_isoform_FASTA(vector<SpliceGraph> const &graphs,
		vector<PyGraph> const &py_graphs, vector<Fasta> &isoforms);

inline ulong get_isof_len(Isoform const &isof, GraphInfo const &graph_info) {

	ulong isof_len = 0;

	for (ulong i = 0; i != graph_info.nodes.size(); ++i) {

		if (isof[i] == true) {
			isof_len += graph_info.nodes[i].est_len;
		}

	}

	return isof_len;

}

} // namespace _graph_seq_0

#endif // _GRAPH_SEQ_0_H_
