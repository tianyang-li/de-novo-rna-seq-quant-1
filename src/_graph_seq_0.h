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
#include <boost/functional/hash.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/unordered_set.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/transitive_closure.hpp>
#include <boost/graph/topological_sort.hpp>

#include "_misc_0.h"

// TODO: remove this
#include <iostream>
#include <boost/graph/graph_utility.hpp>

namespace _graph_seq_0 {

using std::vector;
using std::string;

/*
 * makes it easy to get graph information from
 * python
 */

/*
 * assumes that we know what what graph it is
 */
class Node {
public:
	Node(uint node_id_, string seq_) :
			node_id(node_id_), seq(seq_), est_len(seq_.size()) {
	}

	uint node_id;

	// all the alignment positions are according to
	// @seq, but the length used in estimation is
	// @est_len
	string seq;
	uint est_len;

	vector<uint> edges;
};

class PyGraph {
public:
	PyGraph() {
	}

	uint graph_id;
	vector<Node> nodes;
};

/*
 * assumes that we already know which graph it is
 */
class SeqLoc {
public:
	SeqLoc(uint node_id_, uint start_, uint end_) :
			node_id(node_id_), start(start_), end(end_) {
	}

	uint node_id;
	uint start;
	uint end;
};

/*
 * assumes the graph is already known
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
typedef boost::property_map<DirectedGraph, boost::vertex_index_t>::type DGIndexMap;
typedef boost::graph_traits<DirectedGraph>::vertex_descriptor DGVertex;

class SpliceGraph {
public:
	uint graph_id;
	vector<SeqConstraint> read_constraints;
	DirectedGraph graph;
	IsoformSet isoforms;
	DirectedGraph tc; // transitive closure
	vector<DGVertex> topo_sort; // topological sort results

	inline void setup() {
		boost::transitive_closure(graph, tc);
		boost::topological_sort(graph, std::back_inserter(topo_sort));
		std::reverse(topo_sort.begin(), topo_sort.end());
	}
};

/*
 * assumes that the read is already known
 */
template<class RNodeLoc> // Read Node Loc
class ReadGraphLoc {
public:
	ReadGraphLoc(uint graph_id_) :
			graph_id(graph_id_) {
	}

	uint graph_id;
	std::vector<RNodeLoc> locs;
};

template<class RNodeLoc>
class ReadInGraph {
public:
	ReadInGraph() {
	}

	uint read_id;
	std::vector<ReadGraphLoc<RNodeLoc> > graph_locs;
};

/*
 * given read's @read_id
 * the alignment is at
 * py_reads[read_id].graph_locs[graph_index].locs[align_index]
 */
class ReadIndex {
public:
	ReadIndex(uint read_id_, uint graph_index_, uint align_index_) :
			read_id(read_id_), graph_index(graph_index_), align_index(
					align_index_) {
	}

	uint read_id;
	uint graph_index;
	uint align_index;
};

/*
 * keep the reads aligned to a particular graph
 */
class GraphReads {
public:
	uint graph_id;
	std::vector<ReadIndex> reads;
};

void get_isoform_FASTA(vector<SpliceGraph> const &graphs,
		vector<PyGraph> const &py_graphs, vector<Fasta> &isoforms);

// check isoform constraint
inline bool check_isoform_rc(Isoform const &isoform,
		SeqConstraint const &rc /* read constraints */) {
	return isoform == (isoform | rc);
}

}

#endif // _GRAPH_SEQ_0_H_
