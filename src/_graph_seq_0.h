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
#include <boost/functional/hash.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "_misc_0.h"

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

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> DirectedGraph;

class SpliceGraph {
public:
	uint graph_id;
	vector<SeqConstraint> read_constraints;
	DirectedGraph graph;
};

}

#endif // _GRAPH_SEQ_0_H_
