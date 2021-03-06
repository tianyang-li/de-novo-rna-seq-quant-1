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

#include <vector>
#include <exception>
#include <iostream>
#include <boost/unordered_set.hpp>

#include "_single_1.h"
#include "_misc_0.h"
#include "_graph_seq_0.h"
#include "_mcmc_0.h"

namespace _single_1 {

using std::vector;
using _graph_seq_0::SeqConstraint;
using _graph_seq_0::ReadIndex;
using std::cerr;
using std::endl;
using std::exception;

}

namespace _single_1 {

// setup graph read constraints
inline void setup_graph_rc(_graph_seq_0::SpliceGraph &graph,
		_graph_seq_0::PyGraph const &py_graph, GraphReads const &graph_read,
		vector<ReadInGraph<SingleNodeLoc> > &py_reads) {

	boost::unordered_set<SeqConstraint, _graph_seq_0::SeqConstraintHash> rcs_set;

	for (vector<ReadIndex>::const_iterator i = graph_read.reads.begin();
			i != graph_read.reads.end(); ++i) {
		SeqConstraint rc(py_graph.nodes.size());

		for (_graph_seq_0::ReadNodeLoc::const_iterator j =
				py_reads[i->read_id].graph_locs[i->graph_index].locs[i->align_index].begin();
				j
						!= py_reads[i->read_id].graph_locs[i->graph_index].locs[i->align_index].end();
				++j) {
			rc.set(j->node_id);
		}

		rcs_set.insert(rc);
	}

	// read constraints
	vector<SeqConstraint> rcs_vec;

	for (boost::unordered_set<SeqConstraint, _graph_seq_0::SeqConstraintHash>::const_iterator i =
			rcs_set.begin(); i != rcs_set.end(); ++i) {
		rcs_vec.push_back(*i);
	}

	// build a partial order on read constraints
	// count in-degree
	vector<ulong> rc_in(rcs_vec.size(), 0);

	for (ulong i = 1; i < rcs_vec.size(); ++i) {
		for (ulong j = 0; j < i; ++j) {
			SeqConstraint ij_or = rcs_vec[i] | rcs_vec[j];
			if (ij_or == rcs_vec[i]) {
				rc_in[j] += 1;
			} else {
				if (ij_or == rcs_vec[j]) {
					rc_in[i] += 1;
				}
			}
		}
	}

	ulong rc_id = 0;

	for (vector<ulong>::const_iterator i = rc_in.begin(); i != rc_in.end();
			++i, ++rc_id) {
		if (*i == 0) {
			graph.read_constraints.push_back(rcs_vec[rc_id]);
		}
	}
}

inline void setup_graph(_graph_seq_0::SpliceGraph &graph,
		_graph_seq_0::PyGraph &py_graph, GraphReads const &graph_read,
		vector<ReadInGraph<SingleNodeLoc> > &py_reads) {
	// setup boost graph
	// setup read constraints

	// TODO: fix this for final version
	// here no non-branch start/end site detection
	// or trimming is done
	// transcript length is also not estimated

	for (ulong i = 0; i != py_graph.nodes.size(); ++i) {
		boost::add_vertex(graph.graph);
	}
	for (vector<_graph_seq_0::Node>::const_iterator i = py_graph.nodes.begin();
			i != py_graph.nodes.end(); ++i) {
		for (vector<ulong>::const_iterator j = i->edges.begin();
				j != i->edges.end(); ++j) {
			boost::add_edge(i->node_id, *j, graph.graph);
		}
	}

	// TODO: fix graph node when segmentation is done

	graph.setup();

	setup_graph_rc(graph, py_graph, graph_read, py_reads);

}

void _get_isoforms(vector<_graph_seq_0::PyGraph> *py_graphs,
		vector<ReadInGraph<SingleNodeLoc> > *py_reads,
		vector<_graph_seq_0::Fasta> *isoforms, ulong max_run) {

#ifdef DEBUG
	cerr << "enter extension" << endl;
#endif

	try {

		// the reads that a graph has
		vector<GraphReads> graph_reads(py_graphs->size());

		ulong graph_id = 0;

		for (vector<GraphReads>::iterator i = graph_reads.begin();
				i != graph_reads.end(); ++i, ++graph_id) {
			i->graph_id = graph_id;
		}

		ulong read_id = 0;

		for (vector<ReadInGraph<SingleNodeLoc> >::const_iterator i =
				py_reads->begin(); i != py_reads->end(); ++i, ++read_id) {
			ulong graph_index = 0;

			for (vector<ReadGraphLoc<SingleNodeLoc> >::const_iterator j =
					i->graph_locs.begin(); j != i->graph_locs.end();
					++j, ++graph_index) {
				ulong align_index = 0;

				for (vector<_graph_seq_0::ReadNodeLoc>::const_iterator k =
						j->locs.begin(); k != j->locs.end();
						++k, ++align_index) {
					graph_reads[j->graph_id].reads.push_back(
							ReadIndex(read_id, graph_index, align_index));
				}
			}
		}

		vector<_graph_seq_0::SpliceGraph> graphs(py_graphs->size());

		// setup graphs

		graph_id = 0;
		vector<GraphReads>::const_iterator graph_read = graph_reads.begin();
		vector<_graph_seq_0::PyGraph>::iterator py_graph = py_graphs->begin();

		for (vector<_graph_seq_0::SpliceGraph>::iterator i = graphs.begin();
				i != graphs.end(); ++i, ++graph_id, ++graph_read, ++py_graph) {
			i->graph_id = graph_id;

			setup_graph(*i, *py_graph, *graph_read, *py_reads);

		}

		// for graph information:
		//     @py_graphs sequences, seq_len, est_len -> @graph_infos
		// 	   @graphs structure and isoforms

		// for read to graph alignment information:
		//     @py_reads would have been modified and -> @read_in_graph
		//     will continue to be used

		// for each graph's alignments information:
		//     @graph_reads will continue to be used

		try {
			_mcmc_0::isoform_main<SingleNodeLoc>(*py_graphs, graphs, *py_reads,
					graph_reads, max_run);
		} catch (exception &e) {
			cerr << e.what() << endl;
			throw;
		}

		try {
			_graph_seq_0::get_isoform_FASTA(graphs, *py_graphs, *isoforms);
		} catch (exception &e) {
			cerr << e.what() << endl;
			throw;
		}

	} catch (exception &e) {
		cerr << e.what() << endl;
		throw;
	}

#ifdef DEBUG
	cerr << "exit extension" << endl;
#endif

}

} // namespace _single_1

