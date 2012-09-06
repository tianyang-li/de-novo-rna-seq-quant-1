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

#ifndef _MCMC_0_H_
#define _MCMC_0_H_

#include <vector>
#include <utility>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/dynamic_bitset.hpp>
#include <gsl/gsl_rng.h>

#include "_graph_seq_0.h"
#include "_misc_0.h"

// TODO: remove this
#include <iostream>

namespace _mcmc_0 {

using std::pair;
using _graph_seq_0::SpliceGraph;
using _graph_seq_0::GraphReads;
using std::vector;
using _graph_seq_0::ReadInGraph;
using _graph_seq_0::check_isoform_rc;
using _graph_seq_0::IsoformSet;
using _graph_seq_0::Isoform;
using _graph_seq_0::DirectedGraph;
using boost::dynamic_bitset;

typedef _graph_seq_0::PyGraph GraphInfo;

typedef boost::property_map<DirectedGraph, boost::vertex_index_t>::type DGIndexMap;
typedef boost::graph_traits<DirectedGraph>::vertex_iterator DGVertexIter;
typedef boost::graph_traits<DirectedGraph>::adjacency_iterator DGAdjIter;
typedef boost::graph_traits<DirectedGraph>::in_edge_iterator DGInEdgeIter;

// calculate the probability that a read
// is from a the transcripts (isoforms)
template<class RNodeLoc>
class ReadFromTransProb {
	inline ldbl operator()(ReadInGraph<RNodeLoc> const &r) {
		return 0.0;
	}
};

// given read constraint
// and graph, get a random isoform
// that satisfies the constraint
inline void rand_rc_isof(SpliceGraph const &graph, Isoform &isof, gsl_rng *rn) {
	// size of @isof is already set
}

template<class RNodeLoc>
void isoform_main(vector<GraphInfo> const &graph_info,
		vector<SpliceGraph> &graphs,
		vector<ReadInGraph<RNodeLoc> > const &read_in_graph,
		vector<GraphReads> const &graph_reads,
		ReadFromTransProb<RNodeLoc> &r_prob, uint max_run) {

	// all the nodes with in-dgree == 0
	vector<vector<uint> > start_nodes(graphs.size());
	vector<vector<uint> >::iterator sn_iter = start_nodes.begin();

	for (vector<SpliceGraph>::const_iterator i = graphs.begin();
			i != graphs.end(); ++i, ++sn_iter) {
		for (pair<DGVertexIter, DGVertexIter> j = boost::vertices(i->graph);
				j.first != j.second; ++j.first) {
			DGInEdgeIter in_i, in_end;
			boost::tie(in_i, in_end) = boost::in_edges(i->index[*j.first],
					i->graph);
			if (in_i == in_end) {
				sn_iter->push_back(i->index[*j.first]);
			}
		}
	}

	// TODO: multiple chains pthread parallelization

	gsl_rng *rn = gsl_rng_alloc(gsl_rng_mt19937);

	vector<IsoformSet> isoforms(graphs.size());
	vector<IsoformSet>::iterator isof_iter = isoforms.begin();
	sn_iter = start_nodes.begin();

	for (vector<SpliceGraph>::const_iterator i = graphs.begin();
			i != graphs.end(); ++i, ++isof_iter, ++sn_iter) {
		dynamic_bitset<> satisfied_rc(i->read_constraints.size());
		satisfied_rc.set();
		while (satisfied_rc.any()) {
			uint un_rc = 0; // un-satisfied read constraint index
			while (satisfied_rc[un_rc] == false) {
				++un_rc;
			}
			satisfied_rc[un_rc] = false;

			Isoform isof(boost::num_vertices(i->graph));
			rand_rc_isof(*i, isof, rn);
		}

	}

	uint runs = 0;

	while (runs != max_run) {

		++runs;
	}

	gsl_rng_free(rn);

}

}

#endif // _MCMC_0_H_
