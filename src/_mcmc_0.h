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

#include <cstddef>
#include <vector>
#include <utility>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/tuple/tuple.hpp>
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
using _graph_seq_0::DGVertex;
using _graph_seq_0::SeqConstraint;
using boost::tie;
using boost::in_edges;
using boost::in_degree;
using boost::source;
using boost::out_degree;
using boost::out_edges;
using boost::target;

typedef _graph_seq_0::PyGraph GraphInfo;

typedef boost::property_map<DirectedGraph, boost::vertex_index_t>::type DGIndexMap;
typedef boost::graph_traits<DirectedGraph>::vertex_iterator DGVertexIter;
typedef boost::graph_traits<DirectedGraph>::adjacency_iterator DGAdjIter;
typedef boost::graph_traits<DirectedGraph>::in_edge_iterator DGInEdgeIter;
typedef boost::graph_traits<DirectedGraph>::out_edge_iterator DGOutEdgeIter;
typedef boost::graph_traits<DirectedGraph>::edge_descriptor DGEdge;

// calculate the probability that a read
// is from a the transcripts (isoforms)
template<class RNodeLoc>
class ReadFromTransProb {
	inline ldbl operator()(ReadInGraph<RNodeLoc> const &r) {
		return 0.0;
	}
};

class GSLRngUnifInt {
public:
	GSLRngUnifInt(gsl_rng *rn_) :
			rn(rn_) {
	}

	uint operator()(uint n) {
		return gsl_rng_uniform_int(rn, n);
	}

private:
	gsl_rng *rn;
};

// given read constraint
// and graph, get a random isoform
// that satisfies the constraint
inline void rand_rc_isof(SpliceGraph const &graph, Isoform &isof, uint un_rc,
		gsl_rng *rn) {
	// size of @isof is already set
	// read constraint asumes all 0's doesn't occur

	SeqConstraint const &rc = graph.read_constraints[un_rc];

	// ts - topological sort
	vector<DGVertex>::const_iterator ts_iter = graph.topo_sort.begin();

	while (rc[*ts_iter] == false) {
		++ts_iter;
	}

	{
		// go "up" to find transcript start head

		uint cur_node = *ts_iter;
		isof.set(cur_node);

		DGInEdgeIter in_i, in_end;
		tie(in_i, in_end) = in_edges(cur_node, graph.graph);

		while (in_i != in_end) {
			uint go2node = gsl_rng_uniform_int(rn,
					in_degree(cur_node, graph.graph));

			for (uint i = 0; i != go2node; ++i) {
				++in_i;
			}

			cur_node = source(*in_i, graph.graph);
			isof.set(cur_node);
			tie(in_i, in_end) = in_edges(cur_node, graph.graph);
		}
	}

	// last node that has a constraint
	uint last_rc_node = *ts_iter;
	++ts_iter;

	while (ts_iter != graph.topo_sort.end()) {
		if (rc[*ts_iter] == true) {
			uint cur_node = *ts_iter;

			isof.set(cur_node);

			vector<DGEdge> edges_in;

			DGInEdgeIter in_i, in_end;
			tie(in_i, in_end) = in_edges(cur_node, graph.graph);

			while (in_i != in_end) {
				edges_in.push_back(*in_i);
				++in_i;
			}

			GSLRngUnifInt rn_wrapper(rn);
			std::random_shuffle(edges_in.begin(), edges_in.end(), rn_wrapper);



			last_rc_node = cur_node;
		}
		++ts_iter;
	}

	{
		// extend last node to end,
		// go "down" to end

		uint cur_node = last_rc_node;

		isof.set(cur_node);

		DGOutEdgeIter out_i, out_end;
		tie(out_i, out_end) = out_edges(cur_node, graph.graph);

		while (out_i != out_end) {
			uint go2node = gsl_rng_uniform_int(rn,
					out_degree(cur_node, graph.graph));

			for (uint i = 0; i != go2node; ++i) {
				++out_i;
			}

			cur_node = target(*out_i, graph.graph);
			isof.set(cur_node);
			tie(out_i, out_end) = out_edges(cur_node, graph.graph);
		}
	}

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
			tie(in_i, in_end) = in_edges(*j.first, i->graph);
			if (in_i == in_end) {
				sn_iter->push_back(*j.first);
			}
		}
	}

	// TODO: multiple chains pthread parallelization

	{

		gsl_rng *rn = gsl_rng_alloc(gsl_rng_mt19937);

		// TODO: seed @rn

		vector<IsoformSet> isoforms(graphs.size());
		vector<IsoformSet>::iterator isof_set_iter = isoforms.begin();
		sn_iter = start_nodes.begin();

		for (vector<SpliceGraph>::const_iterator i = graphs.begin();
				i != graphs.end(); ++i, ++isof_set_iter, ++sn_iter) {

			uint rc_size = i->read_constraints.size();

			// 1 - unsatisfied
			// 0 - satisfied
			dynamic_bitset<> satisfied_rc(rc_size);
			satisfied_rc.set();

			while (satisfied_rc.any()) {
				uint un_rc = 0; // un-satisfied read constraint index
				while (satisfied_rc[un_rc] == false) {
					++un_rc;
				}
				satisfied_rc[un_rc] = false;

				Isoform isof(boost::num_vertices(i->graph));
				rand_rc_isof(*i, isof, un_rc, rn);
				isof_set_iter->insert(isof);

				for (uint j = 0; j != rc_size; ++j) {
					if (satisfied_rc[j] == true) {
						if (isof == (isof | i->read_constraints[j])) {
							satisfied_rc[j] = false;
						}
					}
				}
			}

			// TODO: remove this
			for (IsoformSet::const_iterator j = isof_set_iter->begin();
					j != isof_set_iter->end(); ++j) {
				std::cout << *j << std::endl;
			}

		}

		uint runs = 0;

		while (runs != max_run) {

			++runs;
		}

		gsl_rng_free(rn);
	}

}

}

#endif // _MCMC_0_H_
