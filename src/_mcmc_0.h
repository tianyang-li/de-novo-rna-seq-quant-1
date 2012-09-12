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
#include <algorithm>
#include <memory>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/unordered_map.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

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
using boost::edge;
using _graph_seq_0::IsoformMap;
using std::make_pair;
using boost::unordered_map;
using _graph_seq_0::DGVertexIter;
using _graph_seq_0::DGAdjIter;
using _graph_seq_0::DGInEdgeIter;
using _graph_seq_0::DGOutEdgeIter;
using _graph_seq_0::DGEdge;
using _graph_seq_0::GraphInfo;
using std::auto_ptr;

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

	// extend "back up" to
	// connect to the previous contraint
	while (ts_iter != graph.topo_sort.end()) {
		if (rc[*ts_iter] == true) {
			uint cur_node = *ts_iter;

			isof.set(cur_node);

			while (cur_node != last_rc_node) {

				vector<uint> edges_in;

				DGInEdgeIter in_i, in_end;
				tie(in_i, in_end) = in_edges(cur_node, graph.graph);

				while (in_i != in_end) {
					edges_in.push_back(source(*in_i, graph.graph));
					++in_i;
				}

				GSLRngUnifInt rn_wrapper(rn);
				std::random_shuffle(edges_in.begin(), edges_in.end(),
						rn_wrapper);

				for (vector<uint>::const_iterator i = edges_in.begin();
						i != edges_in.end(); ++i) {
					if (edge(last_rc_node, *i, graph.tc).second
							|| *i == last_rc_node) {
						cur_node = *i;
						isof.set(cur_node);
						break;
					}
				}

			}

			last_rc_node = *ts_iter;
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

inline bool _isof_start_ok(DirectedGraph const &graph, uint cur_vert,
		IsoformMap const &isofs, Isoform &isof) {
	isof.set(cur_vert);

	if (isofs.find(isof) == isofs.end()) {
		return true;
	}

	DGOutEdgeIter out_i, out_end;
	tie(out_i, out_end) = out_edges(cur_vert, graph);

	while (out_i != out_end) {
		if (_isof_start_ok(graph, target(*out_i, graph), isofs, isof)) {
			return true;
		}
	}

	return false;

}

inline bool isof_start_ok(DirectedGraph const &graph, uint vert,
		IsoformMap const &isofs) {
	Isoform isof(num_vertices(graph));
	return _isof_start_ok(graph, vert, isofs, isof);
}

inline void isoform_MCMC_init(vector<SpliceGraph> &graphs, gsl_rng *rn,
		vector<IsoformMap> &graph_isoforms) {

	vector<IsoformMap>::iterator isof_set_iter = graph_isoforms.begin();

	for (vector<SpliceGraph>::const_iterator i = graphs.begin();
			i != graphs.end(); ++i, ++isof_set_iter) {

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
			isof_set_iter->insert(make_pair(isof, 0.0L));

			for (uint j = 0; j != rc_size; ++j) {
				if (satisfied_rc[j] == true) {
					if (isof == (isof | i->read_constraints[j])) {
						satisfied_rc[j] = false;
					}
				}
			}

		}

	}

	// graph contraints given isoforms

	vector<IsoformMap>::const_iterator graph_isof_iter = graph_isoforms.begin();
	for (vector<SpliceGraph>::iterator i = graphs.begin(); i != graphs.end();
			++i, ++graph_isof_iter) {
		for (vector<uint>::const_iterator j = i->start_nodes.begin();
				j != i->start_nodes.end(); ++j) {
			i->vert_start_ok[*j] = isof_start_ok(i->graph, *j,
					*graph_isof_iter);
		}

		// TODO: remove
		std::cout << "######\n" << std::endl;
		for (IsoformMap::const_iterator j = graph_isof_iter->begin();
				j != graph_isof_iter->end(); ++j) {
			std::cout << j->first << " ";
		}
		std::cout << std::endl;
		for (vector<bool>::const_iterator j = i->vert_start_ok.begin();
				j != i->vert_start_ok.end(); ++j) {
			std::cout << *j << " ";
		}
		std::cout << std::endl;

	}

	// assign random expression levels according
	// to a dirichlet distribution

	uint isofs_size = 0;
	for (isof_set_iter = graph_isoforms.begin();
			isof_set_iter != graph_isoforms.end(); ++isof_set_iter) {
		isofs_size += isof_set_iter->size();
	}

	double *dir_alpha = new double[isofs_size];

	std::fill(dir_alpha, dir_alpha + isofs_size, 1);

	double *dir_theta = new double[isofs_size];

	gsl_ran_dirichlet(rn, isofs_size, dir_alpha, dir_theta);

	delete[] dir_alpha;

	uint isof_exp_ind = 0;

	for (isof_set_iter = graph_isoforms.begin();
			isof_set_iter != graph_isoforms.end(); ++isof_set_iter) {

		for (IsoformMap::iterator cur_isof = isof_set_iter->begin();
				cur_isof != isof_set_iter->end(); ++cur_isof) {

			cur_isof->second = dir_theta[isof_exp_ind];
			++isof_exp_ind;

		}

	}

	delete[] dir_theta;

}

template<class RNodeLoc>
double get_graph_weight(SpliceGraph const &graph, IsoformMap const &graph_isof,
		vector<ReadInGraph<RNodeLoc> > const &read_in_graph,
		GraphReads const &graph_read);

template<class RNodeLoc>
void isoform_main(vector<GraphInfo> const &graph_infos,
		vector<SpliceGraph> &graphs,
		vector<ReadInGraph<RNodeLoc> > const &read_in_graph,
		vector<GraphReads> const &graph_reads, uint max_run) {

	// TODO: multiple chains pthread parallelization
	{

		gsl_rng *rn = gsl_rng_alloc(gsl_rng_mt19937);

		// TODO: seed @rn

		vector<IsoformMap> graph_isoforms(graphs.size());

		isoform_MCMC_init(graphs, rn, graph_isoforms);

		// main part of MCMC
		// the real stuff is in @isof_jump

		vector<vector<IsoformMap> > mcmc_results;

		// number of graphs
		uint graph_num = graphs.size();

		// used to choose a graph to
		// modify the graph's isoform set
		double *graph_weights = new double[graph_num];

		for (uint runs = 0; runs != max_run; ++runs) {

			uint graph_weight_ind = 0;
			vector<SpliceGraph>::const_iterator graph_iter = graphs.begin();
			vector<IsoformMap>::const_iterator graph_isof_iter =
					graph_isoforms.begin();
			vector<GraphReads>::const_iterator graph_read_iter =
					graph_reads.begin();

			while (graph_weight_ind != graph_num) {

				graph_weights[graph_weight_ind] = get_graph_weight(*graph_iter,
						*graph_isof_iter, read_in_graph, *graph_read_iter);

				++graph_weight_ind;
				++graph_iter;
				++graph_isof_iter;
				++graph_read_iter;
			}

			gsl_ran_discrete_t *choose_graph = gsl_ran_discrete_preproc(
					graph_num, graph_weights);

			uint chosen_graph_ind = gsl_ran_discrete(rn, choose_graph);

			gsl_ran_discrete_free(choose_graph);

			double chosen_graph_amount = 0;
			for (IsoformMap::const_iterator i =
					graph_isoforms[chosen_graph_ind].begin();
					i != graph_isoforms[chosen_graph_ind].end(); ++i) {
				chosen_graph_amount += i->second;
			}

		}

		delete[] graph_weights;

		gsl_rng_free(rn);

	}

// TODO: put @mcmc_results into @graphs

}

}

#endif // _MCMC_0_H_
