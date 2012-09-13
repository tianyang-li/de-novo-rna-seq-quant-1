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
using _graph_seq_0::Node;
using std::fill;
using _graph_seq_0::ReadIndex;

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
		vector<IsoformMap> &graph_isoforms, double * const graph_expr_vals) {
	// TODO: remove
	std::cout << "enter isoform_MCMC_init\n";

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

	double* dir_alpha = new double[isofs_size];

	fill(dir_alpha, dir_alpha + isofs_size, 1);

	double* dir_theta = new double[isofs_size];

	gsl_ran_dirichlet(rn, isofs_size, dir_alpha, dir_theta);

	delete[] dir_alpha;

	uint isof_exp_ind = 0;
	uint graph_ind = 0;

	for (isof_set_iter = graph_isoforms.begin();
			isof_set_iter != graph_isoforms.end();
			++isof_set_iter, ++graph_ind) {

		double graph_expr_val = 0;

		for (IsoformMap::iterator cur_isof = isof_set_iter->begin();
				cur_isof != isof_set_iter->end(); ++cur_isof) {

			cur_isof->second = dir_theta[isof_exp_ind];

			graph_expr_val += dir_theta[isof_exp_ind];

			++isof_exp_ind;

		}

		graph_expr_vals[graph_ind] = graph_expr_val;

	}

	delete[] dir_theta;

	// TODO: remove
	std::cout << "exit isoform_MCMC_init\n";

}

template<class RNodeLoc>
inline double get_graph_weight(SpliceGraph const &graph,
		IsoformMap const &graph_isof,
		vector<ReadInGraph<RNodeLoc> > const &read_in_graph,
		GraphReads const &graph_read);

template<class RNodeLoc>
inline void get_dir_graph_weights(vector<GraphInfo> const &graph_infos,
		vector<SpliceGraph> const &graphs,
		vector<ReadInGraph<RNodeLoc> > const &read_in_graph,
		vector<GraphReads> const &graph_reads,
		double* const dir_graph_weights) {
	// TODO: test this

	// TODO: remove
	std::cout << "enter get_dir_graph_weights\n";

	// used to update expression value
	// of each graph
	// using a dirichlet distribution
	//
	// based on the method used in this paper
	//
	// Ali Mortazavi, Brian A. Williams, Kenneth McCue, Lorian Schaeffer,
	// and Barbara Wold. Mapping and quantifying mammalian
	// transcriptomes by RNA-seq. Nature Methods, 5(7):621–628, May 2008.

	uint graph_num = graphs.size();

	// initial weight a graph gets
	double const kGraphWeightInit = 1.0;
	fill(dir_graph_weights, dir_graph_weights + graph_num, kGraphWeightInit);

	vector<GraphInfo>::const_iterator graph_info_iter = graph_infos.begin();
	vector<SpliceGraph>::const_iterator graph_iter = graphs.begin();
	vector<GraphReads>::const_iterator graph_read_iter = graph_reads.begin();
	uint cur_graph = 0;

	double tot_unique_weight = double(graph_num) * kGraphWeightInit;

	while (cur_graph != graph_num) {

		uint graph_len = 0;
		for (vector<Node>::const_iterator i = graph_info_iter->nodes.begin();
				i != graph_info_iter->nodes.end(); ++i) {
			graph_len += i->est_len;
		}

		// get initial weight for unique reads
		for (vector<ReadIndex>::const_iterator i =
				graph_read_iter->reads.begin();
				i != graph_read_iter->reads.end(); ++i) {

			if (i->get_align_num(read_in_graph) == 1) {
				// weight of a unique read
				double const kUniqueWeight = 1.0;
				dir_graph_weights[cur_graph] += kUniqueWeight;
				tot_unique_weight += kUniqueWeight;
			}

		}

		// adjust weight for multi-reads
		for (vector<ReadIndex>::const_iterator i =
				graph_read_iter->reads.begin();
				i != graph_read_iter->reads.end(); ++i) {

			uint align_num = i->get_align_num(read_in_graph);

			if (align_num != 1) {
				// weight of a unique read
				double const kMultiWeight = 1.0;
				double unique_weight = dir_graph_weights[cur_graph];
				dir_graph_weights[cur_graph] +=
						(kMultiWeight * double(align_num) * unique_weight
								/ tot_unique_weight);
			}

		}

		++cur_graph;
		++graph_iter;
		++graph_info_iter;
		++graph_read_iter;
	}

	for (uint i = 0; i != graph_num; ++i) {
		dir_graph_weights[i] /= 8.0;
	}

	// TODO: remove
	for (uint i = 0; i != graph_num; ++i) {
		std::cout << dir_graph_weights[i] << " ";
	}
	std::cout << std::endl;
	std::cout << "exit get_dir_graph_weights\n";
}

template<class RNodeLoc>
inline uint choose_graph_to_mod(gsl_rng *rn, vector<SpliceGraph> const &graphs,
		vector<IsoformMap> const &graph_isoforms,
		vector<GraphReads> const &graph_reads,
		vector<ReadInGraph<RNodeLoc> > const &read_in_graph,
		double* const graph_weights) {

	uint graph_weight_ind = 0;
	vector<SpliceGraph>::const_iterator graph_iter = graphs.begin();
	vector<IsoformMap>::const_iterator graph_isof_iter = graph_isoforms.begin();
	vector<GraphReads>::const_iterator graph_read_iter = graph_reads.begin();

	uint graph_num = graphs.size();

	while (graph_weight_ind != graph_num) {

		graph_weights[graph_weight_ind] = get_graph_weight(*graph_iter,
				*graph_isof_iter, read_in_graph, *graph_read_iter);

		++graph_weight_ind;
		++graph_iter;
		++graph_isof_iter;
		++graph_read_iter;
	}

	gsl_ran_discrete_t *choose_graph = gsl_ran_discrete_preproc(graph_num,
			graph_weights);

	uint chosen_graph_ind = gsl_ran_discrete(rn, choose_graph);

	gsl_ran_discrete_free(choose_graph);

	return chosen_graph_ind;
}

inline void update_isof_expr_val(vector<IsoformMap> &graph_isoforms,
		double* const graph_expr_vals, double* const new_graph_expr_vals) {

	uint graph_ind = 0;

	for (vector<IsoformMap>::iterator i = graph_isoforms.begin();
			i != graph_isoforms.end(); ++i, ++graph_ind) {
		for (IsoformMap::iterator j = i->begin(); j != i->end(); ++j) {

			// XXX: numerical problems
			j->second *= (new_graph_expr_vals[graph_ind]
					/ graph_expr_vals[graph_ind]);

		}
	}
}

template<class RNodeLoc>
inline void isoform_main(vector<GraphInfo> const &graph_infos,
		vector<SpliceGraph> &graphs,
		vector<ReadInGraph<RNodeLoc> > const &read_in_graph,
		vector<GraphReads> const &graph_reads, uint max_run) {

	// number of graphs
	uint graph_num = graphs.size();

	// used to update expression value
	// of each graph
	// using a dirichlet distribution
	double *dir_graph_weights = new double[graph_num];

	get_dir_graph_weights(graph_infos, graphs, read_in_graph, graph_reads,
			dir_graph_weights);

	// TODO: multiple chains pthread parallelization
	{

		gsl_rng *rn = gsl_rng_alloc(gsl_rng_mt19937);

		// TODO: seed @rn

		vector<IsoformMap> graph_isoforms(graphs.size());

		double* graph_expr_vals = new double[graph_num];

		isoform_MCMC_init(graphs, rn, graph_isoforms, graph_expr_vals);

		// main part of MCMC
		// the real stuff is in @isof_jump

		vector<vector<IsoformMap> > mcmc_results;

		// used to choose a graph to
		// modify the graph's isoform set
		double* graph_weights = new double[graph_num];

		double* new_graph_expr_vals = new double[graph_num];

		for (uint runs = 0; runs != max_run; ++runs) {

			uint chosen_graph_ind = choose_graph_to_mod(rn, graphs,
					graph_isoforms, graph_reads, read_in_graph, graph_weights);

			gsl_ran_dirichlet(rn, graph_num, dir_graph_weights,
					new_graph_expr_vals);

			update_isof_expr_val(graph_isoforms, graph_expr_vals,
					new_graph_expr_vals);

			for (uint i = 0; i != graph_num; ++i) {
				graph_expr_vals[i] = new_graph_expr_vals[i];
			}

		}

		delete[] graph_weights;
		delete[] graph_expr_vals;
		delete[] new_graph_expr_vals;

		gsl_rng_free(rn);

	}

	delete[] dir_graph_weights;

	// TODO: put @mcmc_results into @graphs

}

}

#endif // _MCMC_0_H_
