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
#include <iostream>
#include <cmath>
#include <exception>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/unordered_map.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <boost/unordered_set.hpp>

#include "_graph_seq_0.h"
#include "_misc_0.h"

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
using std::min;
using std::cerr;
using std::endl;
using _graph_seq_0::ReadGraphLoc;
using _graph_seq_0::get_isof_len;
using _graph_seq_0::IsoformHash;
using std::exception;
using _graph_seq_0::SeqConstraintHash;

}

namespace _mcmc_0 {

#ifdef DEBUG
class ReadConstraintError: public exception {
public:
	virtual char const * what() const throw () {
		return "ReadConstraintError:\nRead constraint violated!\n";
	}
};

class IteratorEndError: public exception {
public:
	virtual char const * what() const throw () {
		return "IteratorEndError:\n iterator didn't reach end as expected!\n";
	}
};
#endif

// when choosing discrete probabilities,
// if something is 0 and gets chosen too many times,
// this exception is thrown
class ZeroProbTooManyTimes: public exception {
public:
	virtual char const * what() const throw () {
		return "ZeroProbTooManyTimes:\n"
				"Something with 0 probability "
				"in a discrete distribution "
				"has been chosen "
				"too many times!\n";
	}

	static ulong const kMaxZeroProbTimes = 5;
};

class GSLRngUnifInt {
public:
	GSLRngUnifInt(gsl_rng * const rn_) :
			rn(rn_) {
	}

	ulong operator()(ulong n) {
		return gsl_rng_uniform_int(rn, n);
	}

private:
	gsl_rng * const rn;
};

// given read constraint
// and graph, get a random isoform
// that satisfies the constraint
inline void rand_rc_isof(SpliceGraph const &graph, Isoform &isof, ulong un_rc,
		gsl_rng * const rn) {
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

		ulong cur_node = *ts_iter;
		isof.set(cur_node);

		DGInEdgeIter in_i, in_end;
		tie(in_i, in_end) = in_edges(cur_node, graph.graph);

		while (in_i != in_end) {
			ulong go2node = gsl_rng_uniform_int(rn,
					in_degree(cur_node, graph.graph));

			for (ulong i = 0; i != go2node; ++i) {
				++in_i;
			}

			cur_node = source(*in_i, graph.graph);
			isof.set(cur_node);
			tie(in_i, in_end) = in_edges(cur_node, graph.graph);
		}
	}

	// last node that has a constraint
	ulong last_rc_node = *ts_iter;
	++ts_iter;

	// extend "back up" to
	// connect to the previous contraint
	while (ts_iter != graph.topo_sort.end()) {
		if (rc[*ts_iter] == true) {
			ulong cur_node = *ts_iter;

			isof.set(cur_node);

			while (cur_node != last_rc_node) {

				vector<ulong> edges_in;

				DGInEdgeIter in_i, in_end;
				tie(in_i, in_end) = in_edges(cur_node, graph.graph);

				while (in_i != in_end) {
					edges_in.push_back(source(*in_i, graph.graph));
					++in_i;
				}

				GSLRngUnifInt rn_wrapper(rn);
				std::random_shuffle(edges_in.begin(), edges_in.end(),
						rn_wrapper);

				for (vector<ulong>::const_iterator i = edges_in.begin();
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

		ulong cur_node = last_rc_node;

		isof.set(cur_node);

		DGOutEdgeIter out_i, out_end;
		tie(out_i, out_end) = out_edges(cur_node, graph.graph);

		while (out_i != out_end) {
			ulong go2node = gsl_rng_uniform_int(rn,
					out_degree(cur_node, graph.graph));

			for (ulong i = 0; i != go2node; ++i) {
				++out_i;
			}

			cur_node = target(*out_i, graph.graph);
			isof.set(cur_node);
			tie(out_i, out_end) = out_edges(cur_node, graph.graph);
		}
	}
}

inline bool _isof_start_ok(DirectedGraph const &graph, ulong cur_vert,
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

inline bool isof_start_ok(DirectedGraph const &graph, ulong vert,
		IsoformMap const &isofs) {
	Isoform isof(num_vertices(graph));
	return _isof_start_ok(graph, vert, isofs, isof);
}

// determine if a read is compatible with an isoform
template<class RNodeLoc>
inline bool read_on_isoform(Isoform const &isof, RNodeLoc const &read_loc);

template<class RNodeLoc>
inline void get_prop_graph_ratio(
		vector<ReadInGraph<RNodeLoc> > const &read_in_graph,
		GraphReads const &graph_read, GraphInfo const &graph_info,
		SpliceGraph const &graph, IsoformMap const &graph_isoform,
		IsoformMap &prop_graph_ratio /* empty map */,
		unordered_map<Isoform, ulong, IsoformHash> &isof_lens) {

	// rescue method
	//
	// Ali Mortazavi, Brian A. Williams, Kenneth McCue, Lorian Schaeffer,
	// and Barbara Wold. Mapping and quantifying mammalian
	// transcriptomes by RNA-seq. Nature Methods, 5(7):621–628, May 2008.

	// TODO: change how weights are chosen?
	// XXX: are weights correct?

	double const kIsofInitWeight = 1.0;

	for (IsoformMap::const_iterator i = graph_isoform.begin();
			i != graph_isoform.end(); ++i) {

		prop_graph_ratio.insert(make_pair(i->first, kIsofInitWeight));

		if (isof_lens.find(i->first) == isof_lens.end()) {
			isof_lens.insert(
					make_pair(i->first, get_isof_len(i->first, graph_info)));
		}

	}

	double tot_uniq_weight = kIsofInitWeight * double(graph_isoform.size());

	// set unique weight
	for (vector<ReadIndex>::const_iterator i = graph_read.reads.begin();
			i != graph_read.reads.end(); ++i) {

		// # of isoforms this read is in
		ulong num_read_isofs = 0;

		IsoformMap::iterator isof_iter = prop_graph_ratio.begin();

		// an isoform that this read is on
		IsoformMap::iterator read_isof;

		while (isof_iter != prop_graph_ratio.end()) {

			if (read_on_isoform(isof_iter->first,
					read_in_graph[i->read_id].graph_locs[i->graph_index].locs[i->align_index])) {

				read_isof = isof_iter;
				++num_read_isofs;

				if (num_read_isofs == 2) {
					break;
				}

			}

			++isof_iter;

		}

		if (num_read_isofs == 1) {

			double const kUniqueWeight = 1.0;

			read_isof->second += kUniqueWeight
					/ double(isof_lens.find(read_isof->first)->second);

			tot_uniq_weight += read_isof->second;

		}

	}

	// set multi-read weights
	for (vector<ReadIndex>::const_iterator i = graph_read.reads.begin();
			i != graph_read.reads.end(); ++i) {

		// isoforms that this read is in
		vector<IsoformMap::iterator> read_isofs;

		for (IsoformMap::iterator j = prop_graph_ratio.begin();
				j != prop_graph_ratio.end(); ++j) {

			if (read_on_isoform(j->first,
					read_in_graph[i->read_id].graph_locs[i->graph_index].locs[i->align_index])) {
				read_isofs.push_back(j);
			}

		}

		if (read_isofs.size() >= 2) {
			for (vector<IsoformMap::iterator>::iterator j = read_isofs.begin();
					j != read_isofs.end(); ++j) {

				double const kMultiReadWeight = 1.0;

				double uniq_weight = (*j)->second;

				(*j)->second += uniq_weight / tot_uniq_weight * kMultiReadWeight
						/ double(isof_lens.find((*j)->first)->second);

			}
		}

	}

	for (IsoformMap::iterator i = prop_graph_ratio.begin();
			i != prop_graph_ratio.end(); ++i) {
		i->second /= double(graph_read.reads.size());
	}

}

template<class RNodeLoc>
inline void isoform_MCMC_init(
		vector<ReadInGraph<RNodeLoc> > const &read_in_graph,
		vector<GraphReads> const &graph_reads,
		vector<GraphInfo> const &graph_infos, vector<SpliceGraph> const &graphs,
		gsl_rng * const rn, vector<IsoformMap> &graph_isoforms,
		vector<IsoformMap> &prop_graph_ratios /* each entry is an empty map */,
		vector<unordered_map<Isoform, ulong, IsoformHash> > &graph_isof_lens,
		vector<unordered_map<SeqConstraint, ulong, SeqConstraintHash> > &rc_isof_counts,
		vector<double *> &vec_vert_start_probs,
		vector<double *> &vec_isof_del_probs) {

#ifdef DEBUG
	cerr << "enter isoform_MCMC_init\n";
#endif

	try {
		// get isoforms that satisfy the read constraints

		vector<IsoformMap>::iterator isof_set_iter = graph_isoforms.begin();

		for (vector<SpliceGraph>::const_iterator i = graphs.begin();
				i != graphs.end(); ++i, ++isof_set_iter) {

			ulong rc_size = i->read_constraints.size();

			// 1 - unsatisfied
			// 0 - satisfied
			dynamic_bitset<> satisfied_rc(rc_size);
			satisfied_rc.set();

			while (satisfied_rc.any()) {
				ulong un_rc = 0; // un-satisfied read constraint index
				while (satisfied_rc[un_rc] == false) {
					++un_rc;
				}
				satisfied_rc[un_rc] = false;

				Isoform isof(boost::num_vertices(i->graph));
				rand_rc_isof(*i, isof, un_rc, rn);
				isof_set_iter->insert(make_pair(isof, 0.0L));

				for (ulong j = 0; j != rc_size; ++j) {
					if (satisfied_rc[j] == true) {
						if (isof == (isof | i->read_constraints[j])) {
							satisfied_rc[j] = false;
						}
					}
				}

			}

		}

#ifdef DEBUG
		if (isof_set_iter != graph_isoforms.end()) {
			cerr << "isof_set_iter" << endl;
			throw IteratorEndError();
		}
#endif

	} catch (exception &e) {
		cerr << e.what() << endl;
		throw;
	}

	try {

		// set @rc_isof_counts
		// count the number of isoforms that satisfy a
		// sequence constraint for each sequence constraint

		vector<SpliceGraph>::const_iterator graph_iter = graphs.begin();
		vector<IsoformMap>::const_iterator graph_isof_iter =
				graph_isoforms.begin();

		for (vector<unordered_map<SeqConstraint, ulong, SeqConstraintHash> >::iterator i =
				rc_isof_counts.begin(); i != rc_isof_counts.end();
				++i, ++graph_iter, ++graph_isof_iter) {

			for (vector<SeqConstraint>::const_iterator j =
					graph_iter->read_constraints.begin();
					j != graph_iter->read_constraints.end(); ++j) {

				ulong rc_isof_count = 0;

				for (IsoformMap::const_iterator k = graph_isof_iter->begin();
						k != graph_isof_iter->end(); ++k) {

					if (k->first == (k->first | (*j))) {
						++rc_isof_count;
					}

				}

				i->insert(make_pair(*j, rc_isof_count));

			}

		}

#ifdef DEBUG
		if (graph_iter != graphs.end()) {
			cerr << "graph_iter" << endl;
			throw IteratorEndError();
		}
		if (graph_isof_iter != graph_isoforms.end()) {
			cerr << "graph_isof_iter" << endl;
			throw IteratorEndError();
		}
#endif

	} catch (exception &e) {
		cerr << e.what() << endl;
	}

	try {

		// assign random expression levels according
		// to a dirichlet distribution

		ulong isofs_size = 0;
		for (vector<IsoformMap>::const_iterator isof_set_iter =
				graph_isoforms.begin(); isof_set_iter != graph_isoforms.end();
				++isof_set_iter) {
			isofs_size += isof_set_iter->size();
		}

		if (isofs_size != 1) {

			double *dir_alpha = new double[isofs_size];

			fill(dir_alpha, dir_alpha + isofs_size, 1);

			double *dir_theta = new double[isofs_size];

			gsl_ran_dirichlet(rn, isofs_size, dir_alpha, dir_theta);

			delete[] dir_alpha;

			ulong isof_exp_ind = 0;
			ulong graph_ind = 0;

			for (vector<IsoformMap>::iterator isof_set_iter =
					graph_isoforms.begin();
					isof_set_iter != graph_isoforms.end();
					++isof_set_iter, ++graph_ind) {

				for (IsoformMap::iterator cur_isof = isof_set_iter->begin();
						cur_isof != isof_set_iter->end(); ++cur_isof) {

					cur_isof->second = dir_theta[isof_exp_ind];

					++isof_exp_ind;

				}

			}

#ifdef DEBUG
			if (graph_ind != graphs.size()) {
				cerr << "graph_ind" << endl;
				throw IteratorEndError();
			}
			if (isof_exp_ind != isofs_size) {
				cerr << "isof_exp_ind" << endl;
				throw IteratorEndError();
			}
#endif

			delete[] dir_theta;

		} else {

			// special case when there's only 1 isoform

			graph_isoforms[0].begin()->second = 1.0;

		}

	} catch (exception &e) {
		cerr << e.what() << endl;
		throw;
	}

	try {
		// get proposal ratios for graphs to
		// use in proposal distributions

		vector<IsoformMap>::const_iterator graph_isof_iter =
				graph_isoforms.begin();
		vector<SpliceGraph>::const_iterator graph_iter = graphs.begin();
		vector<GraphInfo>::const_iterator graph_info_iter = graph_infos.begin();
		vector<GraphReads>::const_iterator graph_read_iter =
				graph_reads.begin();
		vector<unordered_map<Isoform, ulong, IsoformHash> >::iterator isof_lens_iter =
				graph_isof_lens.begin();

		for (vector<IsoformMap>::iterator i = prop_graph_ratios.begin();
				i != prop_graph_ratios.end();
				++i, ++graph_iter, ++graph_isof_iter, ++graph_info_iter, ++graph_read_iter, ++isof_lens_iter) {

			get_prop_graph_ratio(read_in_graph, *graph_read_iter,
					*graph_info_iter, *graph_iter, *graph_isof_iter, *i,
					*isof_lens_iter);

		}

#ifdef DEBUG
		if (graph_isof_iter != graph_isoforms.end()) {
			cerr << "graph_isof_iter" << endl;
			throw IteratorEndError();
		}
		if (graph_iter != graphs.end()) {
			cerr << "graph_iter" << endl;
			throw IteratorEndError();
		}
		if (graph_info_iter != graph_infos.end()) {
			cerr << "graph_info_iter" << endl;
			throw IteratorEndError();
		}
		if (graph_read_iter != graph_reads.end()) {
			cerr << "graph_read_iter" << endl;
			throw IteratorEndError();
		}
		if (isof_lens_iter != graph_isof_lens.end()) {
			cerr << "isof_lens_iter" << endl;
			throw IteratorEndError();
		}
#endif

	} catch (exception &e) {
		cerr << e.what() << endl;
		throw;
	}

	try { // TODO

		{
			vector<SpliceGraph>::const_iterator graph_iter = graphs.begin();

			for (vector<double *>::iterator i = vec_vert_start_probs.begin();
					i != vec_vert_start_probs.end(); ++i, ++graph_iter) {

				(*i) = new double[num_vertices(graph_iter->graph)];

			}

#ifdef DEBUG
			if (graph_iter != graphs.end()) {
				cerr << "graph_iter" << endl;
				throw IteratorEndError();
			}
#endif
		}

		{
			vector<IsoformMap>::const_iterator graph_isof_iter =
					graph_isoforms.begin();

			for (vector<double *>::iterator i = vec_isof_del_probs.begin();
					i != vec_isof_del_probs.end(); ++i, ++graph_isof_iter) {

				(*i) = new double[graph_isof_iter->size()];

			}

#ifdef DEBUG
			if (graph_isof_iter != graph_isoforms.end()) {
				cerr << "graph_isof_iter" << endl;
				throw IteratorEndError();
			}
#endif
		}

	} catch (exception &e) {
		cerr << e.what() << endl;
		throw;
	}

#ifdef DEBUG
	cerr << "exit isoform_MCMC_init\n";
#endif

}

// get graph weights for choosing
// which graph to modify
template<class RNodeLoc>
inline double get_graph_weight(SpliceGraph const &graph,
		IsoformMap const &graph_isof,
		vector<ReadInGraph<RNodeLoc> > const &read_in_graph,
		GraphReads const &graph_read);

// get graph weights for
// updating graph expression values
// when modifying a graph
template<class RNodeLoc>
inline void get_dir_graph_weights(vector<GraphInfo> const &graph_infos,
		vector<SpliceGraph> const &graphs,
		vector<ReadInGraph<RNodeLoc> > const &read_in_graph,
		vector<GraphReads> const &graph_reads,
		double * const dir_graph_weights) {
	// TODO: test this

#ifdef DEBUG
	cerr << "enter get_dir_graph_weights\n";
#endif

	// used to update expression value
	// of each graph
	// using a dirichlet distribution
	//
	// based on the method used in this paper
	//
	// Ali Mortazavi, Brian A. Williams, Kenneth McCue, Lorian Schaeffer,
	// and Barbara Wold. Mapping and quantifying mammalian
	// transcriptomes by RNA-seq. Nature Methods, 5(7):621–628, May 2008.

	// TODO: change how weights are chosen???

	ulong graph_num = graphs.size();

	if (graph_num != 1) {

		// initial weight a graph gets
		double const kGraphWeightInit = 1.0;
		fill(dir_graph_weights, dir_graph_weights + graph_num,
				kGraphWeightInit);

		vector<GraphInfo>::const_iterator graph_info_iter = graph_infos.begin();
		vector<SpliceGraph>::const_iterator graph_iter = graphs.begin();
		vector<GraphReads>::const_iterator graph_read_iter =
				graph_reads.begin();
		ulong cur_graph = 0;

		double tot_unique_weight = double(graph_num) * kGraphWeightInit;

		while (cur_graph != graph_num) {

			// get initial weight for unique reads
			for (vector<ReadIndex>::const_iterator i =
					graph_read_iter->reads.begin();
					i != graph_read_iter->reads.end(); ++i) {

				if (i->get_align_num(read_in_graph) == 1) {
					// weight of a unique read
					double const kUniqueWeight = 1.0;

					dir_graph_weights[cur_graph] += kUniqueWeight
							/ double(graph_info_iter->get_gene_len());

				}

			}

			tot_unique_weight += dir_graph_weights[cur_graph];

			++cur_graph;
			++graph_iter;
			++graph_info_iter;
			++graph_read_iter;
		}

#ifdef DEBUG
		if (cur_graph != graphs.size()) {
			cerr << "cur_graph" << endl;
			throw IteratorEndError();
		}
		if (graph_info_iter != graph_infos.end()) {
			cerr << "graph_info_iter" << endl;
			throw IteratorEndError();
		}
		if (graph_iter != graphs.end()) {
			cerr << "graph_iter" << endl;
			throw IteratorEndError();
		}
		if (graph_read_iter != graph_reads.end()) {
			cerr << "graph_read_iter" << endl;
			throw IteratorEndError();
		}
#endif

		graph_info_iter = graph_infos.begin();
		graph_iter = graphs.begin();
		graph_read_iter = graph_reads.begin();
		cur_graph = 0;

		while (cur_graph != graph_num) {

			// adjust weight for multi-reads
			for (vector<ReadIndex>::const_iterator i =
					graph_read_iter->reads.begin();
					i != graph_read_iter->reads.end(); ++i) {

				ulong align_num = i->get_align_num(read_in_graph);

				if (align_num != 1) {
					// weight of a unique read
					double const kMultiWeight = 1.0;

					double unique_weight = dir_graph_weights[cur_graph];

					dir_graph_weights[cur_graph] += kMultiWeight
							* double(align_num)
							* (unique_weight / tot_unique_weight)
							/ double(graph_info_iter->get_gene_len());
				}

			}

			dir_graph_weights[cur_graph] /= double(read_in_graph.size());

			++cur_graph;
			++graph_iter;
			++graph_info_iter;
			++graph_read_iter;
		}

#ifdef DEBUG
		if (cur_graph != graphs.size()) {
			cerr << "cur_graph" << endl;
			throw IteratorEndError();
		}
		if (graph_info_iter != graph_infos.end()) {
			cerr << "graph_info_iter" << endl;
			throw IteratorEndError();
		}
		if (graph_iter != graphs.end()) {
			cerr << "graph_iter" << endl;
			throw IteratorEndError();
		}
		if (graph_read_iter != graph_reads.end()) {
			cerr << "graph_read_iter" << endl;
			throw IteratorEndError();
		}
#endif

	} else {

		dir_graph_weights[0] = 1;

	}

#ifdef DEBUG
	for (ulong i = 0; i != graph_num; ++i) {
		cerr << dir_graph_weights[i] << " ";
	}
	cerr << endl;
	cerr << "exit get_dir_graph_weights\n";
#endif

}

template<class RNodeLoc>
inline ulong choose_graph_to_mod(gsl_rng * const rn,
		vector<SpliceGraph> const &graphs,
		vector<IsoformMap> const &graph_isoforms,
		vector<GraphReads> const &graph_reads,
		vector<ReadInGraph<RNodeLoc> > const &read_in_graph,
		double * const graph_weights) {

	ulong graph_weight_ind = 0;
	vector<SpliceGraph>::const_iterator graph_iter = graphs.begin();
	vector<IsoformMap>::const_iterator graph_isof_iter = graph_isoforms.begin();
	vector<GraphReads>::const_iterator graph_read_iter = graph_reads.begin();

	ulong graph_num = graphs.size();

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

	ulong chosen_graph_ind = gsl_ran_discrete(rn, choose_graph);

	gsl_ran_discrete_free(choose_graph);

	return chosen_graph_ind;
}

// set the weight for a vertex
// based on its @dist_from_start
inline double dist_from_start_weight(ulong cur_vert, SpliceGraph const &graph) {

	// dist from start weight
	double dfs_w = 1.0;

	// TODO

	return dfs_w;

}

// do DFS setting the weight for an vertex
// about whether isoforms can start from the vertex
// for @get_possible_isoform_weigh
template<class RNodeLoc>
inline void get_possible_isoform_weigh_DFS(Isoform &search_isof,
		double &poss_isof_w, IsoformMap const &prop_graph_ratio,
		GraphInfo const &graph_info, SpliceGraph const &graph,
		GraphReads const &graph_read,
		vector<ReadInGraph<RNodeLoc> > const &read_in_graph, ulong cur_vert,
		bool &can_add_isof, vector<bool> &visited) {

	if (!visited[cur_vert] /* not visited in DFS */
			|| (visited[cur_vert] && !can_add_isof) /* visited, but hasn't added an isoform yet */) {

		search_isof.set(cur_vert);

		if (!can_add_isof) {
			if (prop_graph_ratio.find(search_isof) == prop_graph_ratio.end()) {
				can_add_isof = true;
			}
		}

		DGOutEdgeIter out_i, out_end;
		tie(out_i, out_end) = out_edges(cur_vert, graph.graph);

		while (out_i != out_end) {

			get_possible_isoform_weigh_DFS(search_isof, poss_isof_w,
					prop_graph_ratio, graph_info, graph, graph_read,
					read_in_graph, target(*out_i, graph.graph), can_add_isof,
					visited);

			++out_i;
		}

		if (can_add_isof) {
			// TODO: add weight
		}

		search_isof.reset(cur_vert);

	}
}

// given current isoform set @graph_isoform
// get the weight for vertex @start_vert
// that there will be isoforms added starting from
// @start_vert
//
// the weight are on the same scale as @isof2del_weight
template<class RNodeLoc>
inline double get_possible_isoform_weight(IsoformMap const &prop_graph_ratio,
		GraphInfo const &graph_info, SpliceGraph const &graph,
		GraphReads const &graph_read,
		vector<ReadInGraph<RNodeLoc> > const &read_in_graph, ulong start_vert) {

	double poss_isof_w = 0;

	Isoform search_isof(num_vertices(graph.graph));

	// set flag for search
	// because when @can_id_isof is @false
	// more search than a DFS is needed
	bool can_add_isof = false;

	vector<bool> visited(num_vertices(graph.graph), false);

	get_possible_isoform_weigh_DFS(search_isof, poss_isof_w, prop_graph_ratio,
			graph_info, graph, graph_read, read_in_graph, start_vert,
			can_add_isof, visited);

	return poss_isof_w;
}

// get the probability distribution on the nodes
// for the starting position of an isoform
//
// the weights are on the same scale as @get_isof_del_info
template<class RNodeLoc>
inline void get_vert_start_info(IsoformMap const &prop_graph_ratio,
		GraphInfo const &graph_info, SpliceGraph const &graph,
		GraphReads const &graph_read,
		vector<ReadInGraph<RNodeLoc> > const &read_in_graph,
		double * const vert_start_probs /* filled with 0's */) {

	for (vector<DGVertex>::const_iterator i = graph.topo_sort.begin();
			i != graph.topo_sort.end(); ++i) {

		// weight of possible isoforms starting
		// from vertex @i
		double poss_isof_w = get_possible_isoform_weight(prop_graph_ratio,
				graph_info, graph, graph_read, read_in_graph, *i);

		vert_start_probs[*i] = poss_isof_w * dist_from_start_weight(*i, graph);

	}

}

// check whether or not this @isoform can be removed
// @graph_isoform
inline bool can_remove_isoform(Isoform const &isoform,
		unordered_map<SeqConstraint, ulong, SeqConstraintHash> const &rc_isof_count) {

	for (unordered_map<SeqConstraint, ulong, SeqConstraintHash>::const_iterator i =
			rc_isof_count.begin(); i != rc_isof_count.end(); ++i) {

#ifdef DEBUG
		if (i->second == 0) {
			throw ReadConstraintError();
		}
#endif

		if (isoform == (isoform | i->first)) {
			if (i->second == 1) {
				return false;
			}
		}

	}

	return true;

}

// weight of an isoform if it's to be deleted
// this isoform already makes @inline bool can_remove_isoform
// return @true
//
// weights are on the same scale as @get_possible_isoform_weight
template<class RNodeLoc>
inline double isof2del_weight(Isoform const &isoform,
		IsoformMap const &prop_graph_ratio, GraphInfo const &graph_info,
		SpliceGraph const &graph, GraphReads const &graph_read,
		vector<ReadInGraph<RNodeLoc> > const &read_in_graph);

// get the probability distribution on the
// isoforms for removing an isoform
//
// the weights are on the same scale as @get_vert_start_info
template<class RNodeLoc>
inline void get_isof_del_info(IsoformMap const &prop_graph_ratio,
		GraphInfo const &graph_info, SpliceGraph const &graph,
		GraphReads const &graph_read,
		vector<ReadInGraph<RNodeLoc> > const &read_in_graph,
		double * const isof_del_probs /* filled with 0's */,
		unordered_map<SeqConstraint, ulong, SeqConstraintHash> const &rc_isof_count,
		vector<IsoformMap::const_iterator> &isof_del_probs_ind /* empty vector */) {

	ulong isof_del_prob_ind = 0;

	for (IsoformMap::const_iterator i = prop_graph_ratio.begin();
			i != prop_graph_ratio.end(); ++i, ++isof_del_prob_ind) {

		if (can_remove_isoform(i->first, rc_isof_count)) {

			isof_del_probs[isof_del_prob_ind] = isof2del_weight(i->first,
					prop_graph_ratio, graph_info, graph, graph_read,
					read_in_graph);

		}

		isof_del_probs_ind.push_back(i);

	}

}

// grow an isoform from a starting node
// return the probability of growing that isoform like
// that
template<class RNodeLoc>
inline double grow_added_isof_prob(IsoformMap const &graph_isoform,
		GraphInfo const &graph_info, SpliceGraph const &graph,
		GraphReads const &graph_read,
		vector<ReadInGraph<RNodeLoc> > const &read_in_graph, gsl_rng * const rn,
		Isoform &added_isof) {

	double grow_prob = 1;

	// TODO:

	return grow_prob;
}

// add an isoform and return the corresponding
// MCMC ratio
template<class RNodeLoc>
inline double add_isof_ratio(IsoformMap const &graph_isoform,
		IsoformMap const &prop_graph_ratio, GraphInfo const &graph_info,
		SpliceGraph const &graph, GraphReads const &graph_read,
		vector<ReadInGraph<RNodeLoc> > const &read_in_graph, gsl_rng * const rn,
		IsoformMap &new_graph_isof /* empty, if accepted @graph_isoform <- */,
		IsoformMap &new_prop_ratio /* empty, if accepted @prop_graph_ratio <- */,
		unordered_map<Isoform, ulong, IsoformHash> &isof_lens,
		double * const vert_start_probs /* already set */,
		unordered_map<SeqConstraint, ulong, SeqConstraintHash> &rc_isof_count /* this will not be updated here */) {

	ulong num_graph_vert = num_vertices(graph.graph);

	double model_graph_ratio = 1;

	gsl_ran_discrete_t *vert_start_probs_gsl = gsl_ran_discrete_preproc(
			num_graph_vert, vert_start_probs);

	ulong vert_start = gsl_ran_discrete(rn, vert_start_probs_gsl);
	{
		ulong chose0prob_times = 0;
		while (vert_start_probs[vert_start] == 0) {
			if (chose0prob_times == ZeroProbTooManyTimes::kMaxZeroProbTimes) {
				throw ZeroProbTooManyTimes();
			}
			vert_start = gsl_ran_discrete(rn, vert_start_probs_gsl);
			++chose0prob_times;
		}
	}

	gsl_ran_discrete_free(vert_start_probs_gsl);

	Isoform added_isof(num_graph_vert);

	double added_isof_prob = grow_added_isof_prob(graph_isoform, graph_info,
			graph, graph_read, read_in_graph, rn, added_isof);

	new_graph_isof = graph_isoform;
	new_graph_isof.insert(make_pair(added_isof, 0));

	get_prop_graph_ratio(read_in_graph, graph_read, graph_info, graph,
			new_graph_isof, new_prop_ratio, isof_lens);

	ulong num_new_graph_isof = new_graph_isof.size();
	double *new_isof_del_probs = new double[num_new_graph_isof];
	fill(new_isof_del_probs, new_isof_del_probs + num_new_graph_isof, 0);

	for (unordered_map<SeqConstraint, ulong, SeqConstraintHash>::iterator i =
			rc_isof_count.begin(); i != rc_isof_count.end(); ++i) {
		if (added_isof == (added_isof | i->first)) {
			++i->second;
		}
	}

	vector<IsoformMap::const_iterator> new_isof_del_probs_ind;

	get_isof_del_info(new_prop_ratio, graph_info, graph, graph_read,
			read_in_graph, new_isof_del_probs, rc_isof_count,
			new_isof_del_probs_ind);

	for (unordered_map<SeqConstraint, ulong, SeqConstraintHash>::iterator i =
			rc_isof_count.begin(); i != rc_isof_count.end(); ++i) {
		if (added_isof == (added_isof | i->first)) {

#ifdef DEBUG
			if (i->second == 1) {
				throw ReadConstraintError();
			}
#endif

			--i->second;

		}
	}

	// update @model_graph_ratio
	// TODO

	delete[] vert_start_probs;

	delete[] new_isof_del_probs;

	return model_graph_ratio;
}

// remove an isoform and return the corresponding
// MCMC ratio
template<class RNodeLoc>
inline double del_isof_ratio(IsoformMap const &graph_isoform,
		IsoformMap const &prop_graph_ratio, GraphInfo const &graph_info,
		SpliceGraph const &graph, GraphReads const &graph_read,
		vector<ReadInGraph<RNodeLoc> > const &read_in_graph, gsl_rng * const rn,
		IsoformMap &new_graph_isof /* empty, if accepted @graph_isoform <- */,
		IsoformMap &new_prop_ratio /* empty, if accepted @prop_graph_ratio <- */,
		unordered_map<Isoform, ulong, IsoformHash> &isof_lens,
		double * const isof_del_probs /* this is already set */,
		vector<IsoformMap::const_iterator> const &isof_del_probs_ind) {

	ulong num_graph_vert = num_vertices(graph.graph);

	double model_graph_ratio = 1;

	ulong num_graph_isoform = graph_isoform.size();

	gsl_ran_discrete_t *isof_del_probs_gsl = gsl_ran_discrete_preproc(
			num_graph_isoform, isof_del_probs);

	ulong isof_del = gsl_ran_discrete(rn, isof_del_probs_gsl);
	{
		ulong chose0prob_times = 0;
		while (isof_del_probs[isof_del] == 0) {
			if (chose0prob_times == ZeroProbTooManyTimes::kMaxZeroProbTimes) {
				throw ZeroProbTooManyTimes();
			}
			isof_del = gsl_ran_discrete(rn, isof_del_probs_gsl);
			++chose0prob_times;
		}
	}

	gsl_ran_discrete_free(isof_del_probs_gsl);

	new_graph_isof = graph_isoform;
	new_graph_isof.erase(isof_del_probs_ind[isof_del]->first);

	get_prop_graph_ratio(read_in_graph, graph_read, graph_info, graph,
			new_graph_isof, new_prop_ratio, isof_lens);

	double *new_vert_start_probs = new double[num_graph_vert];
	fill(new_vert_start_probs, new_vert_start_probs + num_graph_vert, 0);

	get_vert_start_info(new_prop_ratio, graph_info, graph, graph_read,
			read_in_graph, new_vert_start_probs);

	// update @model_graph_ratio
	// TODO

	delete[] new_vert_start_probs;

	delete[] isof_del_probs;

	return model_graph_ratio;
}

// (partial) MCMC ratio for adding or deleting an isoform from
// this graph
template<class RNodeLoc>
inline double update_chosen_graph_isoform(IsoformMap const &graph_isoform,
		IsoformMap const &prop_graph_ratio, GraphInfo const &graph_info,
		SpliceGraph const &graph, GraphReads const &graph_read,
		vector<ReadInGraph<RNodeLoc> > const &read_in_graph, gsl_rng * const rn,
		IsoformMap &new_graph_isof /* empty, if accepted @graph_isoform <- */,
		IsoformMap &new_prop_ratio /* empty, if accepted @prop_graph_ratio <- */,
		unordered_map<Isoform, ulong, IsoformHash> &isof_lens,
		unordered_map<SeqConstraint, ulong, SeqConstraintHash> &rc_isof_count) {

	double partial_mcmc_ratio;

	try {

		// return the ratio of adding or removing the isoform, or 1

		double graph_expr_val = 0;
		for (IsoformMap::const_iterator i = graph_isoform.begin();
				i != graph_isoform.end(); ++i) {
			graph_expr_val += i->second;
		}

		double *vert_start_probs = new double[num_vertices(graph.graph)];
		fill(vert_start_probs, vert_start_probs + num_vertices(graph.graph), 0);

		get_vert_start_info(prop_graph_ratio, graph_info, graph, graph_read,
				read_in_graph, vert_start_probs);

		// @_add_isof_weight / (@_add_isof_weight + @_del_isof_weight)
		// is the probability of adding an isoform to the current set
		double add_isof_weight = 0;
		for (ulong i = 0; i != num_vertices(graph.graph); ++i) {
			add_isof_weight += vert_start_probs[i];
		}

		// @_del_isof_weight / (@_add_isof_weight + @_del_isof_weight)
		// is the probability of deleting an isoform from the current set
		double *isof_del_probs = new double[graph_isoform.size()];
		fill(isof_del_probs, isof_del_probs + graph_isoform.size(), 0);

		vector<IsoformMap::const_iterator> isof_del_probs_ind;

		get_isof_del_info(prop_graph_ratio, graph_info, graph, graph_read,
				read_in_graph, isof_del_probs, rc_isof_count,
				isof_del_probs_ind);

		double del_isof_weight = 0;
		for (ulong i = 0; i != graph_isoform.size(); ++i) {
			del_isof_weight += isof_del_probs[i];
		}

		if (add_isof_weight == 0 && del_isof_weight == 0) {
			partial_mcmc_ratio = 1;
		} else {

			enum Action {
				ADD, DEL,
			};

			Action action;

			double action_prob = 1.0;

			if (add_isof_weight == 0) {

				action = ADD;

			} else {

				if (del_isof_weight == 0) {
					action = DEL;
				} else {

					double add_prob = add_isof_weight
							/ (add_isof_weight + del_isof_weight);

					if (gsl_rng_uniform(rn) <= add_prob) {

						action = ADD;

						action_prob = add_prob;

					} else {

						action = DEL;

						action_prob = del_isof_weight
								/ (add_isof_weight + del_isof_weight);

					}

				}

			}

			switch (action) {

			case ADD:
				try {
					partial_mcmc_ratio = add_isof_ratio(graph_isoform,
							prop_graph_ratio, graph_info, graph, graph_read,
							read_in_graph, rn, new_graph_isof, new_prop_ratio,
							isof_lens, vert_start_probs, rc_isof_count);
				} catch (exception &e) {
					cerr << e.what() << endl;
					throw;
				}
				break;

			case DEL:
				try {
					partial_mcmc_ratio = del_isof_ratio(graph_isoform,
							prop_graph_ratio, graph_info, graph, graph_read,
							read_in_graph, rn, new_graph_isof, new_prop_ratio,
							isof_lens, isof_del_probs, isof_del_probs_ind);
				} catch (exception &e) {
					cerr << e.what() << endl;
				}
				break;

			}

			partial_mcmc_ratio *= action_prob;

		}

		delete[] vert_start_probs;
		delete[] isof_del_probs;

	} catch (exception &e) {
		cerr << e.what() << endl;
		throw;
	}

	return partial_mcmc_ratio;

}

// update weights for choosing which graph to modify
// this function should be simple
template<class RNodeLoc>
inline void update_graph_weights(vector<GraphInfo> const &graph_infos,
		vector<SpliceGraph> const &graphs,
		vector<ReadInGraph<RNodeLoc> > const &read_in_graph,
		vector<GraphReads> const &graph_reads, double * const graph_weights,
		vector<IsoformMap> const &graph_isoforms);

template<class RNodeLoc>
inline void isoform_main(vector<GraphInfo> const &graph_infos,
		vector<SpliceGraph> &graphs,
		vector<ReadInGraph<RNodeLoc> > const &read_in_graph,
		vector<GraphReads> const &graph_reads, ulong max_run, ulong num_thread =
				1 /* TODO: add threading */) {

	try {

		// number of graphs
		ulong graph_num = graphs.size();

		// used to update expression value
		// of each graph
		// using a dirichlet distribution
		double *dir_graph_weights = new double[graph_num];

		get_dir_graph_weights(graph_infos, graphs, read_in_graph, graph_reads,
				dir_graph_weights);

		double tot_dir_graph_weight = 0;
		for (ulong i = 0; i != graph_num; ++i) {
			tot_dir_graph_weight += dir_graph_weights[i];
		}

		try {

			gsl_rng *rn = gsl_rng_alloc(gsl_rng_mt19937);

			// TODO: seed @rn

			// for each graph, the expression values of
			// each isoform in that graph
			vector<IsoformMap> graph_isoforms(graph_num);

			// for each graph, the proposal dirichlet parameters
			// for the isoforms
			vector<IsoformMap> prop_graph_ratios(graph_num);

			vector<unordered_map<Isoform, ulong, IsoformHash> > graph_isof_lens(
					graph_num);

			// count the number of isoforms that satisfy a sequence
			// constraint for each sequence constraint
			vector<unordered_map<SeqConstraint, ulong, SeqConstraintHash> > rc_isof_counts(
					graph_num);

			vector<double *> vec_vert_start_probs(graph_num);
			vector<double *> vec_isof_del_probs(graph_num);

			isoform_MCMC_init(read_in_graph, graph_reads, graph_infos, graphs,
					rn, graph_isoforms, prop_graph_ratios, graph_isof_lens,
					rc_isof_counts, vec_vert_start_probs, vec_isof_del_probs);

			vector<vector<IsoformMap> > mcmc_results;

			// used to choose a graph to
			// modify the graph's isoform set
			double *graph_weights = new double[graph_num];

			for (ulong runs = 0; runs != max_run; ++runs) {

				ulong chosen_graph_ind = choose_graph_to_mod(rn, graphs,
						graph_isoforms, graph_reads, read_in_graph,
						graph_weights);

				// if accepted @graph_isoform ->
				IsoformMap new_graph_isof;

				// if accepted @prop_graph_ratio ->
				IsoformMap new_prop_ratio;

				double *new_vert_start_probs;
				double *new_isof_del_probs;

				// the ratio of adding or deleting an isoform
				double model_graph_ratio = update_chosen_graph_isoform(
						graph_isoforms[chosen_graph_ind],
						prop_graph_ratios[chosen_graph_ind],
						graph_infos[chosen_graph_ind], graphs[chosen_graph_ind],
						graph_reads[chosen_graph_ind], read_in_graph, rn,
						new_graph_isof, new_prop_ratio,
						graph_isof_lens[chosen_graph_ind],
						rc_isof_counts[chosen_graph_ind]);

				double new_chosen_graph_portion = 1.0;
				if (graph_num != 1) {

					new_chosen_graph_portion = gsl_ran_beta(rn,
							dir_graph_weights[chosen_graph_ind],
							tot_dir_graph_weight
									- dir_graph_weights[chosen_graph_ind]);

				}

				// TODO: calculate @model_graph_ratio

				if (gsl_rng_uniform(rn) <= min(1.0, model_graph_ratio)) {
					// update to accepted state

					// TODO

					update_graph_weights(graph_infos, graphs, read_in_graph,
							graph_reads, graph_weights, graph_isoforms);

					if (graph_num != 1) {

					}

				}

			}

			delete[] graph_weights;

			gsl_rng_free(rn);

			for (vector<double *>::iterator i = vec_vert_start_probs.begin();
					i != vec_vert_start_probs.begin(); ++i) {
				delete[] *i;
			}

			for (vector<double *>::iterator i = vec_isof_del_probs.begin();
					i != vec_isof_del_probs.end(); ++i) {
				delete[] *i;
			}

		} catch (exception &e) {
			cerr << e.what() << endl;
			throw;
		}

		delete[] dir_graph_weights;

		// TODO: put @mcmc_results into @graphs

	} catch (exception &e) {
		cerr << e.what() << endl;
		throw;
	}

}

} // namespace _mcmc_0

#endif // _MCMC_0_H_
