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

#ifndef _SINGLE_1_H_
#define _SINGLE_1_H_

#include <vector>
#include <boost/dynamic_bitset.hpp>
#include <gsl/gsl_rng.h>

#include "_misc_0.h"
#include "_graph_seq_0.h"
#include "_mcmc_0.h"

namespace _single_1 {

using _graph_seq_0::ReadInGraph;
using _graph_seq_0::ReadGraphLoc;
using _graph_seq_0::GraphReads;
using _graph_seq_0::GraphInfo;

typedef _graph_seq_0::ReadNodeLoc SingleNodeLoc;

/*
 * used to interface with python
 */
void _get_isoforms(std::vector<_graph_seq_0::PyGraph> *py_graphs,
		std::vector<ReadInGraph<SingleNodeLoc> > *py_reads,
		std::vector<_graph_seq_0::Fasta> *isoforms, uint max_run);

}

namespace _mcmc_0 {

using _single_1::SingleNodeLoc;
using _graph_seq_0::ReadInGraph;

template<>
inline double get_graph_weight<SingleNodeLoc>(SpliceGraph const &graph,
		IsoformMap const &graph_isof,
		vector<ReadInGraph<SingleNodeLoc> > const &read_in_graph,
		GraphReads const &graph_read) {

	return 0;

}

}

#endif // _SINGLE_1_H_
