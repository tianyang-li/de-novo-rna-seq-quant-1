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

#include "_graph_seq_0.h"
#include "_misc_0.h"

namespace _mcmc_0 {

using _graph_seq_0::SpliceGraph;
using _graph_seq_0::GraphReads;
using std::vector;
using _graph_seq_0::ReadInGraph;
using _graph_seq_0::check_isoform_read_constraint;

typedef _graph_seq_0::PyGraph GraphInfo;

// calculate the probability that a read
// is from a the transcripts (isoforms)
template<class RNodeLoc>
class ReadFromTransProb {
	inline ldbl operator()(ReadInGraph<RNodeLoc> const &r) {
		return 0.0;
	}
};

template<class RNodeLoc>
void isoform_main(vector<GraphInfo> const &graph_info,
		vector<SpliceGraph> &graphs,
		vector<ReadInGraph<RNodeLoc> > const &read_in_graph,
		vector<GraphReads> const &graph_reads,
		ReadFromTransProb<RNodeLoc> &r_prob, uint max_run) {
	uint runs = 0;

	while (runs != max_run) {

		++runs;
	}
}

}

#endif // _MCMC_0_H_
