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
#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_set.hpp>

#include "_single_1.h"
#include "_misc_0.h"
#include "_graph_seq_0.h"

namespace _single_1 {

inline void get_read_constraints(_graph_seq_0::PyGraph const &py_graph,
		GraphReads const &graph_read,
		std::vector<PyReadInGraph> const &py_reads,
		std::vector<boost::dynamic_bitset<> > &rcs /* read constraints */) {

	// non-redundant read constraints
	boost::unordered_set<boost::dynamic_bitset<> > nr_rcs;

	for (std::vector<ReadIndex>::const_iterator i = graph_read.reads.begin();
			i != graph_read.reads.end(); ++i) {
		boost::dynamic_bitset<> rc(py_graph.nodes.size());
		for (_graph_seq_0::PyReadNodeLoc::const_iterator j =
				py_reads[i->read_id].graph_locs[i->graph_index].locs[i->align_index].begin();
				j
						!= py_reads[i->read_id].graph_locs[i->graph_index].locs[i->align_index].end();
				++j) {
		}
	}

}

void _get_isoforms(std::vector<_graph_seq_0::PyGraph> *py_graphs,
		std::vector<PyReadInGraph> *py_reads,
		std::vector<_graph_seq_0::Isoform> *isoforms) {

	// the reads that a graph has
	std::vector<GraphReads> graph_reads(py_graphs->size());

	{
		uint graph_id = 0;

		for (std::vector<GraphReads>::iterator i = graph_reads.begin();
				i != graph_reads.end(); ++i, ++graph_id) {
			i->graph_id = graph_id;
		}
	}

	{
		uint read_id = 0;

		for (std::vector<PyReadInGraph>::const_iterator i = py_reads->begin();
				i != py_reads->end(); ++i, ++read_id) {
			uint graph_index = 0;

			for (std::vector<PyReadGraphLoc>::const_iterator j =
					i->graph_locs.begin(); j != i->graph_locs.end();
					++j, ++graph_index) {
				uint align_index = 0;

				for (std::vector<_graph_seq_0::PyReadNodeLoc>::const_iterator k =
						j->locs.begin(); k != j->locs.end();
						++k, ++align_index) {
					graph_reads[j->graph_id].reads.push_back(
							ReadIndex(read_id, graph_index, align_index));
				}
			}
		}
	}

	std::vector<_graph_seq_0::SpliceGraph> graphs(py_graphs->size());
	{
		uint graph_id = 0;
		std::vector<GraphReads>::const_iterator graph_read =
				graph_reads.begin();
		std::vector<_graph_seq_0::PyGraph>::const_iterator py_graph =
				py_graphs->begin();

		for (std::vector<_graph_seq_0::SpliceGraph>::iterator i =
				graphs.begin(); i != graphs.end();
				++i, ++graph_id, ++graph_read, ++py_graph) {
			i->graph_id = graph_id;
			get_read_constraints(*py_graph, *graph_read, *py_reads,
					i->read_constraints);
		}
	}

}

}

