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

#ifndef GRAPH_SEQ_0_H_
#define GRAPH_SEQ_0_H_

namespace graph_seq_0 {

/*
 * assumes that we already know which graph it is
 */
class SeqLoc {
public:
	SeqLoc(int node_id_, unsigned int start_, unsigned int end_) :
			node_id(node_id_), start(start_), end(end_) {
	}

	int node_id;
	unsigned int start;
	unsigned int end;
};

}

#endif // GRAPH_SEQ_0_H_
