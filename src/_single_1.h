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

#include "_misc_0.h"

namespace _single_1 {

/*
 * makes it easy to get graph information from
 * python
 */
class PyNode {
public:
	PyNode(uint node_id_) :
			node_id(node_id_) {
	}

	PyNode &operator=(PyNode const &x) {
		if (this != &x) {
			node_id = x.node_id;
			edges = x.edges;
		}
		return *this;
	}

	PyNode(PyNode const &x) :
			node_id(x.node_id), edges(x.edges) {
	}

	~PyNode() {
	}

	uint node_id;
	std::vector<uint> edges;
};

class PyGraph {
public:
	PyGraph() {
	}

	uint graph_id;
	std::vector<PyNode> nodes;
};

}

#endif // _SINGLE_1_H_
