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
 */

#ifndef _MISC_0_H_
#define _MISC_0_H_

#include <string>
#include <cstddef>
#include <boost/functional/hash.hpp>
#include <boost/dynamic_bitset.hpp>

typedef unsigned int uint;
typedef long double ldbl;

template<typename Block, typename Allocator>
size_t hash_value(boost::dynamic_bitset<Block, Allocator> const &x) {
	boost::hash<std::string> hasher;
	std::string x_str = x.to_string();
	return hasher(x_str);
}

#endif // _MISC_0_H_
