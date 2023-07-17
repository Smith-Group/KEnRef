/*
 * gmxkenrefinitializer.cpp
 *
 *  Created on: Apr 4, 2023
 *      Author: amr
 */

#include "gmxkenrefinitializer.h"
#include "gromacs/topology/index.h"

GmxKEnRefInitializer::GmxKEnRefInitializer() = default;
GmxKEnRefInitializer::~GmxKEnRefInitializer() = default;
GmxKEnRefInitializer::GmxKEnRefInitializer(const GmxKEnRefInitializer &other) {}
GmxKEnRefInitializer::GmxKEnRefInitializer(GmxKEnRefInitializer &&other)  noexcept {}

//GmxKEnRefInitializer& GmxKEnRefInitializer::operator=(const GmxKEnRefInitializer &other) {}
//GmxKEnRefInitializer& GmxKEnRefInitializer::operator=(GmxKEnRefInitializer &&other) {}

std::map<std::string, std::vector<int>>
GmxKEnRefInitializer::loadGmxIndexFile(const std::string& indexFileName){
    std::map<std::string, std::vector<int>> groupAtoms;
    /* Read index file */
    std::vector<IndexGroup> indexGroups = init_index(indexFileName.c_str());
    for(const auto& group: indexGroups){
    	auto gName = group.name;
    	groupAtoms[gName] = group.particleIndices;
    }
    return groupAtoms;
}

std::vector<int>
GmxKEnRefInitializer::loadGmxIndexGroup(const std::string& groupName, const std::string& indexFileName){
	auto groups = loadGmxIndexFile(indexFileName);
	return groups[groupName];
}
