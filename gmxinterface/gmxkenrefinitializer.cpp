/*
 * gmxkenrefinitializer.cpp
 *
 *  Created on: Apr 4, 2023
 *      Author: amr
 */

#include "gmxkenrefinitializer.h"
#include "gromacs/topology/index.h"

GmxKEnRefInitializer::GmxKEnRefInitializer() {
	// TODO Auto-generated constructor stub

}

GmxKEnRefInitializer::~GmxKEnRefInitializer() {
	// TODO Auto-generated destructor stub
}

GmxKEnRefInitializer::GmxKEnRefInitializer(const GmxKEnRefInitializer &other) {
	// TODO Auto-generated constructor stub

}

GmxKEnRefInitializer::GmxKEnRefInitializer(GmxKEnRefInitializer &&other) {
	// TODO Auto-generated constructor stub

}

//GmxKEnRefInitializer& GmxKEnRefInitializer::operator=(
//		const GmxKEnRefInitializer &other) {
//	// TODO Auto-generated method stub
//}
//GmxKEnRefInitializer& GmxKEnRefInitializer::operator=(
//		GmxKEnRefInitializer &&other) {
//	// TODO Auto-generated method stub
//}

std::map<std::string, std::vector<int>>
GmxKEnRefInitializer::loadGmxIndexFile(std::string indexFileName){
    std::map<std::string, std::vector<int>> groupAtoms;
    /* Read index file */
    std::vector<IndexGroup> indexGroups = init_index(indexFileName.c_str());
    for(auto group: indexGroups){
    	auto gName = group.name;
    	groupAtoms[gName] = group.particleIndices;
    }
    return groupAtoms;
}

std::vector<int>
GmxKEnRefInitializer::loadGmxIndexGroup(std::string groupName, std::string indexFileName){
	auto groups = loadGmxIndexFile(indexFileName);
	return groups[groupName];
}
