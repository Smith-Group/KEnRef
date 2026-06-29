/*
 * gmxkenrefinitializer.cpp
 *
 *  Created on: Apr 4, 2023
 *      Author: amr
 */

#include "gmxinterface/gmxkenrefinitializer.h"
#include "gromacs/topology/index.h"

#include <filesystem>
#include <stdexcept>

GmxKEnRefInitializer::GmxKEnRefInitializer() = default;
GmxKEnRefInitializer::~GmxKEnRefInitializer() = default;
GmxKEnRefInitializer::GmxKEnRefInitializer(const GmxKEnRefInitializer &other) {}
GmxKEnRefInitializer::GmxKEnRefInitializer(GmxKEnRefInitializer &&other)  noexcept {}

//GmxKEnRefInitializer& GmxKEnRefInitializer::operator=(const GmxKEnRefInitializer &other) {}
//GmxKEnRefInitializer& GmxKEnRefInitializer::operator=(GmxKEnRefInitializer &&other) {}

std::map<std::string, std::vector<int>>
GmxKEnRefInitializer::loadGmxIndexFile(const std::string& indexFileName){
    std::map<std::string, std::vector<int>> groupAtoms;
    /* Read index file. GROMACS' init_index() does not null-check the FILE*, so an unreadable path
     * segfaults inside get_a_line()/fgets(). Validate first and fail with a clear, actionable error
     * (the path is resolved to absolute so a wrong working-directory is obvious). */
    {
        std::error_code ec;
        const std::filesystem::path p(indexFileName);
        if (!std::filesystem::is_regular_file(p, ec)) {
            throw std::runtime_error(
                "Index file not found or not readable: '" + indexFileName + "' (resolved to '" +
                std::filesystem::absolute(p, ec).string() + "'). Check the 'index' path / working directory.");
        }
    }
    std::vector<IndexGroup> indexGroups = init_index(indexFileName.c_str());
    for(const auto& group: indexGroups){
    	auto gName = group.name;
    	groupAtoms[gName] = group.particleIndices;
    }
    return groupAtoms;
}

//returns ZERO indexed index file
std::vector<int>
GmxKEnRefInitializer::loadGmxIndexGroup(const std::string& groupName, const std::string& indexFileName){
	auto groups = loadGmxIndexFile(indexFileName);
	return groups[groupName];
}
