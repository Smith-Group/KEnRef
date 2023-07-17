/*
 * gmxkeniefinitializer.h
 *
 *  Created on: Apr 4, 2023
 *      Author: amr
 */

#ifndef GMXKENREFINITIALIZER_H_
#define GMXKENREFINITIALIZER_H_

#include "../core/kenrefinitializer.h"
#include <vector>
#include <map>

class GmxKEnRefInitializer: public DefaultKEnRefInitializer {
public:
	GmxKEnRefInitializer();
	~GmxKEnRefInitializer() override;
	GmxKEnRefInitializer(const GmxKEnRefInitializer &other);
	GmxKEnRefInitializer(GmxKEnRefInitializer &&other) noexcept ;
//	GmxKEnRefInitializer& operator=(const GmxKEnRefInitializer &other);
//	GmxKEnRefInitializer& operator=(GmxKEnRefInitializer &&other);
	static std::map<std::string, std::vector<int>> loadGmxIndexFile(const std::string& indexFileName);
	static std::vector<int> loadGmxIndexGroup(const std::string& groupName, const std::string& indexFileName="KEnRefAtomIndex.ndx");
};

#endif /* GMXKENREFINITIALIZER_H_ */
