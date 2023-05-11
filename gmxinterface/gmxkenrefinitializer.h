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
	virtual ~GmxKEnRefInitializer();
	GmxKEnRefInitializer(const GmxKEnRefInitializer &other);
	GmxKEnRefInitializer(GmxKEnRefInitializer &&other);
	GmxKEnRefInitializer& operator=(const GmxKEnRefInitializer &other);
	GmxKEnRefInitializer& operator=(GmxKEnRefInitializer &&other);
	static std::map<std::string, std::vector<int>> loadGmxIndexFile(std::string indexFileName);
	static std::vector<int> loadGmxIndexGroup(std::string groupName, std::string indexFileName="KEnRefAtomIndex.ndx");
};

#endif /* GMXKENREFINITIALIZER_H_ */
