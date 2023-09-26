/*
 * gmxwrapper.h
 *
 *  Created on: Nov 29, 2022
 *      Author: amr
 */

#ifndef GMXWRAPPER_H_
#define GMXWRAPPER_H_

#include "../config/KEnRefConfig.h"

class GMXWrapper { //This class should replace GROMACS main class
public:
	GMXWrapper();
	virtual ~GMXWrapper();
	GMXWrapper(const GMXWrapper &other);
	GMXWrapper(GMXWrapper &&other) noexcept ;
	GMXWrapper& operator=(const GMXWrapper &other);
	GMXWrapper& operator=(GMXWrapper &&other) noexcept;
};

#endif /* GMXWRAPPER_H_ */
