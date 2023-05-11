/*
 * KEnRefRunner.h
 *
 *  Created on: Nov 28, 2022
 *      Author: amr
 */

#ifndef KENREFRUNNER_H_
#define KENREFRUNNER_H_

#include "../config/KEnRefConfig.h"
#include <filesystem>


class KEnRefRunner {
	std::filesystem::path cwd;
	std::filesystem::path model;

public:
	KEnRefRunner();
	virtual ~KEnRefRunner();
	KEnRefRunner(const KEnRefRunner &other);
	KEnRefRunner(KEnRefRunner &&other);
//	KEnRefRunner& operator=(KEnRefRunner &&other);
//	KEnRefRunner& operator=(const KEnRefRunner &other);
	virtual int run() = 0;
};


#endif /* KENREFRUNNER_H_ */
