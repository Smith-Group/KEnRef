/*
 * GmxKEnRefRunner.h
 *
 *  Created on: Feb 14, 2023
 *      Author: amr
 */

#ifndef GMXKENREFRUNNER_H_
#define GMXKENREFRUNNER_H_

#include "../core/kenrefrunner.h"

class GmxKEnRefRunner: public KEnRefRunner {
public:
	GmxKEnRefRunner();
	~GmxKEnRefRunner() override;
	GmxKEnRefRunner(const GmxKEnRefRunner &other);
	GmxKEnRefRunner(GmxKEnRefRunner &&other) noexcept;
	int run() override;
	GmxKEnRefRunner& operator=(const GmxKEnRefRunner &other);
	GmxKEnRefRunner& operator=(GmxKEnRefRunner &&other) noexcept;
};

#endif /* GMXKENREFRUNNER_H_ */
