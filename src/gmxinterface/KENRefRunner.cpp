/*
 * KENRefRunner.cpp
 *
 *  Created on: Nov 28, 2022
 *      Author: amr
 */

#include "core/kenrefrunner.h"


KEnRefRunner::KEnRefRunner() = default;
KEnRefRunner::~KEnRefRunner() = default;
KEnRefRunner::KEnRefRunner(const KEnRefRunner &other) = default;
KEnRefRunner::KEnRefRunner(KEnRefRunner &&other) noexcept= default;
KEnRefRunner& KEnRefRunner::operator=(KEnRefRunner &&other)  noexcept = default;
KEnRefRunner& KEnRefRunner::operator=(const KEnRefRunner &other) = default;

int KEnRefRunner::run(){
	return 0;
}
