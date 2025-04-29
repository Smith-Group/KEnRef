/*
 * config.h
 *
 *  Created on: Nov 29, 2022
 *      Author: amr
 */

#ifndef KENREFCONFIG_H_
#define KENREFCONFIG_H_


/**This is a very bad way of handing multiple configurations.
 * However, it works as a quick and dirty one one the proof of concept now.
 * TODO Later, use CMAke to configure the project
 */
#define GROMACS_WRAPPER true
#ifdef TEST
#	define _TEST_ true
#else
#	define _TEST_ false
#endif

#endif /* KENREFCONFIG_H_ */
