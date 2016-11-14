/*
 * KMerCounterUtils.h
 *
 *  Created on: Nov 11, 2016
 *      Author: jayanga
 */

#pragma once

#ifdef DEBUG_BUILD
#define DEBUG(k, v) do { printf("%s %"PRIu64"\n", k, v); } while (0)
#else
#define DEBUG(k, v)
#endif
