/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef VERIFY_H
#define VERIFY_H

#include <inttypes.h>

#include "defs.h"
#include "types.h"

extern uint64_t sub_cnt[2][5];

void verify(void);
void init_verification_work(void);
//void delete_intermediate_slices(void);
//void cleanup_verification(void);

#endif
