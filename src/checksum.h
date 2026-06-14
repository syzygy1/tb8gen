/*
  Copyright (c) 2011-2013 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef CHECKSUM_H
#define CHECKSUM_H

void add_cityhash(const char *name);
void add_xxhash(const char *name);
void verify_checksum(const char *name, bool new);
void print_checksum(const char *name, char *sum);

#endif

