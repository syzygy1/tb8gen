#ifndef TB8GEN_H
#define TB8GEN_H

extern bool one_sided, wins_only;
extern int one_sided_stm;

extern size_t kslice_size, kslice_sub_size[MAX_SETS];
extern char *g_tablename;
extern bool symmetric;
extern Position g_pos;
extern uint64_t *work_g, *work_capt[MAX_SETS];
extern uint64_t g_stats[2][MAX_STATS];
extern bool g_use_rans;

static constexpr int g_num_pawns = 0;

#endif
