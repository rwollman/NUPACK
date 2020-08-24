#pragma once

// Zed shaw's awesome debug macros

#include <stdio.h>
#include <errno.h>
#include <string.h>

enum ERROR_RETURN
{
  ERR_OK              = 0, // explicit for clarity
  ERR_NOCONVERGE      = 1,
  ERR_OVERFLOW        = 2,
  ERR_BADN            = 3,
  ERR_BADK            = 4,
  ERR_NOINPUT         = 5,
  ERR_HELP            = 6,
  ERR_HELP_OPEN       = 7,
  ERR_NOCMP           = 8,
  ERR_SOLVDENS        = 9,
  ERR_BADFREEENERGY   = 10,
  ERR_BADROWINP       = 11,
  ERR_CON             = 12,
  ERR_LOG             = 13,
  ERR_EQ              = 14,
  ERR_OOM             = 15,
  ERR_INITIAL         = 16,
  ERR_INVALID_STATE   = 17,
  ERR_INVALID_INPUT   = 18,
  ERR_OTHER           = 19
};

enum NUPACK_DESIGN_STEP
{
  NP_STEP_NOOPT = 0,
  NP_STEP_LEAFOPT=1,
  NP_STEP_TREEOPT=2,
  NP_STEP_TUBEOPT=3,
  NP_STEP_COUNT
};

// old, generic debugging macros
#ifdef NDEBUG
#define debug(M, ...)
#else
#define debug(M, ...) fprintf(stderr, "DEBUG %s:%d: " M "\n", __FILE__, __LINE__, ##__VA_ARGS__)
#endif

#define clean_errno() (errno == 0 ? "None" : strerror(errno))

#define log_err(M, ...) fprintf(stderr, "[ERROR] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, clean_errno(), ##__VA_ARGS__)

#define log_warn(M, ...) fprintf(stderr, "[WARN] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, clean_errno(), ##__VA_ARGS__)

#define log_info(M, ...) fprintf(stderr, "[INFO] (%s:%d) " M "\n", __FILE__, __LINE__, ##__VA_ARGS__)

#define zed_check(A, M, ...) if(!(A)) { log_err(M, ##__VA_ARGS__); errno=0; goto error; }

#define sentinel(M, ...)  { log_err(M, ##__VA_ARGS__); errno=0; goto error; }

#define check_mem(A) zed_check((A), "Out of memory.")

#ifdef NDEBUG
#define check_debug(A, M, ...)
#else
#define check_debug(A, M, ...) if(!(A)) { debug(M, ##__VA_ARGS__); errno=0; goto error; }
#endif

