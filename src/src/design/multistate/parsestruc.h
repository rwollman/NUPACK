#pragma once

#include "pathway_debug.h"

#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Initialize LOC. */
#define LOCATION_RESET(Loc)                  \
  (Loc).first_column = (Loc).first_line = 1;  \
  (Loc).last_column =  (Loc).last_line = 1;

/* Advance of NUM lines. */
#define LOCATION_LINES(Loc, Num)             \
  (Loc).last_column = 1;                      \
  (Loc).last_line += Num;

/* Restart: move the first cursor to the last position. */
#define LOCATION_STEP(Loc)                   \
  (Loc).first_column = (Loc).last_column;     \
  (Loc).first_line = (Loc).last_line;

/* Output LOC on the stream OUT. */
#define LOCATION_PRINT(Out, Loc)                               \
  if ((Loc).first_line != (Loc).last_line)                      \
    fprintf (Out, "%d.%d-%d.%d",                                \
             (Loc).first_line, (Loc).first_column,              \
             (Loc).last_line, (Loc).last_column - 1);           \
  else if ((Loc).first_column < (Loc).last_column - 1)          \
    fprintf (Out, "%d.%d-%d", (Loc).first_line,                 \
             (Loc).first_column, (Loc).last_column - 1);        \
  else                                                          \
    fprintf (Out, "%d.%d", (Loc).first_line, (Loc).first_column)


typedef enum VALUE_TYPE_T_ {
  NP_TT_INTEGER,
  NP_TT_FLOAT,
  NP_TT_NAME,
  NP_TT_STRING,
  NP_TT_DEFINITION,
  NP_TT_GLOBAL,
  NP_TT_PROPERTY,
  NP_TT_LIST
} value_t;

typedef struct VALUE_STRUC_T_ {
   value_t type;
   int intval;
   double doubleval;
   char * strval;
   int strlen;
   int strcap;
   
   struct VALUE_STRUC_T_ * list;
   
   struct VALUE_STRUC_T_ * stype;
   struct VALUE_STRUC_T_ * id;
   struct VALUE_STRUC_T_ * def;
   
   struct VALUE_STRUC_T_ * next;
} value_struc_t;

char * np_tt_strdup(const char * inp);
value_struc_t * np_tt_make_double(double val);
value_struc_t * np_tt_make_int(int val);
value_struc_t * np_tt_make_string(char * str);
value_struc_t * np_tt_append_string(value_struc_t * cur, char * str);
value_struc_t * np_tt_make_list();
value_struc_t * np_tt_make_definition(value_struc_t * stype, 
    value_struc_t * id,  value_struc_t * def);
value_struc_t * np_tt_append_value(value_struc_t * cur, value_struc_t * app);
value_struc_t * np_tt_list_append(value_struc_t * cur, value_struc_t * app);
void np_tt_print_ast(value_struc_t * root, int indent);

void np_tt_destroy_value_struc(value_struc_t * cur);
value_struc_t * get_ast(const char * parse_string);

#ifdef __cplusplus
}
#endif /* __cplusplus */

