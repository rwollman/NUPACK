%{
#include "parsestruc.h"
#include "pathway_parser.h" 
#include "pathway_lexer.h"

int LINE_NUMBER = 1;
int yyerror(YYLTYPE * yyllocp, value_struc_t ** val, yyscan_t scanner, const char * msg);
%}

%code requires {
#ifndef YY_TYPEDEF_YY_SCANNER_T
#define YY_TYPEDEF_YY_SCANNER_T
typedef void* yyscan_t;
#endif
}

%define api.value.type {value_struc_t *}
%define parse.error verbose
%define api.pure  full
%define parse.lac full
%define locations
%locations
%lex-param {yyscan_t scanner} 
%parse-param {value_struc_t ** expression}
%parse-param {yyscan_t scanner}

%token NEWLINE

%token TOK_PERCENT
%token EQUALS
%token COMMA
%token PERIOD
%token LBRACKET
%token RBRACKET
%token LBRACE
%token RBRACE


%token PLUS
%token MINUS
%token LPAREN
%token RPAREN
%token DUPLEX
%token UNPAIRED

%token TOK_NAME
%token TOK_RESERVED
%token TOK_STRUCTURE
%token TOK_DOMAIN
%token TOK_STRAND
%token TOK_TUBE
%token TOK_STAR
%token TOK_INTEGER
%token TOK_FLOAT
%token TOK_NUC
%token TOK_SEED
%token TOK_MAXSIZE
%token TOK_MATERIAL
%token TOK_SODIUM
%token TOK_MAGNESIUM
%token TOK_TEMPERATURE
%token TOK_MBAD
%token TOK_MREOPT
%token TOK_MRESEED
%token TOK_FSPLIT
%token TOK_NSPLIT
%token TOK_HSPLIT
%token TOK_DANGLES
%token TOK_DGCLAMP
%token TOK_CONCDEF
%token TOK_STOPDEF
%token TOK_GLOBAL_STOPDEF
%token TOK_TRUE
%token TOK_FALSE
%token TOK_PRINTLEAVES
%token TOK_PRINTSTEPS
%token TOK_REDECOMPOSE
%token TOK_FSTRINGENT
%token TOK_FREDECOMP
%token TOK_FREFOCUS
%token TOK_FPASSIVE
%token TOK_GC_INIT
%token TOK_ALLOWWOBBLE
%token TOK_ALLOWMISMATCH
%token TOK_DISABLEMUTWEIGHTS
%token TOK_MINPAIR
%token TOK_TRIALS
%token TOK_OPTTIME
%token TOK_PREVENT
%token TOK_LIBRARY
%token TOK_LIBSEQ
%token TOK_POPULATION
%token TOK_SYMMETRY_MIN
%token TOK_WORDSIZE

%token TOK_SOURCE
%token TOK_WINDOW
%token TOK_EXCLUDE
%token TOK_ACCESSIBLE
%token TOK_SIMILARITY

%token TOK_COMPLEMENTARY
%token TOK_IDENTICAL
%token TOK_MATCH
%token TOK_MATCHRANGE

%token TOK_WEIGHT

%token TOK_COMPLEX
%token TOK_OFFTARGETS
%%

input
    : definitions 
    {
        *expression = $1;
        $$ = $1;
    }
    ;

definitions
    : NEWLINE {
        LINE_NUMBER += 1;
        $$ = np_tt_make_string("");
    }
    | definition NEWLINE {
        LINE_NUMBER += 1; 
        $$ = $1;
    }
    | definitions definition NEWLINE {
        LINE_NUMBER += 1;
        $$ = np_tt_append_value($1, $2);
    }
    | definitions NEWLINE {
        LINE_NUMBER += 1;
        $$ = np_tt_append_value($1, np_tt_make_string(""));
    }
    ;

definition
    : structure_def
    | domain_def
    | tube_def
    | global_prop_def
    | dot_prop_def
    | library_def
    | domainlist_define_type
    | source_def
    | namelist_relation_def
    | match_def
    | offtarget_def
    ;

structure_def
    : TOK_STRUCTURE domain_name EQUALS structure
    {
        $$ = np_tt_make_definition(np_tt_make_int(TOK_STRUCTURE), $2, $4);
    }
    | domain_name PERIOD TOK_STRUCTURE EQUALS structure
    {
        $$ = np_tt_make_definition(np_tt_make_int(TOK_STRUCTURE), $1, $5);
    }
    ;

structure
    : dpp_struc_list 
    | dup_struc_list 
    ;

dpp_struc_list
    : dpp_struc_list dpp_struc_el 
    {
        np_tt_append_string($1, $2->strval);
        np_tt_destroy_value_struc($2);
        $$ = $1;
    }
    | dpp_struc_el
    ;

dpp_struc_el
    : PLUS 
    {
        $$ = np_tt_make_string("+");
    }
    | LPAREN 
    {
        $$ = np_tt_make_string("(");
    }
    | RPAREN 
    {
        $$ = np_tt_make_string(")");
    }
    | PERIOD
    {
        $$ = np_tt_make_string(".");
    }
    ;

dup_struc_list
    : dup_struc_list PLUS dup_struc
    {
        np_tt_append_string($1, "+");
        np_tt_append_string($1, $3->strval);
        np_tt_destroy_value_struc($3);
        $$ = $1;
    }
    | dup_struc
    ;

empt_dup_struc_list
    : dup_struc
    | PLUS
    {
        $$ = np_tt_make_string("+");
    }
    | PLUS dup_struc
    {
        $$ = np_tt_make_string("+");
        np_tt_append_string($$, $2->strval);
        np_tt_destroy_value_struc($2);
    }
    | dup_struc PLUS
    {
        np_tt_append_string($1, "+");
        $$ = $1;
    }
    | dup_struc PLUS dup_struc
    {
        np_tt_append_string($1, "+");
        np_tt_append_string($1, $3->strval);
        np_tt_destroy_value_struc($3);
        $$ = $1;
    }
    | dup_struc PLUS dup_struc PLUS
    {
        yyerror(&yylloc, NULL, scanner, 
            "Disconnected structure");
        YYERROR;
    }

    ;

dup_struc
    : dup_struc dup_single_el
    {
        np_tt_append_string($1, $2->strval);
        np_tt_destroy_value_struc($2);
        $$ = $1;
    }
    | dup_single_el ;

dup_single_el
    : duplex_el
    | unpaired_el
    ;

dup_inside_duplex
    : duplex_el
    | unpaired_el
    | PLUS
    { $$ = np_tt_make_string("+"); }
    | LPAREN empt_dup_struc_list RPAREN
    { $$ = $2; }
    ;

unpaired_el
    : UNPAIRED TOK_INTEGER
    {
        $$ = np_tt_make_string("");
        for (int el_count = 0; el_count < $2->intval; el_count++) {
            np_tt_append_string($$, ".");
        }
        np_tt_destroy_value_struc($2);
    }
    ;

duplex_el
    : DUPLEX TOK_INTEGER dup_inside_duplex
    {
        $$ = np_tt_make_string("");
        for (int el_count = 0; el_count < $2->intval; el_count++) {
            np_tt_append_string($$, "(");
        }
        np_tt_append_string($$, $3->strval);
        for (int el_count = 0; el_count < $2->intval; el_count++) {
            np_tt_append_string($$, ")");
        }
        
        np_tt_destroy_value_struc($2);
        np_tt_destroy_value_struc($3);
    }

library_def
    : TOK_LIBRARY units domain_name EQUALS seqlist
    {

        value_struc_t * lib_int = np_tt_make_int(TOK_LIBRARY);
        np_tt_append_value(lib_int, $2);

        $$ = np_tt_make_definition(lib_int, $3, $5);
    }
    ;

domain_def
    : TOK_DOMAIN domain_name EQUALS domain 
    {
        $$ = np_tt_make_definition(np_tt_make_int(TOK_DOMAIN),
                $2, $4);
    }
    ;

domain
    : domain_el
    {
        $$ = $1;
    }
    | domain domain_el
    {
        np_tt_append_string($1, $2->strval);
        np_tt_destroy_value_struc($2);
        $$ = $1;
    }
    ;

domain_el
    : TOK_NUC
    | TOK_NUC TOK_INTEGER 
    {
        int repeat = $2->intval;
        np_tt_destroy_value_struc($2);
        char * curstr = np_tt_strdup($1->strval);
        for (int i = 0; i < repeat - 1; i++) {
            np_tt_append_string($1, curstr);
        }
        free(curstr);
        $$ = $1;
    }
    ;

domain_list
    : domain_name
    | domain_list domain_name
    {
        $$ = np_tt_append_value($1, $2);
    }
    ;

domain_name
    : TOK_NAME
    | TOK_NAME TOK_STAR
    {
        $$ = np_tt_append_string($1, "*");
    }
    ;

tube_def
    : TOK_TUBE TOK_NAME EQUALS name_list
    {
        $$ = np_tt_make_definition(np_tt_make_int(TOK_TUBE),
                $2, $4);
    }
    ;

name_list
    : domain_name
    | name_list domain_name
    {
        $$ = np_tt_append_value($1, $2);
    }
    ;

dot_prop_def
    : numerical_property
    | namelist_property
    | rangelist_property
    ;

rangelist_property
    : dot_namelist rangelist_property_name units name_or_blank EQUALS rangelist
    { 
        np_tt_append_value($2, $3);
        np_tt_append_value($2, $4);
        $$ = np_tt_make_definition($2, $1, $6);
    }
    ;

name_or_blank
    : domain_name
    | %empty
    { $$ = NULL; };

rangelist
    : range
    | rangelist range
    { 
        np_tt_append_value($1, $2); 
        $$ = $1;
    }
    

range
    : LBRACKET TOK_FLOAT COMMA TOK_FLOAT RBRACKET
    {
        np_tt_append_value($2, $4);
        $$ = $2;
    } 
    | TOK_FLOAT
    {
        value_struc_t * tmp_value = np_tt_make_double(($1)->doubleval);
        np_tt_append_value($1, tmp_value);
        $$ = $1;
    }

rangelist_property_name
    : TOK_EXCLUDE
    { $$ = np_tt_make_int(TOK_EXCLUDE); }
    | TOK_ACCESSIBLE
    { $$ = np_tt_make_int(TOK_ACCESSIBLE); }
    | TOK_SIMILARITY
    { $$ = np_tt_make_int(TOK_MATCHRANGE); }
/*    | TOK_MATCHRANGE
    { $$ = np_tt_make_int(TOK_MATCHRANGE); } */
    | TOK_WORDSIZE
    { $$ = np_tt_make_int(TOK_WORDSIZE); }
    ;


namelist_property
    : spacelist_property
    | commalist_property
    ;

spacelist_property
    : dot_namelist spacelist_property_name units EQUALS namelist
    {
        np_tt_append_value($2, $3);
        $$ = np_tt_make_definition($2, $1, $5);
    }
    ;
    
commalist_property
    : dot_namelist commalist_property_name EQUALS commalist
    {
        $$ = np_tt_make_definition($2, $1, $4);
    }
    ;

spacelist_property_name
    : TOK_LIBSEQ
    { $$ = np_tt_make_int(TOK_LIBSEQ); }
    ;

commalist_property_name
    : TOK_SOURCE
    { $$ = np_tt_make_int(TOK_SOURCE); }
    ;

numerical_property
    : dot_namelist numerical_property_name units EQUALS TOK_FLOAT
    {
        np_tt_append_value($2, $3);
        $$ = np_tt_make_definition($2, $1, $5);
    }
    ;
    
numerical_property_name
    : TOK_CONCDEF
    { $$ = np_tt_make_int(TOK_CONCDEF); }
    | TOK_MAXSIZE
    { $$ = np_tt_make_int(TOK_MAXSIZE); }
    | TOK_STOPDEF
    { $$ = np_tt_make_int(TOK_STOPDEF); }
    | TOK_WEIGHT
    { $$ = np_tt_make_int(TOK_WEIGHT); }
    ;


dot_namelist
    : domain_name PERIOD
    | dot_namelist domain_name PERIOD
    {
        $$ = np_tt_append_value($1, $3);
    }
    ;

namelist
    : domain_name
    | namelist domain_name
    {
        $$ = np_tt_append_value($1, $2);
    }
    ;
    
commalist
    : domain_name
    | commalist COMMA domain_name
    {
        $$ = np_tt_append_value($1, $3);
    }
    ;

domainlist_define_type
    : domainlist_type domain_name EQUALS domain_list
    { $$ = np_tt_make_definition($1, $2, $4); }
    ;

domainlist_type
    : TOK_SYMMETRY_MIN
    { $$ = np_tt_make_int(TOK_SYMMETRY_MIN); }
    | TOK_STRAND
    { $$ = np_tt_make_int(TOK_STRAND); }
    | TOK_WINDOW
    { $$ = np_tt_make_int(TOK_WINDOW); }
    | TOK_COMPLEX
    { $$ = np_tt_make_int(TOK_COMPLEX); }
    ;



global_prop_def
    : numerical_global
    | boolean_global
    | namelist_global
    | seqlist_global
    ;

namelist_global
    : namelist_global_name EQUALS namelist
    {
        $$ = np_tt_make_definition($1,
                NULL, $3);
    }
    ;

namelist_global_name
    : TOK_MATERIAL
    { $$ = np_tt_make_int(TOK_MATERIAL); }
    | TOK_DANGLES
    { $$ = np_tt_make_int(TOK_DANGLES); }
    ;

seqlist_global
    : seqlist_global_name EQUALS seqlist
    { $$ = np_tt_make_definition($1, NULL, $3); }
    | seqlist_global_name commalist EQUALS seqlist
    { $$ = np_tt_make_definition($1, $2, $4); }
    ;

seqlist_global_name
    : TOK_PREVENT 
    { $$ = np_tt_make_int(TOK_PREVENT); }
    ;

numerical_global
    : numerical_global_name units EQUALS TOK_FLOAT
    {
        
        np_tt_append_value($1, $2);
        $$ = np_tt_make_definition($1,
                NULL, $4);
    }
    ;

numerical_global_name
    : TOK_MAGNESIUM
    { $$ = np_tt_make_int(TOK_MAGNESIUM); }
    | TOK_STOPDEF
    { $$ = np_tt_make_int(TOK_GLOBAL_STOPDEF); }
    | TOK_SODIUM
    { $$ = np_tt_make_int(TOK_SODIUM); }
    | TOK_TEMPERATURE
    { $$ = np_tt_make_int(TOK_TEMPERATURE); }
    | TOK_SEED
    { $$ = np_tt_make_int(TOK_SEED); }
    | TOK_MBAD
    { $$ = np_tt_make_int(TOK_MBAD); }
    | TOK_MREOPT
    { $$ = np_tt_make_int(TOK_MREOPT); }
    | TOK_MRESEED
    { $$ = np_tt_make_int(TOK_MRESEED); }
    | TOK_FSPLIT
    { $$ = np_tt_make_int(TOK_FSPLIT); }
    | TOK_NSPLIT
    { $$ = np_tt_make_int(TOK_NSPLIT); }
    | TOK_HSPLIT
    { $$ = np_tt_make_int(TOK_HSPLIT); }
    | TOK_DGCLAMP
    { $$ = np_tt_make_int(TOK_DGCLAMP); }
    | TOK_FSTRINGENT
    { $$ = np_tt_make_int(TOK_FSTRINGENT); }
    | TOK_FREDECOMP
    { $$ = np_tt_make_int(TOK_FREDECOMP); }
    | TOK_FREFOCUS
    { $$ = np_tt_make_int(TOK_FREFOCUS); }
    | TOK_FPASSIVE
    { $$ = np_tt_make_int(TOK_FPASSIVE); }
    | TOK_GC_INIT
    { $$ = np_tt_make_int(TOK_GC_INIT); }
    | TOK_MINPAIR
    { $$ = np_tt_make_int(TOK_MINPAIR); }
    | TOK_TRIALS
    { $$ = np_tt_make_int(TOK_TRIALS); }
    | TOK_OPTTIME
    { $$ = np_tt_make_int(TOK_OPTTIME); }
    | TOK_POPULATION
    { $$ = np_tt_make_int(TOK_POPULATION); }
    | TOK_PRINTSTEPS
    { $$ = np_tt_make_int(TOK_PRINTSTEPS); }
    ;

seqlist
    : domain
    | seqlist COMMA domain
    {
        $$ = np_tt_append_value($1, $3);
    }
    ;

boolean_global
    : boolean_global_name EQUALS boolean
    {
        $$ = np_tt_make_definition($1, NULL, $3);
    };

boolean_global_name
    : TOK_PRINTLEAVES
    { $$ = np_tt_make_int(TOK_PRINTLEAVES); }
    | TOK_REDECOMPOSE
    { $$ = np_tt_make_int(TOK_REDECOMPOSE); }
    | TOK_ALLOWWOBBLE
    { $$ = np_tt_make_int(TOK_ALLOWWOBBLE); }
    | TOK_ALLOWMISMATCH
    { $$ = np_tt_make_int(TOK_ALLOWMISMATCH); }
    | TOK_DISABLEMUTWEIGHTS
    { $$ = np_tt_make_int(TOK_DISABLEMUTWEIGHTS); }
    ;

boolean
    : TOK_TRUE
    { $$ = np_tt_make_int(1); }
    | TOK_FALSE
    { $$ = np_tt_make_int(0); }
    ;

units
    : %empty
    { $$ = NULL; }
    | LBRACKET TOK_NAME RBRACKET
    { $$ = $2 ; }
    | LBRACKET TOK_PERCENT RBRACKET
    { $$ = np_tt_make_string("%"); }
    ;


source_def
    : TOK_SOURCE TOK_NAME EQUALS domain
    { $$ = np_tt_make_definition(np_tt_make_int(TOK_SOURCE), $2 , $4); }
    ;

namelist_relation_def
    : namelist_relation_type units namelist EQUALS namelist
    { 
        np_tt_append_value($1, $2);
        $$ = np_tt_make_definition($1, $3, $5); 
    }
    ;

namelist_relation_type
    : TOK_COMPLEMENTARY
    { $$ = np_tt_make_int(TOK_COMPLEMENTARY); }
    | TOK_MATCH
    { $$ = np_tt_make_int(TOK_IDENTICAL); }
    ;

match_def
    : TOK_SIMILARITY TOK_NAME EQUALS domain 
    { $$ = np_tt_make_definition(np_tt_make_int(TOK_MATCH), $2, $4); }
    ;


offtarget_def
    : domain_name PERIOD TOK_OFFTARGETS EQUALS offtargets
    { $$ = np_tt_make_definition(np_tt_make_int(TOK_OFFTARGETS), $1, $5); }
    ;

superlist
    : namelist
    { 
      value_struc_t * inner = np_tt_make_list();
      value_struc_t * outer = np_tt_make_list();
      np_tt_list_append(inner, $1);
      $$ = np_tt_list_append(outer, inner);
    }
    | superlist COMMA namelist
    { 
      value_struc_t * inner = np_tt_make_list();
      np_tt_list_append(inner, $3);
      $$ = np_tt_list_append($1, inner);
    }
    ;

set
    : LBRACE superlist RBRACE
    { $$ = $2; }
    ;
    
offtargets
    : maxsize 
    | maxsize PLUS set
    { $$ = np_tt_append_value($1, $3); }
    | maxsize MINUS set
    { 
      np_tt_append_value($1, np_tt_make_list());
      $$ = np_tt_append_value($1, $3);
    }
    | maxsize PLUS set MINUS set
    { 
      np_tt_append_value($1, $3);
      $$ = np_tt_append_value($1, $5);
    }
    | set 
    {
      $$ = np_tt_append_value(np_tt_make_int(0), $1);
    }
    | set MINUS set
    {
      value_struc_t * tmp = np_tt_make_int(0);
      np_tt_append_value(tmp, $1);
      $$ = np_tt_append_value(tmp, $3);
    }
    ;

maxsize
    : LBRACE TOK_MAXSIZE EQUALS TOK_FLOAT RBRACE
    { $$ = $4; }
    ;

%%

int yyerror(YYLTYPE * yyllocp, value_struc_t ** val, yyscan_t scanner, const char * msg) {
  if (yyllocp->first_line) {
    fprintf(stderr, "Line %d\n", LINE_NUMBER);
    fprintf(stderr, "%d.%d-%d.%d: error: ", 
      yyllocp->first_line, yyllocp->first_column,
      yyllocp->last_line, yyllocp->last_column);
  } else {
    fprintf(stderr, "Error: %s", msg);
  }
  fprintf(stderr, "%s", msg);
  fprintf(stderr, "\n");
  return 0;
} 

static value_struc_t * make_val_struc(void) {
  value_struc_t * val = malloc(sizeof(value_struc_t));
  val->type = NP_TT_INTEGER;
  val->intval = 0;
  val->doubleval = 0.0;
  val->list = NULL;
  val->strval = NULL;
  val->strlen = 0;
  val->strcap = 0;
  val->stype = NULL;
  val->id = NULL;
  val->def = NULL;
  val->next = NULL;
  return val;
}

char* np_tt_strdup(const char* str)
{
  size_t len;
  char* copy;

  len = strlen(str) + 1;
  copy = (char*) malloc(len);
  if (!copy)  return 0;
  memcpy(copy, str, len);
  return copy;
}

value_struc_t * np_tt_make_int(int input) {
  value_struc_t * val = make_val_struc();
  val->strval = NULL;
  val->type = NP_TT_INTEGER;
  val->intval = input;
  val->doubleval = (double) input;
  return val;
}

value_struc_t * np_tt_make_double(double input) {
  value_struc_t * val = make_val_struc();
  val->strval = NULL;
  val->type = NP_TT_FLOAT;
  val->doubleval = input;
  val->intval = (int) input;
  return val;
}

value_struc_t * np_tt_make_string(char * input) {
  value_struc_t * val = make_val_struc();
  int l = strlen(input);
  val->strval = (char*) malloc((l + 1) * sizeof(char));
  strncpy(val->strval, input, l + 1);
  val->type = NP_TT_STRING;
  val->strlen = l + 1;
  val->strcap = l + 1;
  return val;
}

value_struc_t * np_tt_append_string(value_struc_t * cur, 
    char * next) {
  int l = strlen(next);
  char * tempchar = cur->strval;
  char * nextdup = np_tt_strdup(next);
  if (l + cur->strlen + 1 > cur->strcap) {
    tempchar = (char*) realloc(cur->strval, 
          2 * (l + cur->strlen + 1));
    check_mem(tempchar);
    cur->strval = tempchar;
    cur->strcap = 2 * (l + cur->strlen + 1);
  }
  strncat(tempchar, nextdup, l+1);
  cur->strlen = l + cur->strlen + 1;
  free(nextdup);

  return cur;
error:
  cur->strval = NULL;
  return NULL;
}

value_struc_t * np_tt_make_list() {
  value_struc_t * val = make_val_struc();
  val->type = NP_TT_LIST;  
  return val;
}

value_struc_t * np_tt_make_definition(value_struc_t * stype, value_struc_t * id,
    value_struc_t * def) {
  value_struc_t * val = make_val_struc();
  val->type = NP_TT_DEFINITION;
  val->stype = stype;
  val->id = id;
  val->def = def;
  return val;
}

value_struc_t * np_tt_append_value(value_struc_t * cur, value_struc_t * app) {
  value_struc_t * cv = cur;
  while (cv->next) {cv = cv->next; }
  cv->next = app;
  return cur;
}

value_struc_t * np_tt_list_append(value_struc_t * cur, value_struc_t * app) {
  value_struc_t * tmp = cur->list;
  if (!tmp) {
    cur->list = app;
  } else {
    while (tmp->next) tmp = tmp->next;
    tmp->next = app;
  }
  return cur;  
}


void np_tt_destroy_value_struc(value_struc_t * val) {
  value_struc_t * cur = val;
  value_struc_t * prev = cur;

  while (cur) {
    np_tt_destroy_value_struc(cur->list);
    np_tt_destroy_value_struc(cur->stype);
    np_tt_destroy_value_struc(cur->id);
    np_tt_destroy_value_struc(cur->def);
    if (cur->strval && cur->strcap > 0) {
      free(cur->strval);
    }
    prev = cur;
    cur = cur->next;
    free(prev);
  }
}

void np_tt_print_ast(value_struc_t * root, int indent) {
  value_struc_t * cur = root;
  
  char * indent_string = (char *) malloc((indent + 1) * sizeof(char));
  for (int i = 0; i < indent; i++) {
    indent_string[i] = ' ';
  }
  indent_string[indent] = '\0';
  
  while (cur) {
    if (cur->type == NP_TT_DEFINITION) {
      printf("%sDEFINITION:\n", indent_string);
      printf("%s  stype:\n", indent_string);
      np_tt_print_ast(cur->stype, indent + 4);
      printf("%s  id:\n", indent_string);
      np_tt_print_ast(cur->id, indent + 4);
      printf("%s  def:\n", indent_string);
      np_tt_print_ast(cur->def, indent + 4);
    } else if (cur->type == NP_TT_STRING) {
      printf("%s%s\n", indent_string, cur->strval);
    } else if (cur->type == NP_TT_FLOAT) {
      printf("%s%lf\n", indent_string, cur->doubleval);
    } else if (cur->type == NP_TT_INTEGER) {
      printf("%s%i\n", indent_string, cur->intval);
    }
    cur = cur->next;
  }
  free(indent_string);
}


value_struc_t * get_ast(const char * parse_string) {
  value_struc_t * vals;
  yyscan_t scanner;
  YY_BUFFER_STATE state;

  if (yylex_init(&scanner)) {
    return NULL;
  }

  state = yy_scan_string(parse_string, scanner);

  vals = NULL;
  if (yyparse(&vals, scanner)) {
    return NULL;
  }

  yy_delete_buffer(state, scanner);

  yylex_destroy(scanner);

  return vals;
}
