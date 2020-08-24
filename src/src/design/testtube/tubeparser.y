%{
#include "parsestruc.h"
#include "npparser.h"
#include "nplexer.h" 

int LINE_NUMBER = 1;
int yyerror(
        YYLTYPE * yyllocp, 
        value_struc_t ** val, 
        yyscan_t scanner, 
        const char * msg);
%}

%code requires {
#define YYSTYPE value_struc_t*
#ifndef YY_TYPEDEF_YY_SCANNER_T
#define YY_TYPEDEF_YY_SCANNER_T
typedef void* yyscan_t;
#endif
}

%define parse.error verbose
%define api.pure  full
%define parse.lac full
%define locations
%locations
%lex-param {yyscan_t scanner} 
%parse-param {value_struc_t ** expression}
%parse-param {yyscan_t scanner}
    /*
    %union{
        value_struc_t * valstruc;
    }
    */

%token NEWLINE

%token TOK_PERCENT
%token EQUALS
%token PERIOD
%token LBRACKET
%token RBRACKET

%token PLUS
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
%token TOK_MUNFAVORABLE
%token TOK_MLEAFOPT
%token TOK_MRESEED
%token TOK_FSPLIT
%token TOK_NSPLIT
%token TOK_HSPLIT
%token TOK_DANGLES
%token TOK_SEQ
%token TOK_CONCDEF
%token TOK_STOPDEF
%token TOK_TRUE
%token TOK_FALSE
%token TOK_PRINTLEAVES
%token TOK_PRINTSTEPS
%token TOK_INCLUDE_ALL
%token TOK_FREDECOMP
%token TOK_FREFOCUS
%token TOK_FSTRINGENT
%token TOK_FPASSIVE
%token TOK_GC_INIT
%token TOK_ALLOWWOBBLE
%token TOK_ALLOWMISMATCH
%token TOK_DISABLEMUTWEIGHTS
%token TOK_MINPAIR
%token TOK_TRIALS
%token TOK_OPTTIME
%token TOK_SINGLE_DECOMP
%token TOK_DGCLAMP
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
    {$$ = $1;}
    | domain_def
    {$$ = $1;}
    | strand_def
    {$$ = $1;}
    | tube_def
    {$$ = $1;}
    | global_prop_def
    {$$ = $1;}
    | dot_prop_def
    {$$ = $1;}
    ;

structure_def
    : TOK_STRUCTURE TOK_NAME EQUALS structure
        {
            $$ = np_tt_make_definition(np_tt_make_int(TOK_STRUCTURE),
                $2, $4);
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
        int el_count = 0;
        for (el_count = 0; el_count < $2->intval; el_count++) {
            np_tt_append_string($$, ".");
        }
        np_tt_destroy_value_struc($2);
    }
    ;

duplex_el
    : DUPLEX TOK_INTEGER dup_inside_duplex
    {
        $$ = np_tt_make_string("");
        int el_count = 0;
        for (el_count = 0; el_count < $2->intval; el_count++) {
            np_tt_append_string($$, "(");
        }
        np_tt_append_string($$, $3->strval);
        for (el_count = 0; el_count < $2->intval; el_count++) {
            np_tt_append_string($$, ")");
        }
        
        np_tt_destroy_value_struc($2);
        np_tt_destroy_value_struc($3);
    }

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
    {
        $$ = $1;
    }
    | TOK_NUC TOK_INTEGER 
    {
        int repeat = $2->intval;
        np_tt_destroy_value_struc($2);
        int i;
        char * curstr = np_tt_strdup($1->strval);
        for (i = 0; i < repeat - 1; i++) {
            np_tt_append_string($1, curstr);
        }
        free(curstr);
        $$ = $1;
    }
    ;

strand_def
    : TOK_STRAND TOK_NAME EQUALS domain_list
    {
        $$ = np_tt_make_definition(np_tt_make_int(TOK_STRAND),
                $2, $4);
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
    : TOK_NAME TOK_STAR
    {
        $$ = np_tt_append_string($1, "*");
    }
    | TOK_NAME
    {
        $$ = $1;
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
    : TOK_NAME
    | name_list TOK_NAME
    {
        $$ = np_tt_append_value($1, $2);
    }
    ;

dot_prop_def
    : numerical_property
    | namelist_property
    /* | structure_property */
    ;
/* 
structure_property
    : dot_namelist structure_property_name EQUALS structure
    {
        $$ = np_tt_make_definition($2, $1, $4);
    }
    ;

structure_property_name
    : TOK_ALTSTRUC
    {$$ = np_tt_make_int(TOK_ALTSTRUC);}   
    ;
*/

numerical_property
    : dot_namelist numerical_property_name units 
        EQUALS TOK_FLOAT
    {
        np_tt_append_value($2, $3);
        $$ = np_tt_make_definition($2,
                $1, $5);
    }
    ;

namelist_property
    : dot_namelist namelist_property_name EQUALS namelist
    {
        $$ = np_tt_make_definition($2,
                $1, $4);
    }
    ;
    
numerical_property_name
    : TOK_CONCDEF
    { $$ = np_tt_make_int(TOK_CONCDEF); }
    | TOK_MAXSIZE
    { $$ = np_tt_make_int(TOK_MAXSIZE); }
    | TOK_STOPDEF
    { $$ = np_tt_make_int(TOK_STOPDEF); }
    ;

namelist_property_name
    : TOK_SEQ
    { $$ = np_tt_make_int(TOK_SEQ);}
    ;

dot_namelist
    : TOK_NAME PERIOD
    | dot_namelist TOK_NAME PERIOD
    {
        $$ = np_tt_append_value($1, $3);
    }
    ;

namelist
    : TOK_NAME
    | namelist TOK_NAME
    {
        $$ = np_tt_append_value($1, $2);
    }
    ;

global_prop_def
    : numerical_global
    | boolean_global
    | namelist_global
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
    { $$ = np_tt_make_int(TOK_MATERIAL);}
    | TOK_DANGLES
    { $$ = np_tt_make_int(TOK_DANGLES);}
    ;

numerical_global
    : numerical_global_name units EQUALS TOK_FLOAT
    {
        $$ = np_tt_make_definition($1,
                NULL, $4);
    }
    ;

numerical_global_name
    : TOK_MAGNESIUM
    { $$ = np_tt_make_int(TOK_MAGNESIUM); }
    | TOK_SODIUM
    { $$ = np_tt_make_int(TOK_SODIUM); }
    | TOK_TEMPERATURE
    { $$ = np_tt_make_int(TOK_TEMPERATURE); }
    | TOK_SEED
    { $$ = np_tt_make_int(TOK_SEED); }
    | TOK_MUNFAVORABLE
    { $$ = np_tt_make_int(TOK_MUNFAVORABLE); }
    | TOK_MLEAFOPT
    { $$ = np_tt_make_int(TOK_MLEAFOPT); }
    | TOK_MRESEED
    { $$ = np_tt_make_int(TOK_MRESEED); }
    | TOK_FSPLIT
    { $$ = np_tt_make_int(TOK_FSPLIT); }
    | TOK_NSPLIT
    { $$ = np_tt_make_int(TOK_NSPLIT); }
    | TOK_HSPLIT
    { $$ = np_tt_make_int(TOK_HSPLIT); }
    | TOK_FREDECOMP
    { $$ = np_tt_make_int(TOK_FREDECOMP); }
    | TOK_FREFOCUS
    { $$ = np_tt_make_int(TOK_FREFOCUS); }
    | TOK_FSTRINGENT
    { $$ = np_tt_make_int(TOK_FSTRINGENT); }
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
    | TOK_DGCLAMP
    { $$ = np_tt_make_int(TOK_DGCLAMP); }
    ;

boolean_global
    : boolean_global_name EQUALS boolean
    {
        $$ = np_tt_make_definition($1,
                NULL, $3);
    };

boolean_global_name
    : TOK_PRINTSTEPS
    { $$ = np_tt_make_int(TOK_PRINTSTEPS); }
    | TOK_PRINTLEAVES
    { $$ = np_tt_make_int(TOK_PRINTLEAVES); }
    | TOK_ALLOWWOBBLE
    { $$ = np_tt_make_int(TOK_ALLOWWOBBLE); }
    | TOK_ALLOWMISMATCH
    { $$ = np_tt_make_int(TOK_ALLOWMISMATCH); }
    | TOK_INCLUDE_ALL
    { $$ = np_tt_make_int(TOK_INCLUDE_ALL); }
    | TOK_DISABLEMUTWEIGHTS
    { $$ = np_tt_make_int(TOK_DISABLEMUTWEIGHTS); }
    | TOK_SINGLE_DECOMP
    { $$ = np_tt_make_int(TOK_SINGLE_DECOMP); }
    ;

boolean
    : TOK_TRUE
    { $$ = np_tt_make_int(1);}
    | TOK_FALSE
    { $$ = np_tt_make_int(0);}
    ;

units
    : 
    { $$ = NULL;}
    | LBRACKET TOK_NAME RBRACKET
    { $$ = $2 ;}
    | LBRACKET TOK_PERCENT RBRACKET
    { $$ = np_tt_make_string("%");}
    ;


%%

int
yyerror(
    YYLTYPE * yyllocp, 
    value_struc_t ** val, 
    yyscan_t scanner, 
    const char * msg)
{
    
    if(yyllocp->first_line) {
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
