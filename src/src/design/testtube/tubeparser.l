
%{
#include "parsestruc.h"
#include "npparser.h"
/* #define YY_USER_ACTION yylloc.first_line = yylloc.last_line = yylineno; */
#define YY_USER_ACTION yylloc->last_column += yyleng;
int yyerror(YYLTYPE * yyllocp, value_struc_t ** val, yyscan_t scanner,
            const char *s);

static char RESERVED_NAMES[][20] = {
"seed",
"maxsize",
"material",
"sodium",
"magnesium",
"temperature",
"mbad",
"mreopt",
"mreseed",
"fsplit",
"nsplit",
"hsplit",
"dangles",
"structure",
"tube",
"domain",
"strand",
"seq",
"conc",
"stop",
"false",
"true",
"printleaves",
"printsteps",
"fredecomp",
"frefocus",
"fstringent",
"fpassive",
"initgc",
"allowwobble",
"allowmismatch",
"nomutweights",
"minppairsaved",
"trials",
"maxopttime",
"includeall",
"singledecomp",
"dgclamp",
};

static int RESERVED_IDS[] = {
TOK_SEED,           
TOK_MAXSIZE,        
TOK_MATERIAL,       
TOK_SODIUM,         
TOK_MAGNESIUM,      
TOK_TEMPERATURE,    
TOK_MUNFAVORABLE,   
TOK_MLEAFOPT,       
TOK_MRESEED,       
TOK_FSPLIT,        
TOK_NSPLIT,         
TOK_HSPLIT,         
TOK_DANGLES,        
TOK_STRUCTURE,      
TOK_TUBE,           
TOK_DOMAIN,         
TOK_STRAND,         
TOK_SEQ,         
TOK_CONCDEF,        
TOK_STOPDEF,        
TOK_FALSE,          
TOK_TRUE,           
TOK_PRINTLEAVES,
TOK_PRINTSTEPS,     
TOK_FREDECOMP,
TOK_FREFOCUS,
TOK_FSTRINGENT,
TOK_FPASSIVE,
TOK_GC_INIT,
TOK_ALLOWWOBBLE,
TOK_ALLOWMISMATCH,
TOK_DISABLEMUTWEIGHTS,
TOK_MINPAIR,
TOK_TRIALS,
TOK_OPTTIME,
TOK_INCLUDE_ALL,
TOK_SINGLE_DECOMP,
TOK_DGCLAMP,
};

static int NUM_RESERVED = 38;

static int lookup_reserved(char * name);
%}

%option noyywrap never-interactive nounistd nodefault
%option yylineno nounput noinput
%option reentrant bison-bridge bison-locations

%x STRUCTUREMODE SEQUENCEMODE COMMENTMODE
%s PRESEQUENCE PRESTRUCTURE
%%
%{
%}
<*>[ \t\r]     LOCATION_STEP (*yylloc);
<*>\n          {
    LOCATION_LINES (*yylloc, yyleng); 
    LOCATION_STEP (*yylloc); 
    BEGIN(INITIAL);
    return NEWLINE;
}
\% {return TOK_PERCENT; }
[-]?[0-9]*(\.)?[0-9]+([eE][-+]?[0-9]+)? {
    *yylval = np_tt_make_double(atof(yytext));
    return TOK_FLOAT;
}
<*># {BEGIN(COMMENTMODE);}
<COMMENTMODE>. {}
[[] {return LBRACKET;}
[\]] {return RBRACKET;}
<PRESTRUCTURE>= {BEGIN(STRUCTUREMODE); return EQUALS;}
<STRUCTUREMODE>\+ {return PLUS;}
<STRUCTUREMODE>[(] {return LPAREN;}
<STRUCTUREMODE>[)] {return RPAREN;}
<STRUCTUREMODE>[D] {return DUPLEX;}
<STRUCTUREMODE>[U] {return UNPAIRED;}
<PRESEQUENCE>= {BEGIN(SEQUENCEMODE); return EQUALS;}
<SEQUENCEMODE>[acgutmrwsykvhdbnACGUTMRWSYKVHDBN] {
    *yylval = np_tt_make_string(yytext);
    return TOK_NUC;
}
= {return EQUALS;}
<STRUCTUREMODE,SEQUENCEMODE>[0-9]+ {
        int tempint = strtol(yytext, NULL, 10);
        *yylval = np_tt_make_int(tempint);
        return TOK_INTEGER;
    }
\* return TOK_STAR;
<INITIAL,STRUCTUREMODE>\.  {return PERIOD;}
[a-zA-Z_][a-zA-Z0-9_]* {
    int lookupval = lookup_reserved(yytext);
    if (lookupval == TOK_STRUCTURE /* || lookupval == TOK_ALTSTRUC*/) {
        BEGIN(PRESTRUCTURE);
    } else if (lookupval == TOK_DOMAIN) {
        BEGIN(PRESEQUENCE);
    }

    if (lookupval == TOK_NAME) {
        *yylval = np_tt_make_string(yytext);
    }   
    return lookupval;
    }
<*>.|\n fprintf(stderr, "Bad character %c\n", yytext[0]);

%%

static int lookup_reserved(char* name) {
    int rval = TOK_NAME;
    int i, n;
    n = NUM_RESERVED;
    for (i = 0; i < n; i++) {
        if (0 == strncmp(name, RESERVED_NAMES[i], 18)) {
            rval = RESERVED_IDS[i];
            break;
        }
    }   
    return rval;
}

