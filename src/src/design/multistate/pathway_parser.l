
%{
#include "parsestruc.h"
#include "pathway_parser.h"
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
"fsplit",
"nsplit",
"hsplit",
"dangles",
"structure",
"tube",
"domain",
"strand",
"conc",
"stop",
"false",
"true",
"fstringent",
"fredecomp",
"frefocus",
"fpassive",
"allowwobble",
"pairscutoff",
"trials",
"maxopttime",
"prevent",
"library",
"libseq",
"source",
"window",
"similarity",
"complement",
"match",
"weight",
"complex",
"offtargets",
"mreseed",
"dgclamp"
};

static int RESERVED_IDS[] = {
TOK_SEED,           
TOK_MAXSIZE,        
TOK_MATERIAL,       
TOK_SODIUM,         
TOK_MAGNESIUM,      
TOK_TEMPERATURE,    
TOK_MBAD,   
TOK_MREOPT,       
TOK_FSPLIT,        
TOK_NSPLIT,         
TOK_HSPLIT,         
TOK_DANGLES,        
TOK_STRUCTURE,      
TOK_TUBE,           
TOK_DOMAIN,         
TOK_STRAND,         
TOK_CONCDEF,        
TOK_STOPDEF,        
TOK_FALSE,          
TOK_TRUE,           
TOK_FSTRINGENT,
TOK_FREDECOMP,
TOK_FREFOCUS,
TOK_FPASSIVE,
TOK_ALLOWWOBBLE,
TOK_MINPAIR,
TOK_TRIALS,
TOK_OPTTIME,
TOK_PREVENT,
TOK_LIBRARY,
TOK_LIBSEQ,
TOK_SOURCE,         // External sequence constraints
TOK_WINDOW,
TOK_SIMILARITY,
TOK_COMPLEMENTARY,
TOK_MATCH,
TOK_WEIGHT,
TOK_COMPLEX,
TOK_OFFTARGETS,
TOK_MRESEED,
TOK_DGCLAMP
};

static int NUM_RESERVED = sizeof(RESERVED_IDS) / sizeof(RESERVED_IDS[0]);

static int lookup_reserved(char * name);
%}

%option noyywrap never-interactive nounistd nodefault
%option yylineno nounput noinput
%option reentrant bison-bridge bison-locations

%x STRUCTUREMODE SEQUENCEMODE COMMENTMODE
%s PRESEQUENCE PRESTRUCTURE PREENSEMBLE PROPERTYLIST ENSEMBLEMODE
%%

<*>[ \t\r]     LOCATION_STEP (*yylloc);
<*>\\[ ]*\s*\r*\n        {
    LOCATION_LINES (*yylloc, yyleng); 
    LOCATION_STEP (*yylloc); 
}
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
<ENSEMBLEMODE>\{ {return LBRACE;}
<ENSEMBLEMODE>\} {return RBRACE;}
<PRESTRUCTURE>= {BEGIN(STRUCTUREMODE); return EQUALS;}
<PREENSEMBLE>= {BEGIN(ENSEMBLEMODE); return EQUALS;}
<STRUCTUREMODE,ENSEMBLEMODE>\+ {return PLUS;}
<ENSEMBLEMODE>- {return MINUS;}
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
<*>, {return COMMA;};
\* return TOK_STAR;
<INITIAL,PROPERTYLIST>\.  {BEGIN(PROPERTYLIST); return PERIOD;}
<STRUCTUREMODE>\.  {return PERIOD;}
[a-zA-Z_][a-zA-Z0-9_-]* {
    int lookupval = lookup_reserved(yytext);
    if (lookupval == TOK_STRUCTURE) {
        BEGIN(PRESTRUCTURE);
    } else if (lookupval == TOK_OFFTARGETS) {
        BEGIN(PREENSEMBLE);
    } else if (lookupval == TOK_DOMAIN || lookupval == TOK_LIBRARY
            || lookupval == TOK_PREVENT || lookupval == TOK_SOURCE
            || lookupval == TOK_SIMILARITY) {
        if (!((lookupval == TOK_SOURCE || lookupval == TOK_SIMILARITY) && YY_START == PROPERTYLIST)) {
            BEGIN(PRESEQUENCE);
        }
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
