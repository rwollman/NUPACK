#include "design_utils.h"

void testPatternPrevention(void);
void createPPCase(char *seq, char *skel, int violates, int actualViolations);
void createSkeletonCase(char *seq, char *skel, int violates);
void testConstraints(void);
void testParensConverter(void);
void createParensConverterCase(char *parensWithPlus, char *parensWithoutPlus, int *breaks, int numBreaks);
void testNSTheta(void);
DBL_TYPE createNSThetaCase(char *seq, int *pairing, int *fakeBases, int numStrands, int *breaks, DBL_TYPE *pm);
void unitTester(void);
