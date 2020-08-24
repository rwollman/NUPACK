/** \file GetEnergy.h
 *  GetEnergy.h by Robert Dirks.  
 *  This program determines the energies of
 *  substructures, mirroring the loops of Fold.out
 *  and includes the function for determining energies.  03/15/2001  
 */

#ifndef __GETENERGY_H__
#define __GETENERGY_H__

#include "ene.h"

#ifdef __cplusplus
extern "C" {
#endif

//These functions are used to calculate the energy of a structure (linear time algorithm)

//compute the energy of the structure described in *thefold
DBL_TYPE GetEnergy( fold *thefold);

//Calculate the energy of an exterior loop substructure
DBL_TYPE EnergyF( int start, int stop, fold *thefold);
//Energy of Paired substructure
DBL_TYPE EnergyFb( int start, int stop, fold *thefold);
//Energy of pseudoknotted substructure
DBL_TYPE EnergyPk( int i, int j, fold *thefold);
//Energy of gap spanning region
DBL_TYPE EnergyFg( int i, int d, int e, int j, fold *thefold);
//Energy of pseudoknot interior
DBL_TYPE EnergyFz( int start, int stop, fold *thefold);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
