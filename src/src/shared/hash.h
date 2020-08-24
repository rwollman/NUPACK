/* Copyright 2006 David Crawshaw, released under the new BSD license.
 * Version 2, from http://www.zentus.com/c/hash.html 

 * Modified by C. Steenberg <conrad.steenberg@caltech.edu>
   02/12/2008 Make it possible to use keys with \0 chars
 */

#ifndef NUPACK_SHARED_HASH_H
#define NUPACK_SHARED_HASH_H
 
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */




/* Opaque structure used to represent hashtable. */
typedef struct hash hash;

/* Create new hashtable. */
hash * hash_new(unsigned int size);

/* Free hashtable. */
void hash_destroy(hash *h);

/* Add key/value pair. Returns non-zero value on error (eg out-of-memory). */
int hash_add(hash *h, const char *key, const unsigned int len, void *value);

/* Return value matching given key. */
void * hash_get(hash *h, const char *key, const unsigned int len);

/* Remove key from table, returning value. */
void * hash_remove(hash *h, const char *key, const unsigned int len);

/* Returns total number of keys in the hashtable. */
unsigned int hash_size(hash *h);

/* Calculate a hash from a given array of characters that need not be
  zero-terminated
*/
unsigned int vechash(const char *vec, const unsigned int len);

#ifdef __cplusplus
}
#endif 

#endif

