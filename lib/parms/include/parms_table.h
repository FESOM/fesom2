/**
 * @file   parms_table.h
 * @author Zhongze Li
 * @date   Tue Oct 17 11:59:56 2006
 * 
 * @brief  The hash table functions. 
 * 
 * 
 */

#ifndef _PARMS_TABLE_HEADER_H
#define _PARMS_TABLE_HEADER_H

#include "parms_sys.h"

PARMS_CXX_BEGIN

typedef unsigned int (*HashFcn)(int key);
typedef struct parms_Table_ *parms_Table;

/** 
 * Create a hash table.
 *
 * Create a hash table with hf as the hash function.
 * 
 * @param newT  A pointer to the hash table created.              
 * @param hf    The hash function. If hf is null, the default hash
 *              function is used.     
 * @param size  The number of entries stored in the table.
 * 
 * @return 0 on success.
 */
extern int parms_TableCreate(parms_Table *newT, HashFcn hf, int size);

/** 
 * Free the memory for the table object pointed by self.
 * 
 * @param self A pointer to the parms_Table object.  
 * 
 * @return 0 on success.
 */
extern int parms_TableFree(parms_Table *self);

/** 
 * Get the corresponding value for a given key.
 *
 * If return NULL, then the entry with key is not in the table,
 * otherwise return a pointer to the value.
 * 
 * @param self  A parms_Table object. 
 * @param key   The key value.     
 * 
 * @return A pointer to the value.
 */
extern void *parms_TableGet(parms_Table self,  int key);

/** 
 * Put the pair (key, value) into the table.
 * 
 * @param self A parms_Table object.        
 * @param key  The key of the pair.         
 * @param val  The value of the pair.       
 * 
 * @return 0 on success.
 */
extern int parms_TablePut(parms_Table self,  int key, int val);

/** 
 * Put the slot corresponding to the key from the table.
 * 
 * @param table A parms_Table object.        
 * @param key 	The key of the pair to be removed.         
 * 
 * @return 0 on success.
 */
extern int parms_TableRemoveFromLast(parms_Table self, int key);

/** 
 * Swap value of key1 with that of key2.
 * 
 * @param table A parms_Table object.        
 * @param key1, key2 	keys to be swapped.         
 * 
 * @return 0 on success.
 */
extern int parms_TableSwap(parms_Table self, int key1, int key2);

/** 
 * Get the number of entries stored in the table.
 * 
 * @param self A hash table object.
 * 
 * @return The number of entries stored in the table.
 */
extern int parms_TableGetSize(parms_Table self);

PARMS_CXX_END

#endif  
