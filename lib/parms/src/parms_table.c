/*--------------------------------------------------------------------
  parms_TableCreate   : create a hash table.
  parms_TableFree     : free the memory for the table.
  parms_TableGet      : get the corresponding value for a given key.
  parms_TableGetSize  : get the total number of entries in the table.
  parms_TablePut      : put the pair (key, value) into the table.

  A code fragment for using table functions:

  parms_Table table;
  int key, value, *p;

  // create a hash table using the default hash function
  parms_TableCreate(&table, NULL, tsize);
  // put pair(key, value) into the table
  parms_TablePut(table, key, value);
  // get the corresponding value for a given key.
  p = parms_TableGet(table, key);
  if (p != NULL) {
    value = *p;
  }
  // free the memory for the table object.
  parms_TableFree(&table);
  
  $Id: parms_table.c,v 1.1.1.1 2006-11-27 22:28:01 zzli Exp $
  ------------------------------------------------------------------*/

#include "include/parms_table_impl.h"

static unsigned  hashfunc(int key)
{
  unsigned k;

  k = key;
  return k;
}

/** 
 * Create a hash table.
 *
 * Create a hash table with hf as the hash function.
 * 
 * @param newT  A pointer to the hash table created.              
 * @param hf 	The hash function. If hf is null, the default hash
 *              function is used.     
 * @param size  The number of entries stored in the table.
 * 
 * @return 0 on success.
 */
int parms_TableCreate(parms_Table *newT, HashFcn hf, int tsize)
{
  int i, lwr;
  parms_Table table;
  static int primes[] = {389, 389, 769, 1543, 3079, 6151, 12289,
			 24593, 49157, 98317, 196613, 393241, 786433,
			 1572869, 3145739, 6291469, 12582917,
			 25165843, 50331653, 100663319, 201326611,
			 402653189, 805306457, 1610612741}; 

  /* find the maximum prime less than tsize */
  for (i = 1; lwr = (primes[i]*2)/3, lwr < tsize; i++);

  /* malloc the space for the table */
  table = PARMS_ALLOC(sizeof(*table) +
		      sizeof(table->Slots[0])*primes[i-1]); 
  table->space = primes[i-1];
  table->Slots = (struct Slot **)(table+1);
  for (i = 0; i < table->space; i++) {
    table->Slots[i] = NULL;
  }

  table->hf = hf ? hf : hashfunc;
  table->size = 0;
  *newT = table;
  return 0;
}

/** 
 * Get the corresponding value for a given key.
 *
 * If return NULL, then the entry with key is not in the table,
 * otherwise return a pointer to the value.
 * 
 * @param table A parms_Table object. 
 * @param key   The key value.     
 * 
 * @return A pointer to the value.
 */
void *parms_TableGet(parms_Table table, int key)
{
  int index;
  struct Slot *p;
  
  index = table->hf(key) % table->space;
  for (p = table->Slots[index]; p; p = p->link) {
    if (p->key == key) {
      break;
    }
  }
  return p ? &p->value : NULL;
}

/** 
 * Put the pair (key, value) into the table.
 * 
 * @param table A parms_Table object.        
 * @param key 	The key of the pair.         
 * @param val 	The value of the pair.       
 * 
 * @return 0 on success.
 */
int parms_TablePut(parms_Table self, int key, int value)
{
  int index;
  struct Slot *p;

  index = self->hf(key) % self->space;
  for (p = self->Slots[index]; p; p = p->link) {
    if (p->key == key) {
      break;
    }
  }
  if (p == NULL) {		/* new entry */
    PARMS_NEW(p);
    p->key = key;
    p->value = value;
    p->link = self->Slots[index];
    self->Slots[index] = p;
    self->size++;
  }
  else {
    p->key = key;
    p->value = value;
  }
  return 0;
}

/** 
 * Remove the slot corresponding to the key from the table.
 * NOTE: This assumes that the key/value pairs are stored in 
 * a way that needs to be preserved after the removal. For instance,
 * the pARMS matrix object stores the keys as global indices, and 
 * values as local indices. Thus each of these are unique, and the 
 * local indices are increasing with respect to the size of the table
 * or count of entries in the table. Thus removing an entry requires
 * that this relation is preserved. So this removes the slot by first 
 * swapping values with the last entered entry in the table, before 
 * actually removing the slot. 
 * 
 * @param table A parms_Table object.        
 * @param key 	The key of the pair to be removed.         
 * 
 * @return 0 on success.
 */
int parms_TableRemoveFromLast(parms_Table self, int key)
{
  int index, i, size;
  struct Slot *p, *p1, *p2, *p3, *plast;
  
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  index = self->hf(key) % self->space;
  /* p1 points to the physical entry just before slot p */
  p1 = self->Slots[index];
  for (p = self->Slots[index]; p; p = p->link) {    
    if (p->key == key) {
      break;
    }
    p1 = p;
  }
  if(p != NULL)
  {  
    /* Search for last entered entry to swap values with p*/
    /* This entry is most likely at the leftmost entry */
    size = self->size;
    plast = p;/*initialization*/
    for (i = 0; i < self->space; i++) {
      plast = self->Slots[i];
      if(plast && plast->value == size-1) break;
    }
    /* Now swap values of p and plast to maintain table
     * consistency 
     */

    parms_TableSwap(self,p->key,plast->key);

   /* Now break links by letting p1 link to the next physical entry after p */
    p2 = p->link;
    if(p1 == p)/* p is starting slot */
       self->Slots[index] = p2;
    else
     p1->link = p2;      
    /* Now free p to actually remove entry */
     PARMS_FREE(p);
    
  /* now decrement table size */
    self->size--;

    /* Now update last entered entry of table row
     * containing plast, by swapping with other table
     * row entries, till appropriate position is reached.
     * This ensures that entries to the left of a table row are 
     * always of larger value (or are the latest entered entries)
     * of the table row.
     */
    p3 = plast->link;
    for(p3 = plast; p3->link; p3 = p3->link)
    {
      if((p3->link)->value <= p3->value) break;
      /*swap*/
      parms_TableSwap(self,p3->key,(p3->link)->key);
    }

  }

  return 0;
}

/** 
 * Swap value of key1 with that of key2.
 * 
 * @param table A parms_Table object.        
 * @param key1, key2 	keys to be swapped.         
 * 
 * @return 0 on success.
 */
int parms_TableSwap(parms_Table self, int key1, int key2)
{
  int *value1, *value2, temp1, temp2;

  value1 = parms_TableGet(self,key1);
  value2 = parms_TableGet(self,key2);

  if(value1 == NULL || value2 == NULL)
  {
    if(value1==NULL)printf("Cannot swap with a null entry: parms_TableSwap error...key: %d\n", key1);
    if(value2==NULL)printf("Cannot swap with a null entry: parms_TableSwap error...key: %d\n", key2);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  else
  {
     temp1 = *value1;
     temp2 = *value2;
     parms_TablePut(self, key1, temp2);
     parms_TablePut(self, key2, temp1);
  }
  return 0;
}

/** 
 * Free the memory for the table object pointed by self.
 * 
 * @param table A pointer to the parms_Table object.  
 * 
 * @return 0 on success.
 */
int parms_TableFree(parms_Table *self)
{
  struct Slot *p1, *p2;
  int i;

  for (i = 0; i < (*self)->space; i++) {
    for (p1 = (*self)->Slots[i]; p1; p1 = p2) {
      p2 = p1->link;
      PARMS_FREE(p1);
    }
  }
  PARMS_FREE(*self);
  return 0;
}

/** 
 * Get the number of entries stored in the table.
 * 
 * @param table A hash table object.
 * 
 * @return The number of entries stored in the table.
 */
int parms_TableGetSize(parms_Table self)
{
  return self->size;
}
