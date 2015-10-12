/*!
  \file   parms_table_impl.h
  \brief  parms_Table is used for parms_Map object to insert pairs
          (global_index, local_index) into the table. 

  \author zzli
  \date   2006-05-05
*/

#ifndef _PARMS_TABLE_IMPL_H_
#define _PARMS_TABLE_IMPL_H_

#include "parms_mem.h"
#include "parms_table.h"

/*! typedef Slot
*/
/*! \struct Slot
  \brief structure Slot
 */
typedef struct Slot {
  struct Slot *link;		//!< pointer to a next entry in this slot.
  int key;			//!< key in a pair (key,value)
  int value;			//!< value in a pair (key, value)
} Slot;
  
/*! \struct parms_Table_
  \brief parms_Table_ structure.
 */
struct parms_Table_ {
  Slot     **Slots;    //!< array of a pointer array to struct Slot 
  HashFcn  hf;	       //!< function pointer to a hash function.
  int      space;      //!< size of array Slots.
  int      size;       //!< number of pairs in the table.
};

#endif  
