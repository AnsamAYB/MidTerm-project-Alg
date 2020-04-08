// CPCS 324 Algorithms & Data Structures 2
// Outline - Priority queue data structure
// 2019, Dr. Muhammad Al-Hashimi


// -----------------------------------------------------------------------
// Basic design decisions and implementation planning (objects & interfaces)

// initial requirements: to quickly support Dijkstra and second Prim's algorithms, 
// implement minimum priority functionality

// design decisions:
// Reuse the 324 linked list implementation
// how will the PQ ops be implemented?
// <student fill here>

// code plan: start to write your API specifications (JSDOC comments) and think about 
// design consequences (design impact)

// Impact analysis:
// <student fill here>



// -----------------------------------------------------------------------
// Priority queue object constructor (document using JSDOC comments)



function PQueue()
{
	 /** Head pointer of a queue  */
	this.pq = new List();          
	
    // specify (design) methods
    
	/**
      Return true if queue is empty
      @method
    */
	this.isEmpty = isEmptyImpl;                   

	/**
      Remove and return item with minimum priority
      @method
     */
	this.deleteMin = deleteMinImpl;                

	/**
      Insert or update an item with priority
      @method
     */
	this.insert = insertImpl;                   
	
}

// -----------------------------------------------------------------------
// Priority queue node constructor (document using JSDOC comments)


function PQNode(item, key)
{
	/**
      The value of Data item 
    */
	this.item = item;

    /**
      The value of priority 
    */
	this.prior = key;
	
	// specify (design) methods
	
}

// -----------------------------------------------------------------------
// functions used by PQueue() object methods
// specify interface information (JSDOC comments)
// function names should not clash with linklist.js and queue.js
// ....



function isEmptyImpl()
{
    // check if PQ is empty or not
	return (this.pq.isEmpty());
	
}




function deleteMinImpl()
{
       // pointer ponits to the first node in the PQ 
       var pointer = this.pq.first;
       // the node of item with minimum value of key is the first node of PQ (intial value)
       var minimum = this.pq.first;
       // initialize the previous node of node of item with minimum value with null value
       var previous = null ;
       // search for the node of item with the minimum value of key 
       while(pointer != null){
           // check if the key of the current item is smaller than the key of the minimum item
         if((pointer.next != null) && (pointer.next.item.prior < minimum.item.prior)){
           // update the value of node of item with the minimum value of key
           minimum = pointer.next ;
           // update the value of node that has previous item of the minimum item 
           previous = pointer ;
         }
         // to move to the next node in PQ
         pointer = pointer.next;
       }
       // check if node of item with minimum value of key is the first one in PQ
       if(minimum == this.pq.first){ 
          // return the data item of item with the minimum value of key and delete it
          return this.pq.deleteFirst().item;
       }else{ // if node of item with minimum value of key is not the first one in PQ 
          // previous node points to the node next to the node with minimum item 
          previous.next = minimum.next;
          // return the data item of item with the minimum value of key
          return minimum.item.item;
       }
}



function insertImpl(item, key)
{
    // creat new PQNode to store the item and key 
    var item = new PQNode(item,key);
    // pointer ponits to the first item in the PQ
    var pointer = this.pq.first;
    // boolean variable : true in insert case and false in update case
    var toInsert = true;
    // check if the PQ is empty so it must be the insert case 
    if (this.isEmpty()) {
        // insert the item in the first of the PQ
        this.pq.insert(item);
        // change toInsert variable to false to not insert again the first position
        toInsert = false;

    } else { // if PQ is  not empty so it might be either the insert case or update case
        // search for the received item if it exists then it must be the update case 
        while (pointer != null) { 
            // if the data item of the received item equal to some data item for an item in PQ
            if (item.item == pointer.item.item) { 
                // update the priority value 
                pointer.item.prior = key;
                //change toInsert variable to false to not insert the updated item
                toInsert = false;
                // break the loop
                break; 
            }
            // to go to the next item in the PQ
            pointer = pointer.next;
        }

    }
    // insert case
    if (toInsert) {  
        // insert item in PQ
        this.pq.insert(item);
    }
    
}
