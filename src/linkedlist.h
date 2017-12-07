#ifndef LINKEDLIST_H
#define LINKEDLIST_H

typedef struct LinkedList LL;
typedef struct LinkedListNode LLn;

// This list is for integers.
// If we need a list for another type, we'll modify.
// This will be useful for quick sorting of particles
struct LinkedList
{
	LLn * head;
	LLn * tail;
};

struct LinkedListNode
{
	int value;
	LLn * next;
};

// Append to end
LLn * ll_append(LL * list, int newValue);

LLn * ll_append_node(LL * list, LLn * node);

LLn * ll_pop(LL * list);

// Initializes with head and tail pointing to NULL
LL * ll_create_list();

void ll_delete_list(LL * list);

void ll_delete_node(LLn * node);

LLn * ll_next(LLn * node);

int ll_length(LL * list);

int ll_length_helper(LLn * node);

int ll_quick_length(LL * list);

LL * int_array_to_ll(int n);

#endif