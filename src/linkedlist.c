#include <stdio.h>
#include <stdlib.h>
#include "linkedlist.h"

/*
	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	LINKEDLIST

	Authors:
		Ramsey Karim
		Harrison Agrusa
		Julian Marohnic

	Description:
		Linked list implementation that will be useful for quick
			sorting of particles using flexible lists.
		This linked list is for integers.
		If we need a list for another type, we'll modify.
		We did not need a list for another type, so this is great.

	Usage:
		Used as a library, called elsewhere.

	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/

LLn * ll_append(LL * list, int newValue) {
	// Return the pointer to the new node at tail
	if (!list) return NULL;
	LLn * newNode = (LLn*) malloc(sizeof(LLn));
	if (!newNode) return NULL;
	newNode->value = newValue;
	newNode->next = NULL;
	if (!list->head || !list->tail) {
		list->head = newNode;
		list->tail = newNode;
	} else {
		list->tail->next = newNode;
		list->tail = newNode;
	}
	return newNode;
}

LLn * ll_append_node(LL * list, LLn * node) {
	if (!list || !node) return NULL;
	node->next = NULL;
	if (!list->head || !list->tail) {
		list->head = node;
		list->tail = node;
	} else {
		list->tail->next = node;
		list->tail = node;
	}
	return node;
}

LLn * ll_pop(LL * list) {
	// Pop off the FIRST element (head) and return it
	// This does not delete the node, but it removes it from the list
	LLn * returnNode, * newHead;
	returnNode = list->head;
	if (returnNode == NULL) {
		return NULL;
	}
	newHead = returnNode->next;
	list->head = newHead;
	if (list->tail == returnNode) {
		list->tail = NULL;
	}
	return returnNode;
}

LL * ll_create_list() {
	LL * newList = (LL*) malloc(sizeof(LL));
	newList->head = NULL;
	newList->tail = NULL;
	return newList;
}

void ll_delete_list(LL * list) {
	if (list->head) {
		ll_delete_node(list->head);
	}
	free(list);
}

void ll_delete_node(LLn * node) {
	if (node->next) {
		ll_delete_node(node->next);
	}
	free(node);
}

LLn * ll_next(LLn * node) {
	return node->next;
}

int ll_length(LL * list) {
	if (!list->head) {
		return 0;
	} else {
		return ll_length_helper(list->head);
	}
}

int ll_length_helper(LLn * node) {
	if (!node->next) {
		return 1;
	} else {
		return ll_length_helper(node->next) + 1;
	}
}

int ll_quick_length(LL * list) {
	// If the list is longer than 1, return -1
	if (!list->head || !list->tail) {
		return 0;
	} else if (list->head == list->tail) {
		return 1;
	} else return -1;
}

LL * int_array_to_ll(int n) {
	/* Useful if you have n integers (0 to n-1)
		and you want a linked list with each
		integer represented
		Hmmmm wonder how we could use this ;)
	*/
	if (n < 1) {
		printf("Bad input (n = %d) to INT_TO_ARRAY_LL\n", n);
		return NULL;
	}
	LL * list = ll_create_list();
	for (int i = 0; i < n; i++) {
		ll_append(list, i);
	}
	return list;
}
