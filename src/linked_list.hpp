/**
 *  File: linkedlist.h
 *  Author: Allen Sanford (ras9841@rit.edu)
 *  Description:
 *      Outlines the LinkedList class generically.
 *  Date Creaded: 08/03/16
 */

#ifndef __LINKEDLIST__
#define __LINKEDLIST__

#include <iostream>

template <class T>
class LinkedList {
/**
 *  Class outline for a Node object.
 */
        struct Node {
                T data;
                Node *next;
                Node *prev;
        };
        Node *head;
        size_t size;
    public:
        LinkedList<T> ();
        ~LinkedList<T> ();
        T getHead ();
        T popHead ();
        void prepend (T val);
        void append(T val);
        void remove(T val);
        void printList();
};

template <class T>
LinkedList<T>::LinkedList()
{
/**
 * Consturctor for empty LinkedList.
 *
 * @param val - generic value to store.
 */
    this->size = 0;
    this->head = nullptr;
}

template <class T>
LinkedList<T>::~LinkedList()
{
/**
 * Desturctor for LinkedList.
 *
 * @param val - generic value to store.
 */
    Node *n = this->head;
    
    while (n != nullptr)
    {
        if (n->next == nullptr)
        {
            delete n;
            n = nullptr;
        }
        else 
        {
            n = n->next;
            delete n->prev;
        }
    }
    this->size = 0;
}

template <class T>
T LinkedList<T>::getHead()
/**
 * Getter method for the head of a linked list. 
 *
 * @return the list's head
 */
{
    if (this->size) return this->head->data;
    else
    {
        Node n = Node();
        return n.data;
    }
}

template <class T>
T LinkedList<T>::popHead()
/**
 * Getter method for the head of a linked list. The head is removed. 
 *
 * @return the list's head
 */
{
    if (this->size) {
        Node *n = this->head;
        this->head = n->next;
        return n->data;
    }
    else
    {
        return nullptr;
    }
}

template <class T>
void LinkedList<T>::printList()
/**
 * Prints the linked list to standard output starting from
 * the current node (this).
 */
{
    Node *node = this->head;
    std::cout << "[";
    
    if (this->size)
    {
        std::cout << node->data;
        while (node->next != nullptr)
        {
            node = node->next;
            std::cout << ", " << node->data;
        }
    }
    std::cout << "]\n";
}

template <class T>
void LinkedList<T>::prepend(T val)
/**
 * Creates a new node at the head of the linked list
 * with data val.
 *
 * @param val - value to prepend to the list.
 */
{
    Node *n = new Node();
    n->data = val;
    n->next = this->head;
    n->prev = nullptr;

    this->head = n;
}

template <class T>
void LinkedList<T>::append(T val)
/**
 * Creates a new node at the tail of the linked list
 * with data val.
 *
 * @param val - value to append to the list.
 */
{
    Node *n = new Node();
    n->data = val;
    n->next = nullptr;

    if (!this->size)
    {
        n->prev = nullptr;
        this->head = n;
    }
    else
    {
        Node *node = this->head;
        while (node->next != nullptr)
        {
            node = node->next;
        }
        n->prev = node;
        node->next = n;
    }
    this->size++;
}


template <class T>
void LinkedList<T>::remove(T val)
/**
 * Removes the node containig the specified value from the list.
 *
 * @param val - value to be removed from the list.
 */
{
    Node *node = this->head;
    while (node != nullptr)
    {
        if (node->data == val)
        {
            if (node == this->head)
            {
                this->head = node->next;
            }
            if (node->prev != nullptr) //update previous' next node
            {
               node->prev->next = node->next;
            }
            if (node->next != nullptr) //update next's previous node
            {
               node->next->prev = node->prev;
            }
            delete node;
            this->size--;
        }
        node = node->next;
    }
}
#endif
