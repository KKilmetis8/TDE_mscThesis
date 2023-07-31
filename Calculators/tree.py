#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 12:02:08 2023

@author: konstantinos
"""
import numpy as np

class Node:
    def __init__(self, val):
        self.val = val
        self.left = None
        self.right = None
        
def insert(node, val):
    ''' Insert a Node '''
    # Return a new node if the tree is empty
    if node is None:
        return Node(val)

    # Traverse to the right place and insert the node
    if val < node.val:
        node.left = insert(node.left, val)
    else:
        node.right = insert(node.right, val)

    return node

def build_tree(array):
    root = None
    for element in array:
        root = insert(root, element)
    return root        

# def search_tree(element, tree):
#     tol = 1e-4
#     diff = np.abs(tree.val - element)
#     if diff < tol:
#         return tree.val
    
#%%
if __name__ == '__main__':
    import time
    
    array = np.linspace(0, 1, num = 72091)
    barray = np.random.random(1_000_000)
    holder = []
#%% Brute Force   
    start_brute_time = time.time()
    for b in barray:
        diffs = np.abs(b - array)
        idx = np.argmin(diffs)
        holder.append(array[idx])
    end_brute_time = time.time()
    print('Brute Force: ', end_brute_time - start_brute_time)
    #%% Tree
    start_tree_time = time.time()
    tree = build_tree(array)
    for b in barray:
        closest = search_bst(b)
        idx = np.where(array == closest)
        holder.append(array[idx])
    end_tree_time = time.time()
    print('Tree Force: ', end_tree_time - start_tree_time)