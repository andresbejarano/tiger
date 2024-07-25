# -*- coding: utf-8 -*-

class Queue:
    
    def __init__(self):
        from collections import deque
        self.__d = deque()
        self.__count = 0
    
    
    def dequeue(self):
        assert self.__count > 0, "Queue is empty"
        self.__count -= 1
        return self.__d.popleft()
    
    
    def enqueue(self, item):
        self.__d.append(item)
        self.__count += 1
    
    
    def isempty(self) -> bool:
        return self.__count == 0
    
    
    def size(self) -> int:
        return self.__count
