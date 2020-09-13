# -*- coding: utf-8 -*-
"""
@author:xuyuntao
@time:2020/9/10 18:58
@email:xuyuntao@189.cn
"""
import csv
import os
import matplotlib.pyplot as plt
import numpy

data=[]
firstFlag=True
with open(r"../data.csv", "r") as f:
    csvRea=csv.reader(f)
    for _ in csvRea:
        if firstFlag:
            data.append(_)
            firstFlag=False
        else:
            subData=[]
            for __ in _:
               subData.append(float(__))
            data.append(subData)
dataNp=numpy.array(data[1:],dtype=numpy.float32)

plt.plot(dataNp[:,0],dataNp[:,1])
plt.show()