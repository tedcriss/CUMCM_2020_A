# -*- coding: utf-8 -*-
"""
@author:xuyuntao
@time:2020/9/4 21:04
@email:xuyuntao@189.cn
"""
# -- coding:utf-8 --
import numpy
import csv
# import math
import sympy
import matplotlib.pyplot as plt

disInt=0.01#s
v=7/6#cm/s
k=3
tau=108
ovenTemp=[25,175,195,235,255,25]
airIntAdj=[1.5, 0.2, 0.4, 0.2, 2.3]#cm
ovenLength=30.5#cm
airLength=5#cm
dis2Oven=25#cm

p0=[0,ovenTemp[0]]
p1=[dis2Oven + airIntAdj[0] * k, ovenTemp[1]]
p2=[p1[0] - airIntAdj[0] * k + (5 * ovenLength + (5 - 1) * airLength) - airIntAdj[1] / 2 * k, ovenTemp[1]]
p3=[p2[0] + airLength + airIntAdj[1] * k, ovenTemp[2]]
p4=[p3[0] - airIntAdj[1]/2 * k + ovenLength - airIntAdj[2] / 2 * k, ovenTemp[2]]
p5=[p4[0] + airLength + airIntAdj[2] * k, ovenTemp[3]]
p6=[p5[0] - airIntAdj[2]/2 * k + ovenLength - airIntAdj[3] / 2 * k, ovenTemp[3]]
p7=[p6[0] + airLength + airIntAdj[3] * k, ovenTemp[4]]
p8=[p7[0] - airIntAdj[3]/2 * k + ovenLength * 2 + airLength - airIntAdj[4] / 2 * k, ovenTemp[4]]
p9=[p8[0] + airLength + airIntAdj[4] * k, ovenTemp[5]]
p10=[p9[0] - airIntAdj[4]/2 * k + ovenLength * 2 + airLength+dis2Oven, ovenTemp[5]]
changeIndex=[]#切换点索引

temp=numpy.zeros([int(p10[0]/disInt)],dtype=numpy.float32)#oC
dis=temp.copy()#cm
temp0=25*temp.copy()#oC
preTemp0=25

def calDis(startDis,endDis):
    """输入距离，返回索引、距离矩阵"""
    startIndex=int(startDis/disInt)
    endIndex=int(endDis/disInt)
    subDis=numpy.arange(startDis,endDis,disInt)[:endIndex-startIndex]
    return startIndex,endIndex,subDis

firstFlag=True
for _ in range(11-1):
    startIndex,endIndex,subDis=calDis(eval("p"+str(_))[0],eval("p"+str(_+1))[0])
    if firstFlag:
        changeIndex.append(startIndex)
        firstFlag=False
    changeIndex.append(endIndex)
    # print(startIndex,endIndex)
    subTemp=(subDis-eval("p"+str(_+1))[0])/(eval("p"+str(_))[0]-eval("p"+str(_+1))[0]) *\
            (eval("p"+str(_))[1]-eval("p"+str(_+1))[1]) + eval("p"+str(_+1))[1]
    # print("subTemp",subTemp)
    temp[startIndex:endIndex]+=subTemp
    dis[startIndex:endIndex]+=subDis

# plt.plot(dis,temp)
time=dis/v#s

tempT=numpy.zeros(temp.shape,dtype=numpy.float32)
for _ in range(11-1):
    temp0[changeIndex[_]:changeIndex[_ + 1]]=preTemp0
    subTempT = temp0[changeIndex[_]:changeIndex[_+1]] + (temp[changeIndex[_]:changeIndex[_+1]] - temp0[changeIndex[_]:changeIndex[_+1]]) *\
               (1 - numpy.exp(-1 * time[changeIndex[_]:changeIndex[_+1]] / tau))
    preTemp0=subTempT[-1]
    tempT[changeIndex[_]:changeIndex[_ + 1]]=subTempT


# plt.subplot(2,1,1)
plt.plot(dis,temp)
plt.xlabel("distance (cm)")
plt.ylabel("temperature (°C)")
# plt.subplot(2,1,2)
# plt.plot(dis/v,tempT)
# plt.xlabel("time (s)")
# plt.ylabel("temperature (°C)")
plt.title("Reflow Oven Temperature Curve")

# todo:按0.5s间隔取出时间列表近似值
