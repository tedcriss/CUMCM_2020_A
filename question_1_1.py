# -*- coding: utf-8 -*-
"""
@author:xuyuntao
@time:2020/9/13 11:23
@email:xuyuntao@189.cn
"""
import numpy
import csv
import sympy
from judgeTempFunction import judgeTempgA
from calTt import TtF
from calTt import TatF
import matplotlib.pyplot as plt
sympy.init_printing()

fullLength = 25 + 25 + 11 * 30.5 + 10 * 5

vNp = 70 / 60
tNp = numpy.arange(0,fullLength / vNp, 0.5)
tauNp = 45
t1Np = 175
t2Np = 195
t3Np = 235
t4Np = 255

ovenLength=30.5#cm    # 回焊炉温区长度
airLength=5#cm                 # 回焊炉各温区间隔
dis2Oven=25#cm                 # 回焊炉路前后区域长度
coolLength=ovenLength * 2               # 实际冷却区缩减

data=[]
firstFlag=True
with open(r"data.csv", "r") as f:
    csvRea=csv.reader(f)
    for _ in csvRea:
        data.append(_)
        # if firstFlag:
        #     firstFlag=False
        # else:
        #     subData=[]
        #     for __ in _:
        #         subData.append(float(__))
        #     data.append(subData)
givenDataNp=numpy.array(data[1:],dtype=numpy.float32)
# dataNpFull[37:]=givenDataNp[:,1]
givenDataNp_SameLen=numpy.ones([tNp.size,2],dtype=numpy.float32)*30
givenDataNp_SameLen[:,0]=tNp
startTimeIndex=numpy.where(tNp==givenDataNp[0,0])[0][0]

givenDataNp_SameLen[startTimeIndex:,:]=givenDataNp

deltaList=[]     # 误差矩阵
minChar=None    # 最小误差值
minTtNp=None    # 最小误差对应T(t)
minTau=None     # 最小误差对应tau值
useTau=471/10     # 从另一处得到的合适的tau
useChar=None    # 使用tau的误差
useTtNp=None    # 使用tau求出的T(t)
for _ in range(250,1100):    # 遍历所有tau
    tauNp=_/10    # 修改tau为分式，方便计算
    # 提取条件：Tat.args[方程号][1：条件].args[第几个条件]._args[0：t，1：数值]
    # 提取公式：Tat.args[方程号][0：公式].as_poly().coeffs() -> 其中，0为a，1为b
    print("当前tau=",tauNp)

    TtNp=TtF(tNp,vNp,tauNp,t1Np,t2Np,t3Np,t4Np)


    char=numpy.sum((TtNp-givenDataNp_SameLen[:,1])**2)/TtNp.shape[0]   # 计算平方差
    if (minChar==None) or (char<=minChar):
        minChar=char
        minTtNp=TtNp
        minTau=tauNp
    if (tauNp==useTau):
        useChar=char
        useTtNp=TtNp
    deltaList.append([tauNp, char])
    print("char=",char)
    print("--------------------------------")

print("deltaList", deltaList)

plt.figure(1)
plt.plot(tNp,useTtNp)
plt.plot(tNp,givenDataNp_SameLen[:,1])
plt.text(tNp[numpy.where(useTtNp==useTtNp.max())[0][0]],useTtNp.max(),"maximum="+str(round(useTtNp.max(),2)))
plt.scatter(tNp[numpy.where(useTtNp==useTtNp.max())[0][0]],useTtNp.max())
plt.legend(["$T(t)$","Given Data"])
# plt.plot(tNp,TatNp)
# plt.plot(tNp,dataNpFull-minTtNp)
plt.xlabel("time (s)")
plt.ylabel("temperature (°C)")
plt.title("$T(t)$ and Original Temperature, $\\tau={}$".format(useTau))
# plt.savefig("figure1.png",dpi=400)

plt.figure(2)
charListNp=numpy.array(deltaList, dtype=numpy.float32)
plt.plot(charListNp[:,0],charListNp[:,1])
plt.text(minTau,minChar,"min{$\\tau$}="+str(minTau))
plt.text(useTau,useChar,"use $\\tau$="+str(useTau))
plt.scatter(minTau,minChar)
plt.scatter(useTau,useChar)
plt.xlabel("$\\tau$")
plt.ylabel("$\\delta$")
plt.title("$\\delta$")
# plt.savefig("figure2.png",dpi=400)

TatNp=TatF(tNp,vNp,tauNp,t1Np,t2Np,t3Np,t4Np)
plt.figure(3)
plt.plot(tNp,TatNp)
plt.grid(True)     # 显示网格线
plt.xlabel("time (s)")
plt.ylabel("temperature (°C)")
plt.title("$T_a(t)$")
# plt.savefig("figure3.png",dpi=400)