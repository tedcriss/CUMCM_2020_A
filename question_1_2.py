# -*- coding: utf-8 -*-
"""
@author:xuyuntao
@time:2020/9/13 11:55
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

vNp = 78 / 60
tNp = numpy.arange(0,fullLength / vNp, 0.5)
tauNp = 45
t1Np = 173
t2Np = 198
t3Np = 230
t4Np = 257

ovenLength=30.5#cm    # 回焊炉温区长度
airLength=5#cm                 # 回焊炉各温区间隔
dis2Oven=25#cm                 # 回焊炉路前后区域长度
coolLength=ovenLength * 2               # 实际冷却区缩减
timeIntv=0.5                            # 计算时间间隔

validList=[]
minTau=None
# validTtList=[]
minTauTtNp=None
allTau=[]
maxTt=[]
for _ in range(40,1000):
    tauNp=_/10
    allTau.append(tauNp)
    # 提取条件：Tat.args[方程号][1：条件].args[第几个条件]._args[0：t，1：数值]
    # 提取公式：Tat.args[方程号][0：公式].as_poly().coeffs() -> 其中，0为a，1为b
    print("当前tau=",tauNp)

    # TtNpF=sympy.lambdify(t,Tt,"numpy")   # T(t)
    TtNp=TtF(tNp,vNp,tauNp,t1Np,t2Np,t3Np,t4Np)

    # TtDiffF=sympy.lambdify(t,Tt.diff(t),"numpy")
    # TtDiffNp=TtDiffF(tNp)

    maxTt.append(TtNp.max())
    if (judgeTempgA(tNp,TtNp)[0]<=0):
        validList.append(tauNp)
        if (minTau==None) or (tauNp<minTau):
            minTau=tauNp
            minTauTtNp=TtNp
        print("tau合格","，TtNp.max()=",TtNp.max())
    else:
        print("tau不合格","，TtNp.max()=",TtNp.max())
    print("----------------------------------")



if len(validList)>0:
    print("合格的tau有：", validList)
    print("使用tau=",minTau)
    firstLine=None
    firstFlag=True
    with open("result.csv","r") as f:   # 读取result的第一行，写入时不变
        csvRea=csv.reader(f)
        for _ in csvRea:
            if firstFlag:
                firstLine = _
                firstFlag=False
            else:
                break

    with open("result.csv","w",newline="") as f:
        csvWri=csv.writer(f)
        data=[]
        data.append(firstLine)
        resultNp=numpy.zeros([tNp.shape[0],2],dtype=numpy.float32)
        resultNp[:,0]=tNp
        resultNp[:,1]=minTauTtNp
        data+=resultNp.tolist()
        csvWri.writerows(data)

    plt.figure(1)
    plt.plot(tNp,minTauTtNp)
    plt.text(tNp[numpy.where(minTauTtNp == minTauTtNp.max())[0][0]], minTauTtNp.max(), "maximum=" + str(round(minTauTtNp.max(), 2)))
    plt.scatter(tNp[numpy.where(minTauTtNp == minTauTtNp.max())[0][0]], minTauTtNp.max())
    plt.xlabel("time (s)")
    plt.ylabel("temperature (°C)")
    plt.title("$T(t)|\\tau={}$".format(minTau))
    # plt.savefig("figure1.png",dpi=400)

    TatNp = TatF(tNp, vNp, tauNp, t1Np, t2Np, t3Np, t4Np)
    plt.figure(2)
    plt.plot(tNp,TatNp)
    plt.grid(True)  # 显示网格线
    plt.xlabel("time (s)")
    plt.ylabel("temperature (°C)")
    plt.title("$T_a(t)$")
    # plt.savefig("figure3.png",dpi=400)

    plt.figure(3)
    # TtDiffF=sympy.lambdify(t,TtF.diff(t),"numpy")
    TtDiffNp=numpy.diff(minTauTtNp)
    plt.plot(tNp[:-1],TtDiffNp)
    plt.xlabel("time (s)")
    # plt.ylabel("temperature (°C)")
    plt.title("$\\frac{\\mathrm{d}T(t)}{\\mathrm{d}t}|\\tau="+str(minTau)+"$")
    # plt.savefig("figure3.png",dpi=400)

    if len(maxTt)>1:
        maxTtDiff=numpy.diff(numpy.array(maxTt))
        allTauDiff=numpy.diff(numpy.array(allTau))
        tauTDiff=maxTtDiff/allTauDiff
        plt.figure(4)
        plt.plot(allTau,maxTt)
        plt.legend(["Average Gradient = "+str(round(tauTDiff.mean(),2))])
        plt.xlabel("$\\tau$")
        plt.ylabel("max{$T(t)$}")
        plt.title("$T_{\\max}(\\tau)$")


    heatNum=3
    heatTime=dis2Oven+(heatNum-1)*ovenLength+(heatNum-1)*airLength+sympy.Rational(1,2)*ovenLength
    heatIndex=int(round(heatTime/timeIntv))
    print("小温区",str(heatNum),"中点的温度=",str(minTauTtNp[heatIndex]),"°C")

    heatNum=6
    heatTime=dis2Oven+(heatNum-1)*ovenLength+(heatNum-1)*airLength+sympy.Rational(1,2)*ovenLength
    heatIndex=int(round(heatTime/timeIntv))
    print("小温区",str(heatNum),"中点的温度=",str(minTauTtNp[heatIndex]),"°C")

    heatNum=7
    heatTime=dis2Oven+(heatNum-1)*ovenLength+(heatNum-1)*airLength+sympy.Rational(1,2)*ovenLength
    heatIndex=int(round(heatTime/timeIntv))
    print("小温区",str(heatNum),"中点的温度=",str(minTauTtNp[heatIndex]),"°C")

    heatTime=dis2Oven+(8)*ovenLength+(8-1)*airLength+sympy.Rational(1,2)*ovenLength
    heatIndex=int(round(heatTime/timeIntv))
    print("小温区",str(8),"结束的温度=",str(minTauTtNp[heatIndex]),"°C")
else:
    print("没有合格的tau值")