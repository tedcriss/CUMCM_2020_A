# -*- coding: utf-8 -*-
"""
@author:xuyuntao
@time:2020/9/13 10:58
@email:xuyuntao@189.cn
"""
# -*- coding: utf-8 -*-
"""
@author:xuyuntao
@time:2020/9/12 17:11
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
tNp = numpy.arange(0, fullLength / vNp, 0.1)
tauNp = 47.1
t1Np = 182
t2Np = 203
t3Np = 237
t4Np = 254

validList=[]
allV=[]
allJudge=[]
for _ in range(650,1000):
    vNp=_/600
    allV.append(vNp)
    print("当前v=",vNp)
    tNp=numpy.arange(0,fullLength/vNp,0.1) # 时间矩阵，后续数值计算需使用

    TtNp=TtF(tNp,vNp,tauNp,t1Np,t2Np,t3Np,t4Np)

    judge=judgeTempgA(tNp,TtNp)
    allJudge.append(judge)
    if (judge[0]):
        validList.append(vNp)
        print("v合格，为",vNp,"，此时过217面积为",judge[1])
    else:
        print("v不合格，为", vNp,"，此时过217面积为",judge[1])
        # pass
    print("----------------------------------")



if len(validList)>0:
    print("合格的v有：", validList)
    # print("使用tau=",minTau)
    # TatNpF=TatF

    allVNp=numpy.array(allV,dtype=numpy.float32)
    allJudgeNp=numpy.array(allJudge,dtype=numpy.float32)
    if True:
        fig1=plt.figure(1,dpi=90)
        ax1=fig1.add_subplot(1,1,1)
        ax1.scatter(allVNp,allJudgeNp[:,0]<=0,label = 'Valid',color="orange")
        ax1.set_xlabel("Velocity (cm/s)")
        ax1.set_ylabel("Valid")
        ax1.legend(loc=2)
        ax1_2=ax1.twinx()
        ax1_2.plot(allVNp,allJudgeNp[:,1],label = 'Area')
        ax1_2.set_ylabel("Area")
        ax1_2.legend(loc=1)
        ax1.set_title("Valid Velocity and Area of It")

        fig2=plt.figure(2,dpi=70)
        ax2 = fig2.add_subplot(2, 2, 1)
        ax2.plot(allVNp,allJudgeNp[:,2],label = 'Max{Gradient}')
        ax2.set_title("Max Gradient")
        ax2.set_xlabel("Velocity (cm/s)")
        ax3 = fig2.add_subplot(2, 2, 2)
        ax3.plot(allVNp,allJudgeNp[:,3],label = 'Max{Temperature}')
        ax3.set_title("Max Temperature")
        ax3.set_xlabel("Velocity (cm/s)")
        ax3.set_ylabel("°C")
        ax4 = fig2.add_subplot(2, 2, 3)
        ax4.plot(allVNp,allJudgeNp[:,4],label = 'Time of Temperature Bigger Than 217°C')
        ax4.set_title("Time of Temperature Bigger Than 217°C")
        ax4.set_xlabel("Velocity (cm/s)")
        ax4.set_ylabel("s")
        ax5 = fig2.add_subplot(2, 2, 4)
        ax5.plot(allVNp,allJudgeNp[:,5],label = 'Time of Temperature in 150°C~190°C')
        ax5.set_title("Time of Temperature in 150°C~190°C")
        ax5.set_xlabel("Velocity (cm/s)")
        ax5.set_ylabel("s")

    firstLine=["传送速度","是否合格（布尔值）","是否合格（负为合格，正为不合格）","过217面积","最大斜率","峰值温度","150°C~190°C时间","大于217°C时间"]
    writeDataNp=numpy.zeros([allVNp.shape[0],2+allJudgeNp.shape[1]],dtype=numpy.float32)
    writeDataNp[:,0]=allVNp
    writeDataNp[:,1]=allJudgeNp[:,0]<=0
    writeDataNp[:,2:]=allJudgeNp
    with open("question_2.csv","w",newline="") as f:
        csvWri=csv.writer(f)
        csvWri.writerow(firstLine)
        csvWri.writerows(writeDataNp)

    # plt.figure(1)
    # plt.plot(tNp,minTauTtNp)
    # plt.text(tNp[numpy.where(minTauTtNp == minTauTtNp.max())[0][0]], minTauTtNp.max(), "maximum=" + str(round(minTauTtNp.max(), 2)))
    # plt.scatter(tNp[numpy.where(minTauTtNp == minTauTtNp.max())[0][0]], minTauTtNp.max())
    # plt.xlabel("time (s)")
    # plt.ylabel("temperature (°C)")
    # plt.title("$T(t)|\\tau={}$".format(minTau))
    # # plt.savefig("figure1.png",dpi=400)
    #
    #
    # plt.figure(2)
    # plt.plot(tNp,TatNp)
    # plt.grid(True)  # 显示网格线
    # plt.xlabel("time (s)")
    # plt.ylabel("temperature (°C)")
    # plt.title("$T_a(t)$")
    # # plt.savefig("figure3.png",dpi=400)
    #
    # plt.figure(3)
    # TtDiffF=sympy.lambdify(t,Tt.diff(t),"numpy")
    # TtDiffNp=TtDiffF(tNp)
    # plt.plot(tNp,TtDiffNp)
    # plt.xlabel("time (s)")
    # # plt.ylabel("temperature (°C)")
    # plt.title("$\\frac{\\mathrm{d}T(t)}{\\mathrm{d}t}|\\tau="+str(minTau)+"$")
    # # plt.savefig("figure3.png",dpi=400)
    #
    # maxTtDiff=numpy.diff(maxTt)
    # allTauDiff=numpy.diff(allTau)
    # tauTDiff=maxTtDiff/allTauDiff
    # plt.figure(4)
    # plt.plot(allTau,maxTt)
    # plt.legend(["Average Gradient = "+str(round(tauTDiff.mean(),2))])
    # plt.xlabel("$\\tau$")
    # plt.ylabel("max{$T(t)$}")
    # plt.title("$T_{\\max}(\\tau)$")

else:
    print("没有合格的v值")