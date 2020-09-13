# -*- coding: utf-8 -*-
"""
@author:xuyuntao
@time:2020/9/12 10:38
@email:xuyuntao@189.cn
"""
import numpy
import csv
# import math
import sympy
from judgeTempFunction import judgeTempNp
import matplotlib.pyplot as plt
sympy.init_printing()


# disInt=0.01#s
vValue=sympy.Rational(78,60)#cm/s
k=sympy.S(0)
tauValue=sympy.S(108)
ovenTemp=[sympy.S(_) for _ in [25,173,198,230,257,25]]
airIntAdj=[sympy.Rational(3,2), sympy.Rational(2,10), sympy.Rational(4,10), sympy.Rational(2,10), sympy.Rational(150,10)]#cm
ovenLength=sympy.Rational(305,10)#cm
airLength=sympy.S(5)#cm
dis2Oven=sympy.S(25)#cm
coolLength=ovenLength * 2
timeIntv=0.5

# 转换点
p0=[0,ovenTemp[0]]
p1=[dis2Oven + airIntAdj[0] * k, ovenTemp[1]]
p2=[p1[0] - airIntAdj[0] * k + (5 * ovenLength + (5 - 1) * airLength) - airIntAdj[1] / 2 * k, ovenTemp[1]]
p3=[p2[0] + airLength + airIntAdj[1] * k, ovenTemp[2]]
p4=[p3[0] - airIntAdj[1]/2 * k + ovenLength - airIntAdj[2] / 2 * k, ovenTemp[2]]
p5=[p4[0] + airLength + airIntAdj[2] * k, ovenTemp[3]]
p6=[p5[0] - airIntAdj[2]/2 * k + ovenLength - airIntAdj[3] / 2 * k, ovenTemp[3]]
p7=[p6[0] + airLength + airIntAdj[3] * k, ovenTemp[4]]
p8=[p7[0] - airIntAdj[3]/2 * k + ovenLength * 2 + airLength - airIntAdj[4] / 2 * k, ovenTemp[4]]
p9=[p8[0] + airLength + airIntAdj[4] * k + coolLength, ovenTemp[5]]
p10=[p9[0] - airIntAdj[4]/2 * k + ovenLength * 2 + airLength + dis2Oven - coolLength, ovenTemp[5]]
# changeIndex=[]#切换点索引

t=sympy.var("t")
v=sympy.var("v")
tau=sympy.var("tau")
dis = sympy.var("dis")
TatList = []
for _ in range(11 - 1):
    dis1=eval("p"+str(_))[0]
    dis2=eval("p"+str(_+1))[0]
    temp1=eval("p"+str(_))[1]
    temp2=eval("p"+str(_+1))[1]
    # print(dis1,dis2,temp1,temp2)
    equa=(dis-dis2)*sympy.Rational(1,dis1-dis2)*(temp1-temp2)+temp2
    condition=sympy.And(dis>=dis1,dis<=dis2)
    # sympy.pprint(equa)
    TatList.append((equa,condition))
TatList.append((sympy.S(0),True))
Tat=sympy.Piecewise(*TatList).subs(dis,v*t)
Tat=Tat.subs(v,vValue) # Ta(t)
tDomain=sympy.Interval(0,Tat.args[-2][1].args[1]._args[1],True,True)  # 时间定义域

tNp=numpy.arange(0,float(tDomain.args[1].evalf(15)),timeIntv,dtype=numpy.float32) # 时间矩阵，后续数值计算需使用

validList=[]
minTau=None
# validTtList=[]
minTauTtNp=None
allTau=[]
maxTt=[]
for _ in range(450,455,10):
    tau=tau.subs(tau,sympy.Rational(_,10))
    allTau.append(float(tau))
    # 提取条件：Tat.args[方程号][1：条件].args[第几个条件]._args[0：t，1：数值]
    # 提取公式：Tat.args[方程号][0：公式].as_poly().coeffs() -> 其中，0为a，1为b
    print("当前tau=",tau.evalf(4))

    T0=25
    C0=((T0-Tat.args[0][0].as_poly().coeffs()[1])+\
       Tat.args[0][0].as_poly().coeffs()[0]*(tau-Tat.args[0][1].args[0]._args[1]))\
       /sympy.exp(-1*Tat.args[0][1].args[0]._args[1]/tau)
    temp0=C0*sympy.exp(-1*Tat.args[0][1].args[1]._args[1]/tau)+\
        Tat.args[0][0].as_poly().coeffs()[1]\
        -Tat.args[0][0].as_poly().coeffs()[0]*(tau-Tat.args[0][1].args[1]._args[1])

    C1=(temp0-Tat.args[1][0])\
       /sympy.exp(-1*Tat.args[1][1].args[0]._args[1]/tau)
    temp1=C1*sympy.exp(-1*Tat.args[1][1].args[1]._args[1]/tau)+\
          Tat.args[1][0]

    C2=((temp1-Tat.args[2][0].as_poly().coeffs()[1])+\
       Tat.args[2][0].as_poly().coeffs()[0]*(tau-Tat.args[2][1].args[0]._args[1]))\
       /sympy.exp(-1*Tat.args[2][1].args[0]._args[1]/tau)
    temp2=C2*sympy.exp(-1*Tat.args[2][1].args[1]._args[1]/tau)+\
         Tat.args[2][0].as_poly().coeffs()[1]\
         -Tat.args[2][0].as_poly().coeffs()[0]*(tau-Tat.args[2][1].args[1]._args[1])

    C3=(temp2-Tat.args[3][0])\
       /sympy.exp(-1*Tat.args[3][1].args[0]._args[1]/tau)
    temp3=C3*sympy.exp(-1*Tat.args[3][1].args[1]._args[1]/tau)+\
          Tat.args[3][0]

    C4=((temp3-Tat.args[4][0].as_poly().coeffs()[1])+\
       Tat.args[4][0].as_poly().coeffs()[0]*(tau-Tat.args[4][1].args[0]._args[1]))\
       /sympy.exp(-1*Tat.args[4][1].args[0]._args[1]/tau)
    temp4=C4*sympy.exp(-1*Tat.args[4][1].args[1]._args[1]/tau)+\
         Tat.args[4][0].as_poly().coeffs()[1]\
         -Tat.args[4][0].as_poly().coeffs()[0]*(tau-Tat.args[4][1].args[1]._args[1])

    C5=(temp4-Tat.args[5][0])\
       /sympy.exp(-1*Tat.args[5][1].args[0]._args[1]/tau)
    temp5=C5*sympy.exp(-1*Tat.args[5][1].args[1]._args[1]/tau)+\
          Tat.args[5][0]

    C6=((temp5-Tat.args[6][0].as_poly().coeffs()[1])+\
       Tat.args[6][0].as_poly().coeffs()[0]*(tau-Tat.args[6][1].args[0]._args[1]))\
       /sympy.exp(-1*Tat.args[6][1].args[0]._args[1]/tau)
    temp6=C6*sympy.exp(-1*Tat.args[6][1].args[1]._args[1]/tau)+\
          Tat.args[6][0].as_poly().coeffs()[1]\
          -Tat.args[6][0].as_poly().coeffs()[0]*(tau-Tat.args[6][1].args[1]._args[1])

    C7=(temp6-Tat.args[7][0])\
       /sympy.exp(-1*Tat.args[7][1].args[0]._args[1]/tau)
    temp7=C7*sympy.exp(-1*Tat.args[7][1].args[1]._args[1]/tau)+\
          Tat.args[7][0]

    C8=((temp7-Tat.args[8][0].as_poly().coeffs()[1])+\
       Tat.args[8][0].as_poly().coeffs()[0]*(tau-Tat.args[8][1].args[0]._args[1]))\
       /sympy.exp(-1*Tat.args[8][1].args[0]._args[1]/tau)
    temp8=C8*sympy.exp(-1*Tat.args[8][1].args[1]._args[1]/tau)+\
          Tat.args[8][0].as_poly().coeffs()[1]\
          -Tat.args[8][0].as_poly().coeffs()[0]*(tau-Tat.args[8][1].args[1]._args[1])

    C9=(temp8-Tat.args[9][0])\
       /sympy.exp(-1*Tat.args[9][1].args[0]._args[1]/tau)


    t0=C0*sympy.exp(-t/tau)+Tat.args[0][0].as_poly().coeffs()[1]-Tat.args[0][0].as_poly().coeffs()[0]*(tau-t)
    t1=C1*sympy.exp(-t/tau)+Tat.args[1][0]
    t2=C2*sympy.exp(-t/tau)+Tat.args[2][0].as_poly().coeffs()[1]-Tat.args[2][0].as_poly().coeffs()[0]*(tau-t)
    t3=C3*sympy.exp(-t/tau)+Tat.args[3][0]
    t4=C4*sympy.exp(-t/tau)+Tat.args[4][0].as_poly().coeffs()[1]-Tat.args[4][0].as_poly().coeffs()[0]*(tau-t)
    t5=C5*sympy.exp(-t/tau)+Tat.args[5][0]
    t6=C6*sympy.exp(-t/tau)+Tat.args[6][0].as_poly().coeffs()[1]-Tat.args[6][0].as_poly().coeffs()[0]*(tau-t)
    t7=C7*sympy.exp(-t/tau)+Tat.args[7][0]
    t8=C8*sympy.exp(-t/tau)+Tat.args[8][0].as_poly().coeffs()[1]-Tat.args[8][0].as_poly().coeffs()[0]*(tau-t)
    t9=C9*sympy.exp(-t/tau)+Tat.args[9][0]

    TtList=[]

    for _ in range(10):
        TtList.append((eval("t"+str(_)),TatList[_][1]))

    Tt=sympy.Piecewise(*TtList).subs(dis,v*t).subs(v,vValue)
    # Tt.evalf(8)

    print("开始计算数值解")
    TtNp=None#tNp.copy()

    TtNpF=sympy.lambdify(t,Tt,"numpy")   # T(t)
    TtNp=TtNpF(tNp)

    # TtDiffF=sympy.lambdify(t,Tt.diff(t),"numpy")
    # TtDiffNp=TtDiffF(tNp)

    maxTt.append(TtNp.max())
    if (judgeTempNp(tNp,TtNp)[0]):
        validList.append(tau)
        if (minTau==None) or (tau<minTau):
            minTau=tau
            minTauTtNp=TtNp
        print("tau合格","，TtNp.max()=",TtNp.max())
    else:
        print("tau不合格","，TtNp.max()=",TtNp.max())
    print("----------------------------------")



if len(validList)>0:
    print("合格的tau有：", validList)
    print("使用tau=",minTau)
    TatNpF=sympy.lambdify(t,Tat,"numpy")   # Ta(t)
    TatNp=TatNpF(tNp)
    firstLine=None
    firstFlag=True
    with open("../result.csv", "r") as f:
        csvRea=csv.reader(f)
        for _ in csvRea:
            if firstFlag:
                firstLine = _
                firstFlag=False
            else:
                break

    with open("../result.csv", "w", newline="") as f:
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


    plt.figure(2)
    plt.plot(tNp,TatNp)
    plt.grid(True)  # 显示网格线
    plt.xlabel("time (s)")
    plt.ylabel("temperature (°C)")
    plt.title("$T_a(t)$")
    # plt.savefig("figure3.png",dpi=400)

    plt.figure(3)
    TtDiffF=sympy.lambdify(t,Tt.diff(t),"numpy")
    TtDiffNp=TtDiffF(tNp)
    plt.plot(tNp,TtDiffNp)
    plt.xlabel("time (s)")
    # plt.ylabel("temperature (°C)")
    plt.title("$\\frac{\\mathrm{d}T(t)}{\\mathrm{d}t}|\\tau="+str(minTau)+"$")
    # plt.savefig("figure3.png",dpi=400)

    maxTtDiff=numpy.diff(maxTt)
    allTauDiff=numpy.diff(allTau)
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
    print("小温区",str(heatNum),"中点的温度=",str(TtNp[heatIndex]),"°C")

    heatNum=6
    heatTime=dis2Oven+(heatNum-1)*ovenLength+(heatNum-1)*airLength+sympy.Rational(1,2)*ovenLength
    heatIndex=int(round(heatTime/timeIntv))
    print("小温区",str(heatNum),"中点的温度=",str(TtNp[heatIndex]),"°C")

    heatNum=7
    heatTime=dis2Oven+(heatNum-1)*ovenLength+(heatNum-1)*airLength+sympy.Rational(1,2)*ovenLength
    heatIndex=int(round(heatTime/timeIntv))
    print("小温区",str(heatNum),"中点的温度=",str(TtNp[heatIndex]),"°C")

    heatTime=dis2Oven+(8)*ovenLength+(8-1)*airLength+sympy.Rational(1,2)*ovenLength
    heatIndex=int(round(heatTime/timeIntv))
    print("小温区",str(8),"结束的温度=",str(TtNp[heatIndex]),"°C")
else:
    print("没有合格的tau值")