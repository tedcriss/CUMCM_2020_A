# -*- coding: utf-8 -*-
"""
@author:xuyuntao
@time:2020/9/12 7:54
@email:xuyuntao@189.cn
"""
import numpy
import csv
# import math
import sympy
from judgeTempFunction import judgeTempNp
import matplotlib.pyplot as plt
sympy.init_printing()


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
dataNpFull=numpy.ones([746],dtype=numpy.float32)*30
dataNpFull[37:]=dataNp[:,1]

# disInt=0.01#s
vValue=sympy.Rational(7,6)#cm/s       # 速度值
k=sympy.S(0)                          # 扩缩倍数
ovenTemp=[sympy.S(_) for _ in [25,175,195,235,255,25]]       # 各温区温度
airIntAdj=[sympy.Rational(3,2), sympy.Rational(2,10), sympy.Rational(4,10),\
           sympy.Rational(2,10), sympy.Rational(150,10)]#cm   # 扩展数值
ovenLength=sympy.Rational(305,10)#cm    # 回焊炉温区长度
airLength=sympy.S(5)#cm                 # 回焊炉各温区间隔
dis2Oven=sympy.S(25)#cm                 # 回焊炉路前后区域长度
coolLength=ovenLength * 2               # 实际冷却区缩减
timeIntv=0.5                            # 计算时间间隔

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

t=sympy.var("t")         # 时间
v=sympy.var("v")         # 速度
tau=sympy.var("tau")     # ？
dis = sympy.var("dis")   # 距离

print("开始计算空气温度")
TatList = []
for _ in range(11 - 1):
    dis1=eval("p"+str(_))[0]
    dis2=eval("p"+str(_+1))[0]
    temp1=eval("p"+str(_))[1]
    temp2=eval("p"+str(_+1))[1]
    # print(dis1,dis2,temp1,temp2)
    equa=(dis-dis2)*sympy.Rational(1,dis1-dis2)*(temp1-temp2)+temp2   # 线性方程，自变量为距离，因变量为空气温度，equa表示空气温度
    condition=sympy.And(dis>=dis1,dis<=dis2)                          # 上方程定义域
    # sympy.pprint(equa)
    TatList.append((equa,condition))
TatList.append((sympy.S(0),True))                                     # 包括所有情况
Tat=sympy.Piecewise(*TatList).subs(dis,v*t)                           # 空气 温度-距离 方程，分段表示
Tat=Tat.subs(v,vValue) # Ta(t)                                        # 变换为空气 温度-时间 方程，分段表示
tDomain=sympy.Interval(0,Tat.args[-2][1].args[1]._args[1],True,True)  # 时间定义域

tNp=numpy.arange(0,float(tDomain.args[1].evalf(15)),timeIntv,dtype=numpy.float32)[:dataNpFull.shape[0]] # 时间矩阵，后续数值计算需使用

charList=[]     # 误差矩阵
minChar=None    # 最小误差值
minTtNp=None    # 最小误差对应T(t)
minTau=None     # 最小误差对应tau值
useTau=sympy.Rational(450,10).round(0)     # 从另一处得到的合适的tau
useChar=None    # 使用tau的误差
useTtNp=None    # 使用tau求出的T(t)
for _ in range(450,455,10):    # 遍历所有tau
    tau=tau.subs(tau,sympy.Rational(_,10))    # 修改tau为分式，方便计算
    # 提取条件：Tat.args[方程号][1：条件].args[第几个条件]._args[0：t，1：数值]
    # 提取公式：Tat.args[方程号][0：公式].as_poly().coeffs() -> 其中，0为a，1为b
    print("当前tau=",tau.evalf(4))

    T0=25           # 给定环境温度，即电路板初始温度
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

    TtList=[]    # T(t)函数生成列表

    for _ in range(10):
        TtList.append((eval("t"+str(_)),TatList[_][1]))

    Tt=sympy.Piecewise(*TtList).subs(dis,v*t).subs(v,vValue)    # 电路板温度-时间 T(t) 函数，分段
    # Tt.evalf(8)

    print("开始计算数值解")
    TtNp=None#tNp.copy()

    TtNpF=sympy.lambdify(t,Tt,"numpy")   # T(t)
    TtNp=TtNpF(tNp)

    TatNpF=sympy.lambdify(t,Tat,"numpy")   # Ta(t)
    TatNp=TatNpF(tNp)

    char=numpy.sum((TtNp-dataNpFull)**2)   # 计算平方差
    if (minChar==None) or (char<=minChar):
        minChar=char
        minTtNp=TtNp
        minTau=tau.evalf(12)
    if (tau==useTau):
        useChar=char
        useTtNp=TtNp
    charList.append([tau,char])
    print("char=",char)
    print("--------------------------------")

print("charList",charList)

plt.figure(1)
plt.plot(tNp,useTtNp)
plt.plot(tNp,dataNpFull)
plt.text(tNp[numpy.where(useTtNp==useTtNp.max())[0][0]],useTtNp.max(),"maximum="+str(round(useTtNp.max(),2)))
plt.scatter(tNp[numpy.where(useTtNp==useTtNp.max())[0][0]],useTtNp.max())
plt.legend(["$T(t)$","Given Data"])
# plt.plot(tNp,TatNp)
# plt.plot(tNp,dataNpFull-minTtNp)
plt.xlabel("time (s)")
plt.ylabel("temperature (°C)")
plt.title("$T(t)$ and Original Temperature, $\\tau={}$".format(useTau.evalf(2)))
# plt.savefig("figure1.png",dpi=400)

plt.figure(2)
charListNp=numpy.array(charList,dtype=numpy.float32)
plt.plot(charListNp[:,0],charListNp[:,1])
plt.text(minTau,minChar,"min{$\\tau$}="+str(minTau.evalf(2)))
plt.text(useTau,useChar,"use $\\tau$="+str(useTau.evalf(2)))
plt.scatter(minTau,minChar)
plt.scatter(useTau,useChar)
plt.xlabel("$\\tau$")
plt.ylabel("char")
plt.title("char")
# plt.savefig("figure2.png",dpi=400)

plt.figure(3)
plt.plot(tNp,TatNp)
plt.grid(True)     # 显示网格线
plt.xlabel("time (s)")
plt.ylabel("temperature (°C)")
plt.title("$T_a(t)$")
# plt.savefig("figure3.png",dpi=400)