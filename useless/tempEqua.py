# -*- coding: utf-8 -*-
"""
@author:xuyuntao
@time:2020/9/11 18:48
@email:xuyuntao@189.cn
"""
import numpy
import csv
# import math
import sympy
import matplotlib.pyplot as plt
sympy.init_printing()

# disInt=0.01#s
vValue=sympy.Rational(7,6)#cm/s
k=sympy.S(0)
tauValue=sympy.S(108)
ovenTemp=[sympy.S(_) for _ in [25,175,195,235,255,25]]
airIntAdj=[sympy.Rational(3,2), sympy.Rational(2,10), sympy.Rational(4,10), sympy.Rational(2,10), sympy.Rational(23,10)]#cm
ovenLength=sympy.Rational(305,10)#cm
airLength=sympy.S(5)#cm
dis2Oven=sympy.S(25)#cm

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
p9=[p8[0] + airLength + airIntAdj[4] * k, ovenTemp[5]]
p10=[p9[0] - airIntAdj[4]/2 * k + ovenLength * 2 + airLength+dis2Oven, ovenTemp[5]]
# changeIndex=[]#切换点索引

# temp=numpy.zeros([int(p10[0]/disInt)],dtype=numpy.float32)#oC
# dis=temp.copy()#cm
# temp0=25*temp.copy()#oC
# preTemp0=25
# T=sympy.Function("T")
# Ta=sympy.Function("Ta")
t=sympy.var("t")
v=sympy.var("v")
tau=sympy.var("tau")
dis = sympy.var("dis")
# dis1=sympy.var("dis1")
# dis2=sympy.var("dis2")
# temp1=sympy.var("temp1")
# temp2=sympy.var("temp2")
# equa=(dis-dis2)/(dis1-dis2)*(temp1-temp2)+temp2
tempTList = []
for _ in range(11 - 1):
    dis1=eval("p"+str(_))[0]
    dis2=eval("p"+str(_+1))[0]
    temp1=eval("p"+str(_))[1]
    temp2=eval("p"+str(_+1))[1]
    # print(dis1,dis2,temp1,temp2)
    equa=(dis-dis2)*sympy.Rational(1,dis1-dis2)*(temp1-temp2)+temp2
    condition=sympy.And(dis>=dis1,dis<=dis2)
    # sympy.pprint(equa)
    tempTList.append((equa,condition))
tempTList.append((sympy.S(0),True))
tempT=sympy.Piecewise(*tempTList).subs(dis,v*t)
tempT=tempT.subs(v,vValue) # Ta(t)

# tempT.args[方程号][0：方程，1：条件].args[第几个条件]._args[0：t，1：数值]

# 数值化
# timeNp=numpy.arange(0,373.2,0.1)
# tempTNpF=sympy.lambdify(t,tempT,"numpy")
# tempTNp=tempTNpF(timeNp)
# plt.plot(timeNp,tempTNp)
# plt.xlabel("time (s)")
# plt.ylabel("temperature (°C)")
# plt.title("$T_a(t)$ Curve")

# 暂时无用
# dsolve(Eq((Ta(t)-T(t)),tau*T(t).diff(t)))
# equa=sympy.Eq((Ta(t)-T(t)),tau*T(t).diff(t))
# result=sympy.dsolve(equa,T(t))
# result=result.subs({Ta(t):tempT,tau:tauValue}).doit()