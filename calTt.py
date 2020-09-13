# -*- coding: utf-8 -*-
"""
@author:xuyuntao
@time:2020/9/12 22:53
@email:xuyuntao@189.cn
"""
import numpy
import csv
import sympy
from judgeTempFunction import judgeTempgA
from judgeTempFunction import cvValue
from judgeTempFunction import judgeFunction
import time
sympy.init_printing()

# 四个小温区的温度
T1,T2,T3,T4=sympy.var("T1 T2 T3 T4")
# disInt=0.01#s
k=sympy.S(0)
# tauValue=sympy.S(108)
ovenTemp=[sympy.S(25),T1,T2,T3,T4,sympy.S(25)]
airIntAdj=[sympy.Rational(3,2), sympy.Rational(2,10), sympy.Rational(4,10), sympy.Rational(2,10), sympy.Rational(150,10)]#cm
ovenLength=sympy.Rational(305,10)#cm
airLength=sympy.S(5)#cm
dis2Oven=sympy.S(25)#cm
coolLength=ovenLength * 2


##################转换点#######################
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

#######################计算Ta(t)######################
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
    equa=((dis-dis2)*sympy.Rational(1,dis1-dis2)*(temp1-temp2)+temp2).subs(dis,v*t)
    condition=sympy.And(t>=dis1/v,t<=dis2/v)
    # sympy.pprint(equa)
    TatList.append((equa,condition))
TatList.append((sympy.S(0),True))
Tat=sympy.Piecewise(*TatList) # Ta(t)

######################计算T(t)##########################
# 提取条件：Tat.args[方程号][1：条件].args[第几个条件]._args[0：t，1：数值]
# 提取公式：Tat.args[方程号][0：公式].as_poly(t).coeffs() -> 其中，0为a，1为b
T0 = 25
C0 = ((T0 - Tat.args[0][0].as_poly(t).coeffs()[1]) + \
      Tat.args[0][0].as_poly(t).coeffs()[0] * (tau - Tat.args[0][1].args[0]._args[1])) \
     / sympy.exp(-1 * Tat.args[0][1].args[0]._args[1] / tau)
temp0 = C0 * sympy.exp(-1 * Tat.args[0][1].args[1]._args[1] / tau) + \
        Tat.args[0][0].as_poly(t).coeffs()[1] \
        - Tat.args[0][0].as_poly(t).coeffs()[0] * (tau - Tat.args[0][1].args[1]._args[1])

C1 = (temp0 - Tat.args[1][0]) \
     / sympy.exp(-1 * Tat.args[1][1].args[0]._args[1] / tau)
temp1 = C1 * sympy.exp(-1 * Tat.args[1][1].args[1]._args[1] / tau) + \
        Tat.args[1][0]

C2 = ((temp1 - Tat.args[2][0].as_poly(t).coeffs()[1]) + \
      Tat.args[2][0].as_poly(t).coeffs()[0] * (tau - Tat.args[2][1].args[0]._args[1])) \
     / sympy.exp(-1 * Tat.args[2][1].args[0]._args[1] / tau)
temp2 = C2 * sympy.exp(-1 * Tat.args[2][1].args[1]._args[1] / tau) + \
        Tat.args[2][0].as_poly(t).coeffs()[1] \
        - Tat.args[2][0].as_poly(t).coeffs()[0] * (tau - Tat.args[2][1].args[1]._args[1])

C3 = (temp2 - Tat.args[3][0]) \
     / sympy.exp(-1 * Tat.args[3][1].args[0]._args[1] / tau)
temp3 = C3 * sympy.exp(-1 * Tat.args[3][1].args[1]._args[1] / tau) + \
        Tat.args[3][0]

C4 = ((temp3 - Tat.args[4][0].as_poly(t).coeffs()[1]) + \
      Tat.args[4][0].as_poly(t).coeffs()[0] * (tau - Tat.args[4][1].args[0]._args[1])) \
     / sympy.exp(-1 * Tat.args[4][1].args[0]._args[1] / tau)
temp4 = C4 * sympy.exp(-1 * Tat.args[4][1].args[1]._args[1] / tau) + \
        Tat.args[4][0].as_poly(t).coeffs()[1] \
        - Tat.args[4][0].as_poly(t).coeffs()[0] * (tau - Tat.args[4][1].args[1]._args[1])

C5 = (temp4 - Tat.args[5][0]) \
     / sympy.exp(-1 * Tat.args[5][1].args[0]._args[1] / tau)
temp5 = C5 * sympy.exp(-1 * Tat.args[5][1].args[1]._args[1] / tau) + \
        Tat.args[5][0]

C6 = ((temp5 - Tat.args[6][0].as_poly(t).coeffs()[1]) + \
      Tat.args[6][0].as_poly(t).coeffs()[0] * (tau - Tat.args[6][1].args[0]._args[1])) \
     / sympy.exp(-1 * Tat.args[6][1].args[0]._args[1] / tau)
temp6 = C6 * sympy.exp(-1 * Tat.args[6][1].args[1]._args[1] / tau) + \
        Tat.args[6][0].as_poly(t).coeffs()[1] \
        - Tat.args[6][0].as_poly(t).coeffs()[0] * (tau - Tat.args[6][1].args[1]._args[1])

C7 = (temp6 - Tat.args[7][0]) \
     / sympy.exp(-1 * Tat.args[7][1].args[0]._args[1] / tau)
temp7 = C7 * sympy.exp(-1 * Tat.args[7][1].args[1]._args[1] / tau) + \
        Tat.args[7][0]

C8 = ((temp7 - Tat.args[8][0].as_poly(t).coeffs()[1]) + \
      Tat.args[8][0].as_poly(t).coeffs()[0] * (tau - Tat.args[8][1].args[0]._args[1])) \
     / sympy.exp(-1 * Tat.args[8][1].args[0]._args[1] / tau)
temp8 = C8 * sympy.exp(-1 * Tat.args[8][1].args[1]._args[1] / tau) + \
        Tat.args[8][0].as_poly(t).coeffs()[1] \
        - Tat.args[8][0].as_poly(t).coeffs()[0] * (tau - Tat.args[8][1].args[1]._args[1])

C9 = (temp8 - Tat.args[9][0]) \
     / sympy.exp(-1 * Tat.args[9][1].args[0]._args[1] / tau)

t0 = C0 * sympy.exp(-t / tau) + Tat.args[0][0].as_poly(t).coeffs()[1] - Tat.args[0][0].as_poly(t).coeffs()[0] * (
            tau - t)
t1 = C1 * sympy.exp(-t / tau) + Tat.args[1][0]
t2 = C2 * sympy.exp(-t / tau) + Tat.args[2][0].as_poly(t).coeffs()[1] - Tat.args[2][0].as_poly(t).coeffs()[0] * (
            tau - t)
t3 = C3 * sympy.exp(-t / tau) + Tat.args[3][0]
t4 = C4 * sympy.exp(-t / tau) + Tat.args[4][0].as_poly(t).coeffs()[1] - Tat.args[4][0].as_poly(t).coeffs()[0] * (
            tau - t)
t5 = C5 * sympy.exp(-t / tau) + Tat.args[5][0]
t6 = C6 * sympy.exp(-t / tau) + Tat.args[6][0].as_poly(t).coeffs()[1] - Tat.args[6][0].as_poly(t).coeffs()[0] * (
            tau - t)
t7 = C7 * sympy.exp(-t / tau) + Tat.args[7][0]
t8 = C8 * sympy.exp(-t / tau) + Tat.args[8][0].as_poly(t).coeffs()[1] - Tat.args[8][0].as_poly(t).coeffs()[0] * (
            tau - t)
t9 = C9 * sympy.exp(-t / tau) + Tat.args[9][0]

TtList = []  # T(t)函数生成列表

for _ in range(10):
    TtList.append((eval("t" + str(_)), TatList[_][1]))

############T(t)################
Tt = sympy.Piecewise(*TtList)  # T(t,v,tau,t1,t2,t3,t4)


#############T(T) numpy数值计算函数#################
TatF=sympy.lambdify([t,v,tau,T1,T2,T3,T4],Tat,"numpy")
TtF=sympy.lambdify([t,v,tau,T1,T2,T3,T4],Tt,"numpy")


if __name__=="__main__":
    sympy.print_latex(Tt)
    import matplotlib.pyplot as plt
    fullLength=25+25+11*30.5+10*5
    # [ 79.92383957, 166.75268173, 197.07910538, 225 , 258.30924988]

    vNp=78/60
    tNp=numpy.arange(0,fullLength/vNp,0.1)
    tauNp=45
    t1Np=173
    t2Np=198
    t3Np=230
    t4Np=257

    TtNp=TtF(tNp,vNp,tauNp,t1Np,t2Np,t3Np,t4Np)
    # print(judgeTempgA(tNp,TtNp))
    print(judgeFunction(TtF,vNp,t1Np,t2Np,t3Np,t4Np))
    # 分别为[是否合格，过217面积，最大斜率，峰值温度，150°C~190°C时间，大于217°C时间]
    plt.plot(tNp,TtNp)
    plt.axhline(y=217, ls=":")  # 添加水平直线
    plt.axvline(x=tNp[TtNp==TtNp.max()][0],ls=":")
    plt.text(1,217+1,"T=217°C")
    plt.xlabel("time (s)")
    plt.ylabel("temperature (°C)")
    plt.title("Best Solution of 3rd Question")
