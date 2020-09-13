# -*- coding: utf-8 -*-
"""
@author:xuyuntao
@time:2020/9/11 21:40
@email:xuyuntao@189.cn
"""
import numpy
import csv
# import math
import sympy
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
vValue=sympy.Rational(7,6)#cm/s
k=sympy.S(0)
tauValue=sympy.S(108)
ovenTemp=[sympy.S(_) for _ in [25,175,195,235,255,25]]
airIntAdj=[sympy.Rational(3,2), sympy.Rational(2,10), sympy.Rational(4,10), sympy.Rational(2,10), sympy.Rational(23,10)]#cm
ovenLength=sympy.Rational(305,10)#cm
airLength=sympy.S(5)#cm
dis2Oven=sympy.S(25)#cm

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

# tau=sympy.var("tau")
# tau=100

tNp=numpy.arange(0,373,0.5,dtype=numpy.float32)
charList=[]
# tempT.args[方程号][0：方程，1：条件].args[第几个条件]._args[0：t，1：数值]
for tau in numpy.arange(75.8,75.89,0.1):
    print("nowTau=",tau)
    C0=7*tau
    temp0=C0*sympy.exp(sympy.Rational(-150,7)/tau)+25-(7*(tau+sympy.Rational(-150,7)))

    C1=(temp0-175)/sympy.exp(sympy.Rational(-150,7)/tau)
    temp1=C1*sympy.exp(sympy.Rational(-1185,7)/tau)+175

    C2=(temp1+615+sympy.Rational(14,3)*(tau-sympy.Rational(1185,7)))/sympy.exp(sympy.Rational(-1185,7)/tau)
    temp2=C2*sympy.exp(sympy.Rational(-1215,7)/tau)-615-(sympy.Rational(14,3)*(tau+sympy.Rational(-1215,7)))

    C3=(temp2-195)/sympy.exp(sympy.Rational(-1215,7)/tau)
    temp3=C3*sympy.exp(sympy.Rational(-1398,7)/tau)+195

    C4=(temp3+1669+sympy.Rational(28,3)*(tau-sympy.Rational(1398,7)))/sympy.exp(sympy.Rational(-1398,7)/tau)
    temp4=C4*sympy.exp(-204/tau)-1669-(sympy.Rational(28,3)*(tau-204))

    C5=(temp4-235)/sympy.exp(-204/tau)
    temp5=C5*sympy.exp(sympy.Rational(-1611,7)/tau)+235

    C6=(temp5+839+sympy.Rational(14,3)*(tau-sympy.Rational(1611,7)))/sympy.exp(sympy.Rational(-1611,7)/tau)
    temp6=C6*sympy.exp(sympy.Rational(-1641,7)/tau)-839-(sympy.Rational(14,3)*(tau+sympy.Rational(-1641,7)))

    C7=(temp6-255)/sympy.exp(sympy.Rational(-1641,7)/tau)
    temp7=C7*sympy.exp(-291/tau)+255

    C8=(temp7-15872+sympy.Rational(-161,3)*(tau-291))/sympy.exp(-291/tau)
    temp8=C8*sympy.exp(sympy.Rational(-2067,7)/tau)+15872-(sympy.Rational(-161,3)*(tau+sympy.Rational(-2067,7)))

    C9=(temp8-25)/sympy.exp(sympy.Rational(-2067,7)/tau)

    t0=C0*sympy.exp(-t/tau)+25-7*(tau-t)
    t1=C1*sympy.exp(-t/tau)+175
    t2=C2*sympy.exp(-t/tau)-615-sympy.Rational(14,3)*(tau-t)
    t3=C3*sympy.exp(-t/tau)+195
    t4=C4*sympy.exp(-t/tau)-1669-sympy.Rational(28,3)*(tau-t)
    t5=C5*sympy.exp(-t/tau)+235
    t6=C6*sympy.exp(-t/tau)-839-sympy.Rational(14,3)*(tau-t)
    t7=C7*sympy.exp(-t/tau)+255
    t8=C8*sympy.exp(-t/tau)+15872-sympy.Rational(-161,3)*(tau-t)
    t9=C9*sympy.exp(-t/tau)+25

    TtList=[]

    for _ in range(10):
        TtList.append((eval("t"+str(_)),tempTList[_][1]))

    Tt=sympy.Piecewise(*TtList).subs(dis,v*t).subs(v,vValue)
    Tt.evalf(8)

    print("cal num result")
    TtNp=None#tNp.copy()
    TtNpF=sympy.lambdify(t,Tt,"numpy")
    TtNp=TtNpF(tNp)
    # for _ in range(tNp.shape[0]):
    #     print(_,"/",tNp.shape[0])
    #     TtNp[_]=Tt.subs(t,tNp[_]).evalf(5)
    char=numpy.sum((TtNp-dataNpFull)**2)
    charList.append([tau,char])
    print("char",char)
print(charList)
plt.plot(tNp,TtNp)
plt.plot(tNp,dataNpFull)