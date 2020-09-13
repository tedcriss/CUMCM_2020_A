# -*- coding: utf-8 -*-
"""
@author:xuyuntao
@time:2020/9/12 19:47
@email:xuyuntao@189.cn
"""
import numpy
import csv
# import math
import sympy
from judgeTempFunction import judgeTempNp
import matplotlib.pyplot as plt
sympy.init_printing()

T1Value=175
T2Value=195
T3Value=235
T4Value=255
vValue=sympy.Rational(70,60)#cm/s
tauValue=45
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
timeIntv=0.5

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
Tt_Origin = sympy.Piecewise(*TtList)  # T(t)


############################遍历############################
validList=[]    # 合格var列表
# allV=[]
allJudge=[]     # 所有var对应judge矩阵
allVar=[]       # 所有var
# tau=tau.subs(tau,sympy.Rational(450,10))
varName="T1"
for _ in range(160,190):
    var=_#sympy.Rational(_,60)#T2.subs(T2,_)     # 每次循环修改的变量值
    T1Value=var    # 可修改
    allVar.append(var)     # var添加到allVar
    print("当前var=",var)
    # Tat=Tat_Origin#.subs(v,vValue) # Ta(t)
    # tDomain=sympy.Interval(0,Tat.args[-2][1].args[1]._args[1],True,True)  # 时间定义域

    # sympy.print_latex(Tt)
    ########T(t)#########
    Tt=Tt_Origin.subs({T1:T1Value,T2:T2Value,T3:T3Value,T4:T4Value}).subs(tau,tauValue).subs(v,vValue)    # 电路板温度-时间 T(t) 函数，分段

    print("开始计算数值解")
    TtNp=None
    tNp=numpy.arange(0,float(Tat.args[-2][1].args[1]._args[1].subs(v,vValue).evalf(15)),timeIntv,dtype=numpy.float32) # 时间矩阵，后续数值计算需使用
    TtNpF=sympy.lambdify(t,Tt,"numpy")
    TtNp=TtNpF(tNp)
    # plt.plot(tNp,TtNp)

    judge=judgeTempNp(tNp,TtNp)
    allJudge.append(judge)
    # allVar.append在上面
    if (judge[0]):
        validList.append(var)
        print("var合格，为",var,"，此时过217面积为",judge[1])
    else:
        print("var不合格，为", var,"，此时过217面积为",judge[1])
        # pass
    print("----------------------------------")

allVarNp=numpy.array(allVar)
allJudgeNp=numpy.array(allJudge)
saveData=numpy.zeros([allJudgeNp.shape[0],allJudgeNp.shape[1]+1],dtype=numpy.float32)
saveData[:,0]=allVarNp
saveData[:,1:]=allJudgeNp
firstLine=["var","是否合格","过217面积","最大斜率","峰值温度","150°C~190°C时间","大于217°C时间"]
with open(varName+"-judge.csv","w",newline="") as f:
    csvWri=csv.writer(f)
    csvWri.writerow(firstLine)
    csvWri.writerows(saveData)
    print("写入文件成功")
# fig1 = plt.figure(1, dpi=90)
# ax1 = fig1.add_subplot(1, 1, 1)
# ax1.scatter(allVarNp, allJudgeNp[:, 0], label='Valid')
# ax1.set_xlabel("Variable")
# ax1.set_ylabel("Valid")
# ax1.legend(loc=2)
# ax1_2 = ax1.twinx()
# ax1_2.plot(allVarNp, allJudgeNp[:, 1], label='Area')
# ax1_2.set_ylabel("Area")
# ax1_2.legend(loc=1)
# ax1.set_title("Valid Variable and Area of It")