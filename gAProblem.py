# -*- coding: utf-8 -*-
"""
@author:xuyuntao
@time:2020/9/12 23:34
@email:xuyuntao@189.cn
"""
import numpy
import geatpy as ea
from judgeTempFunction import judgeTempNp
from judgeTempFunction import judgeFunction
from calTt import TtF

class minTemp3(ea.Problem): # 继承Problem父类
    def __init__(self):
        name = 'minTemp3' # 初始化name（函数名称，可以随意设置）
        M = 1 # 初始化M（目标维数）
        maxormins = [1] # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = 5 # 初始化Dim（决策变量维数）
        varTypes = [0,0,0,0,0] # 初始化varTypes（决策变量的类型，元素为0表示对应的变量是连续的；1表示是离散的）
        lb = [65,165,185,225,245] # 决策变量下界
        ub = [100,185,205,245,265] # 决策变量上界
        lbin = [1,1,1,1,1] # 决策变量下边界（0表示不包含该变量的下边界，1表示包含）
        ubin = [1,1,1,1,1] # 决策变量上边界（0表示不包含该变量的上边界，1表示包含）
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)

    def aimFunc(self, pop): # 目标函数
        Vars = pop.Phen # 得到决策变量矩阵
        v = Vars[:, [0]]/60
        t1 = Vars[:, [1]]
        t2 = Vars[:, [2]]
        t3 = Vars[:, [3]]
        t4 = Vars[:, [4]]
        objv=numpy.zeros([Vars.shape[0],1],dtype=numpy.float64)
        judgeArray=numpy.zeros([Vars.shape[0],1],dtype=numpy.float64)    # 这里True表示可行，False表示不可行
        # todo：将judgeFunction返回的布尔值换成负的，约负代表约合格，约正代表约不合格
        for _ in range(Vars.shape[0]):
            judgeArray[_,0],objv[_,0]=judgeFunction(TtF,v[_],t1[_],t2[_],t3[_],t4[_])[0:2]
        # fixme：上面可以并行化
        # judgeArray=judgeArray
        pop.ObjV=objv
        # pop.ObjV = 4*x1 + 2*x2 + x3 # 计算目标函数值，赋值给pop种群对象的ObjV属性
        # 采用可行性法则处理约束
        # pop.CV = numpy.hstack([2*x1 + x2 - 1,
        #                 x1 + 2*x3 - 2,
        #                 numpy.abs(x1 + x2 + x3 - 1)])
        pop.CV=judgeArray

    def calReferObjV(self): # 设定目标数参考值（本问题目标函数参考值设定为理论最优值）
        referenceObjV = numpy.array([[500]])
        return referenceObjV

class minTemp4(ea.Problem): # 继承Problem父类
    def __init__(self):
        name = 'minTemp4' # 初始化name（函数名称，可以随意设置）
        M = 3 # 初始化M（目标维数）
        maxormins = [1,1,1] # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = 5 # 初始化Dim（决策变量维数）
        varTypes = [0,0,0,0,0] # 初始化varTypes（决策变量的类型，元素为0表示对应的变量是连续的；1表示是离散的）
        lb = [65,165,185,225,245] # 决策变量下界
        ub = [100,185,205,245,265] # 决策变量上界
        lbin = [1,1,1,1,1] # 决策变量下边界（0表示不包含该变量的下边界，1表示包含）
        ubin = [1,1,1,1,1] # 决策变量上边界（0表示不包含该变量的上边界，1表示包含）
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)

    def aimFunc(self, pop): # 目标函数
        Vars = pop.Phen # 得到决策变量矩阵
        v = Vars[:, [0]]/60
        t1 = Vars[:, [1]]
        t2 = Vars[:, [2]]
        t3 = Vars[:, [3]]
        t4 = Vars[:, [4]]
        objv=numpy.zeros([Vars.shape[0],3],dtype=numpy.float64)
        judgeArray=numpy.zeros([Vars.shape[0],1],dtype=numpy.float64)    # 这里True表示可行，False表示不可行
        # todo：将judgeFunction返回的布尔值换成负的，约负代表约合格，约正代表约不合格
        for _ in range(Vars.shape[0]):
            judgeArray[_,0],objv[_,2],objv[_,1],objv[_,0]=judgeFunction(TtF,v[_],t1[_],t2[_],t3[_],t4[_])[0:4]
        # fixme：上面可以并行化
        # judgeArray=judgeArray
        pop.ObjV=objv
        pop.CV=judgeArray

    def calReferObjV(self): # 设定目标数参考值（本问题目标函数参考值设定为理论最优值）
        referenceObjV = numpy.array([[500]])
        return referenceObjV