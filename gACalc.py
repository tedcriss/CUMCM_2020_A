# -*- coding: utf-8 -*-
"""
@author:xuyuntao
@time:2020/9/12 23:45
@email:xuyuntao@189.cn
"""
import numpy
import geatpy as ea # import geatpy
import matplotlib.pyplot as plt
from gAProblem import minTemp3 # 导入自定义问题接口
from judgeTempFunction import judgeTempgA
from calTt import TtF

if __name__ == '__main__':
    """================================实例化问题对象============================"""
    problem = minTemp3() # 生成问题对象
    """==================================种群设置==============================="""
    Encoding = 'RI'        # 编码方式，采用实整数编码
    NIND = 50             # 种群规模
    Field = ea.crtfld(Encoding, problem.varTypes, problem.ranges, problem.borders) # 创建区域描述器
    population = ea.Population(Encoding, Field, NIND) # 实例化种群对象（此时种群还没被初始化，仅仅是完成种群对象的实例化）
    """================================算法参数设置============================="""
    myAlgorithm = ea.soea_SEGA_templet(problem, population) # 实例化一个算法模板对象
    myAlgorithm.MAXGEN = 400 # 最大进化代数
    myAlgorithm.mutOper.Pm = 0.5 # 变异概率
    myAlgorithm.drawing = 2 # 设置绘图方式（0：不绘图；1：绘制结果图；2：绘制目标空间过程动画；3：绘制决策空间过程动画）
    """===========================调用算法模板进行种群进化======================="""
    [population, obj_trace, var_trace] = myAlgorithm.run() # 执行算法模板
    population.save() # 把最后一代种群的信息保存到文件中
    """===============================输出结果及绘图============================"""
    # 输出结果
    best_gen = numpy.argmin(problem.maxormins * obj_trace[:, 1]) # 最优代索引值
    best_ObjV = numpy.min(obj_trace[:, 1])    # 最优代 值
    print("最优解：",var_trace[best_gen,:])
    print("最优解的过217面积：",best_ObjV)
    objvDiff=numpy.diff(obj_trace[:, 1])
    objvPer=objvDiff/obj_trace[:-1, 1]

    v_best=var_trace[best_gen,0]/60
    t1_best=var_trace[best_gen,1]
    t2_best=var_trace[best_gen,2]
    t3_best=var_trace[best_gen,3]
    t4_best=var_trace[best_gen,4]
    fullLength = 25 + 25 + 11 * 30.5 + 10 * 5   # 传送带总长度
    tNp = numpy.arange(0, fullLength / v_best, 0.1)   # 时间矩阵
    TtNp = TtF(tNp, v_best, 45, t1_best, t2_best, t3_best, t4_best)  # 最优解的T(t)曲线
    print("最优解检验结果：",judgeTempgA(tNp,TtNp))


    plt.figure(3)
    plt.plot(obj_trace[:, 1])
    plt.xlabel("Generation")
    plt.ylabel("Area")
    plt.title("Area-Gen Curve")

    plt.figure(4)
    plt.plot(numpy.abs(objvPer))
    plt.xlabel("Generation")
    plt.ylabel("Percentage")
    plt.title("S Change Percentage of Each Generation")

    plt.figure(5)
    plt.plot(tNp,TtNp)
    plt.xlabel("time (s)")
    plt.ylabel("T(t) (°C)")
    plt.axhline(y=217, ls=":")  # 添加水平直线
    plt.axvline(x=tNp[TtNp==TtNp.max()][0],ls=":")
    plt.text(1,217+3,"T=217°C")
    plt.title("Best Generation's Temperature Curve")

    plt.show()