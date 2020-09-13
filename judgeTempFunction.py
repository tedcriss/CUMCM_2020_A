# -*- coding: utf-8 -*-
"""
@author:xuyuntao
@time:2020/9/12 17:50
@email:xuyuntao@189.cn
"""
import numpy

def cvValue(slope,peak,t1,t2):
    value=1
    if (t1>=40 and t1<=90):
        if (t2>=60 and t2<=120):
            if (slope<=3):
                if (peak>=240 and peak<=250):
                    value=-(slope-3)*(slope-0)*(peak-250)*(peak-240)*(t1-90)*(t1-40)*(t2-120)*(t2-60)
                else:
                    value=(peak-250)*(peak-240)
            else:
                value=(slope-3)*(slope-0)
        else:
            value=(t2-120)*(t2-60)
    else:
        value=(t1-90)*(t1-40)
    return value

def judgeTempNp(tNp,TtNp):
    judge=True
    TtNpDiff=numpy.diff(TtNp)
    tNpDiff=numpy.diff(tNp)
    diff=TtNpDiff/tNpDiff
    gradient3Valid=numpy.abs(diff).max()<=3
    if not (gradient3Valid):
        print("最大斜率不合格，为",numpy.abs(diff).max())
    judge=judge and (gradient3Valid)    # 上下斜率小于3
    maxTempValid=(TtNp.max()<=250) and (TtNp.max()>=240)
    if not (maxTempValid):
        print("峰值温度不合格，为",TtNp.max())
    judge=judge and maxTempValid   # 峰值温度小于250，大于240

    validSection=numpy.where(numpy.logical_and(TtNp>=150,TtNp<=190))[0]
    maxIndex=numpy.where(TtNp==TtNp.max())[0][0]
    upValidSection=numpy.where(validSection<=maxIndex)
    upTime=tNp[validSection[upValidSection]].max()-tNp[validSection[upValidSection]].min()
    upTimeJudge=(upTime>=60) and (upTime<=120)
    if not (upTimeJudge):
        print("150°C~190°C时间不合格，为",upTime)
    judge=judge and upTimeJudge

    validSection_2=numpy.where(TtNp>217)[0]
    upTime_2 = tNp[validSection_2].max() - tNp[validSection_2].min()
    upTimeJudge_2=(upTime_2<=90) and (upTime_2>=40)
    if not (upTimeJudge_2):
        print("大于217°C时间不合格，为",upTime_2)
    judge=judge and upTimeJudge_2

    area=numpy.sum((TtNpDiff[validSection_2]/2+TtNp[validSection_2])*tNpDiff[validSection_2])

    return judge,area,numpy.abs(diff).max(),TtNp.max(),upTime,upTime_2
    # 分别为[是否合格，过217面积，最大斜率，峰值温度，150°C~190°C时间，大于217°C时间]

def judgeTempgA(tNp,TtNp):
    judge=True
    TtNpDiff = numpy.diff(TtNp)   # T(t)前向差分
    tNpDiff = numpy.diff(tNp)     # t前向差分
    diff = TtNpDiff / tNpDiff     # 差分比，即斜率

    ########下为上升区间内 150-190温度的时间
    validSection_1 = numpy.where(numpy.logical_and(TtNp >= 150, TtNp <= 190))[0]
    maxIndex_1 = numpy.where(TtNp == TtNp.max())[0][0]
    upValidSection_1 = numpy.where(validSection_1 <= maxIndex_1)
    upTime_1 = tNp[validSection_1[upValidSection_1]].max() - tNp[validSection_1[upValidSection_1]].min()

    ########## 下为过217时间
    validSection_2 = numpy.where(TtNp > 217)[0]
    if (validSection_2.size>0):
        upTime_2 = tNp[validSection_2].max() - tNp[validSection_2].min()
        maxValueIndex = numpy.where(TtNp == TtNp.max())[0][0]  # 峰值索引值
        validSection_3 = validSection_2[validSection_2 <= maxValueIndex]  # 取出过217索引内峰值左边的值
        area = numpy.sum((TtNpDiff[validSection_3] / 2 + TtNp[validSection_3] - 217) * tNpDiff[validSection_3])
    else:
        upTime_2=0
        area=0

    return cvValue(numpy.abs(diff).max(),TtNp.max(),upTime_2,upTime_1),area,numpy.abs(diff).max(),TtNp.max(),upTime_1,upTime_2
    # 分别为[是否合格（负为合格，正为不合格），过217面积，最大斜率，峰值温度，150°C~190°C时间，大于217°C时间]

def judgeFunction(TtF,v=7/6,t1=175,t2=195,t3=235,t4=255,tau=47.1,timeIntv=0.1):
    """计算输入参数[是否合格（负为合格，正为不合格），过217面积，峰值对称方差，峰值两边点数差]"""
    fullLength = 25 + 25 + 11 * 30.5 + 10 * 5
    tNp = numpy.arange(0, fullLength / v, timeIntv)
    TtNp=TtF(tNp,v,tau,t1,t2,t3,t4)
    # todo:增加计算方差、个数差的函数
    TtNpDiff = numpy.diff(TtNp)   # T(t)前向差分
    tNpDiff = numpy.diff(tNp)     # t前向差分
    diff = TtNpDiff / tNpDiff     # 差分比，即斜率

    ########下为上升区间内 150-190温度的时间
    validSection_1 = numpy.where(numpy.logical_and(TtNp >= 150, TtNp <= 190))[0]
    maxIndex_1 = numpy.where(TtNp == TtNp.max())[0][0]
    upValidSection_1 = numpy.where(validSection_1 <= maxIndex_1)
    upTime_1 = tNp[validSection_1[upValidSection_1]].max() - tNp[validSection_1[upValidSection_1]].min()

    ########## 下为过217时间
    validSection_2 = numpy.where(TtNp > 217)[0]
    upTime_2 = tNp[validSection_2].max() - tNp[validSection_2].min()

    maxValueIndex = numpy.where(TtNp == TtNp.max())[0][0]    # 峰值索引值

    validSection_3 = validSection_2[validSection_2 <= maxValueIndex]   # 取出过217索引内峰值左边的值

    area = numpy.sum((TtNpDiff[validSection_3] / 2 + TtNp[validSection_3]-217) * tNpDiff[validSection_3])   # 计算面积

    ##########下面计算方差
    rightIndex=validSection_2[validSection_2 > maxValueIndex]
    # print("rightIndex",rightIndex)
    leftIndex=validSection_2[validSection_2 < maxValueIndex]
    # print("leftIndex",leftIndex)
    minPeakIndex = min(rightIndex.size, leftIndex.size)
    # print("minPeakIndex",minPeakIndex)

    sigma=numpy.sum((TtNp[rightIndex[:minPeakIndex]]-TtNp[leftIndex[:-minPeakIndex-1:-1]])**2)/minPeakIndex
    # print("sigma",sigma)

    #########下为个数差
    lenDiff=abs(rightIndex.size-leftIndex.size)

    cv=cvValue(numpy.abs(diff).max(), TtNp.max(), upTime_2, upTime_1)
    print("两边点数差:",abs(rightIndex.size-leftIndex.size),"，约束值：",cv,"，过217面积：",area)
    return cv,area,sigma,lenDiff
    # [是否合格（负为合格，正为不合格），过217面积，峰值对称方差，峰值两边点数差]