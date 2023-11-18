# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 11:28:02 2023

@author: TB
"""

import numpy as np
import sys
import re

#sys.argv[1]DME分子数
#sys.argb[2]单个ion分子的原子数
def distance(a,b,x,y,z): #a、b为两个原子的坐标，x，y，z为当前盒子的大小,处理周期性边界条件的问题
    dx = abs(a[0] - b[0])
    if dx > x/2:
        dx = x - dx
    dy = abs(a[1] - b[1])
    if dy > y/2:
        dy = y - dy
    dz = abs(a[2] - b[2])
    if dz > z/2:
        dz = z - dz
    return np.sqrt(dx**2 + dy**2 + dz**2)

count_array = np.zeros((6,6))  #构建一个4配位情况的数组（DME数，ion数，count）

for step in range(0,4500):   #MD轨迹步数
    print("step: " + str(step))
    with open("./system" + str(step) + ".gro","r") as file:
                
            
        
        
        
        
        DME_lines = lines[2:int(sys.argv[1])*16 + 2]
        MG_lines = lines[int(sys.argv[1])*16 + 2: int(sys.argv[1])*16 + 2 + 75]  #75个镁离子
        ion_lines = lines[int(sys.argv[1])*16 + 2 + 75 : -1]
        atom_number = int(lines[1].split()[0]) #原子数目
        box_size = lines[-1].split()
        box_x = float(box_size[0])
        box_y = float(box_size[1])
        box_z = float(box_size[2])

        #DME信息提取
        DME_O = np.zeros((int(sys.argv[1]),2,3))
        for num1 in range(int(sys.argv[1])):
            coord_O = np.zeros((2,3))
            i = num1 * 16
            coord_O[0,:] = [float(value) for value in DME_lines[i+1][20:].strip().split()[:3]]
            coord_O[1,:] = [float(value) for value in DME_lines[i+4][20:].strip().split()[:3]]
            DME_O[num1] = coord_O
        
        #Mg2+信息提取
        MG = np.zeros((75,3))
        for num2 in range(75): 
            MG[num2,:] = [float(value) for value in MG_lines[num2][20:].strip().split()[:3]]
        
        #ion信息提取
        ion_B = np.zeros((150,3))
        for num3 in range(150):
            j = num3 * int(sys.argv[2])
            ion_B[num3,:] = [float(value) for value in ion_lines[j][20:].strip().split()[:3]]
            
        #距离判断配位数
        for i,j in enumerate(MG):
        #    print(i)
            DME_num = 0
            ion_num = 0
            for k,l in enumerate(DME_O):
                dis1 = distance(j,l[0],box_x,box_y,box_z)
                dis2 = distance(j,l[1],box_x,box_y,box_z)
                if dis1 <= 0.4 or dis2 <= 0.4:
                    DME_num += 1
            for m,n in enumerate(ion_B):
                dis3 = distance(j,n,box_x,box_y,box_z)
                if dis3 <= 0.4:
                    ion_num += 1
            count_array[DME_num][ion_num] += 1

count_array = count_array/75/4500
ion_cn = count_array.sum(axis=0)[1]*1 + count_array.sum(axis=0)[2]*2 + count_array.sum(axis=0)[3]*3 + count_array.sum(axis=0)[4]*4 + count_array.sum(axis=0)[5]*5
DME_cn = count_array.sum(axis=1)[1]*1 + count_array.sum(axis=1)[2]*2 + count_array.sum(axis=1)[3]*3 + count_array.sum(axis=1)[4]*4 + count_array.sum(axis=1)[5]*5
print(count_array)
print("sum: " + str(count_array.sum()))
print("ion_配位数: " + "{:.2f}".format(ion_cn))
print("DME_配位数: " + "{:.2f}".format(DME_cn))
for i in range(6):
    for j in range(6):
        print(str(i)+ "-" + str(j) + ": " + "{:.2f}".format(count_array[i][j]*100))

                    
            
            
            
            
            
            
