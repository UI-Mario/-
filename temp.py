import numpy as np
from math import *

# 读数据方法2
s = np.loadtxt('right.txt')
print(len(s))
# 记得手动删除像点坐标文件里的重复点，
# 懒得写了
q = np.loadtxt('控制点坐标.txt')
# 点号,x,y
Pxy = []
for i in range(len(s)):
	if s[i][0] > 100:
		Pxy.append([s[i][0], s[i][1], s[i][2]])
print(len(Pxy))

# X,Y,Z，与Pxy相对应
ground_list = []
# 用来存点号的
temp = [0]*len(q)
for i in range(len(q)):
	temp[i] = q[i][0]
temp_2 = []
for i in range(len(Pxy)):
	position = temp.index(Pxy[i][0])
	ground_list.append([q[position][2], q[position][3], -q[position][1]])
# print(ground_list)

# x,y
image_list = []
for i in range(len(Pxy)):
	image_list.append([Pxy[i][1], Pxy[i][2]])


# 初值设定
# 内外方位元素(9)+畸变系数(4)
phi = 0.0
omega = 0.0
kappa = 0.0
Xs = 4500
Ys = -150
Zs = -800
# 单位mm
f = 36.0
x0 = 0.0
y0 = 0.0
k1 = 0.0
k2 = 0.0
p1 = 0.0
p2 = 0.0

# 迭代次数
num_ite = 0

# 方程内的矩阵
rotate = np.mat(np.zeros((3, 3)))#旋转矩阵
A = np.mat(np.zeros((len(image_list)*2, 9+4)))
# LL:V=(AX+CX2+DXad-L)=(AX-L),L
L = np.mat(np.zeros((len(image_list)*2,1)))
# 存改正数的
delta = [0]*(9+4)
# 存中误差
m = [0]*(9+4)


# 开始处理
while (True):

	# 旋转矩阵
	a1 = cos(phi)*cos(kappa)-sin(phi)*sin(omega)*sin(kappa)
	a2 = (-1.0) * cos(phi) * sin(kappa) - sin(phi) * sin(omega) * cos(kappa)
	a3 = (-1.0) * sin(phi) * cos(omega)
	b1 = cos(omega) * sin(kappa)
	b2 = cos(omega) * cos(kappa)
	b3 = (-1.0) * sin(omega)
	c1 = sin(phi) * cos(kappa) + cos(phi) * sin(omega) * sin(kappa)
	c2 = (-1.0) * sin(phi) * sin(kappa) + cos(phi) * sin(omega) * cos(kappa)
	c3 = cos(phi) * cos(omega)

	for i in range(0, len(image_list)):
		# 单位转换(像素->毫米)
		x = image_list[i][0]*5.1966/1000
		y = image_list[i][1]*5.1966/1000
		X = ground_list[i][0]
		Y = ground_list[i][1]
		Z = ground_list[i][2]

		rr = (x-x0)*(x-x0)+(y-y0)*(y-y0)

		Xp = a1*(X-Xs) + b1*(Y-Ys) + c1*(Z-Zs);
		Yp = a2*(X-Xs) + b2*(Y-Ys) + c2*(Z-Zs);
		Zp = a3*(X-Xs) + b3*(Y-Ys) + c3*(Z-Zs);


		A[i * 2,0] = 1.0*(a1*f + a3*(x-x0))/Zp
		A[i * 2,1] = 1.0*(b1*f + b3*(x-x0))/Zp
		A[i * 2,2] = 1.0*(c1*f + c3*(x-x0))/Zp
		A[i * 2,3] = (y-y0)*sin(omega) - ((x-x0)*((x-x0)*cos(kappa) - (y-y0)*sin(kappa))/f + f*cos(kappa))*cos(omega)
		A[i * 2,4] = -f*sin(kappa) - (x-x0)*((x-x0)*sin(kappa) + (y-y0)*cos(kappa))/f
		A[i * 2,5] = y-y0

		A[i * 2,6] = (x-x0)/f
		A[i * 2,7] = 1
		A[i * 2,8] = 0

		A[i * 2,9] = -(x-x0)*rr
		A[i * 2,10] = -(x-x0)*rr*rr
		A[i * 2,11] = -rr-2*(x-x0)*(x-x0)
		A[i * 2,12] = -2*(x-x0)*(y-y0)

		A[i * 2 + 1,0] = 1.0*(a2*f + a3*(y-y0))/Zp
		A[i * 2 + 1,1] = 1.0*(b2*f + b3*(y-y0))/Zp
		A[i * 2 + 1,2] = 1.0*(c2*f + c3*(y-y0))/Zp
		A[i * 2 + 1,3] = -(x-x0)*sin(omega) - ((y-y0)*((x-x0)*cos(kappa)-(y-y0)*sin(kappa))/f-f*sin(kappa))*cos(omega)
		A[i * 2 + 1,4] = -f*cos(kappa) - (y-y0)*((x-x0)*sin(kappa) + (y-y0)*cos(kappa))/f
		A[i * 2 + 1,5] = -(x-x0)

		A[i * 2 + 1,6] = (y-y0)/f
		A[i * 2 + 1,7] = 0
		A[i * 2 + 1,8] = 1

		A[i * 2 + 1,9] = -(y-y0)*rr
		A[i * 2 + 1,10] = -(y-y0)*rr*rr
		A[i * 2 + 1,11] = -2*(x-x0)*(y-y0)
		A[i * 2 + 1,12] = -rr-2*(y-y0)*(y-y0)

		detax = (x-x0)*(k1*rr+k2*rr*rr)+p1*(rr+2*(x-x0)*(x-x0))+2*p2*(x-x0)*(y-y0)
		detay = (y-y0)*(k1*rr+k2*rr*rr)+p2*(rr+2*(y-y0)*(y-y0))+2*p1*(x-x0)*(y-y0)

		L[i*2, 0] = x + f*Xp/Zp-x0+detax
		L[i*2+1, 0] = y + f*Yp/Zp-y0+detay

	AT = A.T #转置矩阵
	ATA = AT*A#矩阵相乘
	ATAInv = ATA.I#求逆矩阵
	ATAAT = ATAInv*AT
	delta = ATAAT*L

	# 添加改正数
	Xs += delta[0]
	Ys += delta[1]
	Zs += delta[2]
	phi += delta[3]
	omega += delta[4]
	kappa += delta[5]
	f += delta[6]
	x0 += delta[7]
	y0 += delta[8]
	k1 += delta[9]
	k2 += delta[10]
	p1 += delta[11]
	p2 += delta[12]

	num_ite+=1

	# 限差(看具体情况设置限差)
	if (fabs(delta[3]) < 1e-6) and (fabs(delta[4]) < 1e-6) and (fabs(delta[5]) < 1e-6):
		break
	if num_ite>60:
		print("Error:overtime")
		break

'''
后续处理，各类残差、中误差
'''

print("迭代次数:%f\n"%num_ite)
rotate[0,0] = cos(phi)*cos(kappa) - sin(phi)*sin(omega)*sin(kappa)
rotate[0,1] = (-1.0)*cos(phi)*sin(kappa) - sin(phi)*sin(omega)*cos(kappa)
rotate[0,2] = (-1.0)*sin(phi)*cos(omega)
rotate[1,0] = cos(omega)*sin(kappa)
rotate[1,1] = cos(omega)*cos(kappa)
rotate[1,2] = (-1.0)*sin(omega)
rotate[2,0] = sin(phi)*cos(kappa) + cos(phi)*sin(omega)*sin(kappa)
rotate[2,1] = (-1.0)*sin(phi)*sin(kappa) + cos(phi)*sin(omega)*cos(kappa)
rotate[2,2] = cos(phi)*cos(omega)
 
AX = A*delta
# 残差
v = AX - L
 
vv = 0
for it in v:
    vv += it[0,0]*it[0,0]
m0 = sqrt(vv/(2*len(Pxy)-13))

print('---------------------')
# 存残差的
print_v = np.zeros((len(Pxy), 3))
for i in range(0,len(Pxy)):
	print_v[i][0] = Pxy[i][0]
	print_v[i][1] = v[2*i][0]
	print_v[i][2] = v[2*i+1][0]
	# print("%d   %f   %f"%(Pxy[i][0], v[2*i][0], v[2*i+1][0]))
# 可以排个序来挑出较大残差，这里没有
for i in range(0, len(Pxy)):
	print("点号:%d  vx:%f  vy:%f\n"%(print_v[i][0], print_v[i][1], print_v[i][2]))
print('---------------------')

print("单位权中误差m0=%f\n"%m0)
# print(ATAInv)
for n in range(9+4):
    m[n] = m0 * sqrt(abs(ATAInv[n,n]))

print("\n结果：\n")
print("Xs=%f  m=%f"%(Xs, m[0])+"  单位：mm\n")
print("Ys=%f  m=%f"%(Ys, m[1])+"  单位：mm\n")
print("Zs=%f   m=%f"% (Zs, m[2])+"  单位：mm\n")
print("phi=%f  m=%f"% (phi, m[3])+"  单位：弧度\n")
print("omega=%f    m=%f"% (omega, m[4])+"  单位：弧度\n")
print("kappa=%f    m=%f"%(kappa, m[5])+"  单位：弧度\n")
print("f=%f    m=%f"%(f, m[6])+"  单位：mm\n")
print("x0=%f    m=%f"%(x0, m[7])+"  单位：mm\n")
print("y0=%f    m=%f"%(y0, m[8])+"  单位：mm\n")
print("k1=%e    m=%e\n"%(k1, m[9]))
print("k2=%e    m=%e\n"%(k2, m[10]))
print("p1=%e    m=%e\n"%(p1, m[11]))
print("p2=%e    m=%e\n"%(p2, m[12]))
print("\n旋转矩阵R为：\n")
print(rotate)


# 写文件
with open('后方交会结果.txt','w') as fo:
	fo.writelines("\n像点残差：\n")
	for i in range(0, len(Pxy)):
		fo.writelines("点号:%d  vx:%f  vy:%f\n"%(print_v[i][0], print_v[i][1], print_v[i][2]))
	fo.writelines("   \n")
	fo.writelines("单位权中误差m0=%f\n"%m0)
	fo.writelines("          \n")
	fo.writelines("Xs=%f  m=%f"%(Xs, m[0])+"  单位：mm\n")
	fo.writelines("Ys=%f  m=%f"%(Ys, m[1])+"  单位：mm\n")
	fo.writelines("Zs=%f   m=%f"% (Zs, m[2])+"  单位：mm\n")
	fo.writelines("phi=%f  m=%f"% (phi, m[3])+"  单位：弧度\n")
	fo.writelines("omega=%f    m=%f"% (omega, m[4])+"  单位：弧度\n")
	fo.writelines("kappa=%f    m=%f"%(kappa, m[5])+"  单位：弧度\n")
	fo.writelines("   \n")
	fo.writelines("f=%f    m=%f"%(f, m[6])+"  单位：mm\n")
	fo.writelines("x0=%f    m=%f"%(x0, m[7])+"  单位：mm\n")
	fo.writelines("y0=%f    m=%f"%(y0, m[8])+"  单位：mm\n")
	fo.writelines("     \n")
	fo.writelines("k1=%e    m=%e\n"%(k1, m[9]))
	fo.writelines("k2=%e    m=%e\n"%(k2, m[10]))
	fo.writelines("p1=%e    m=%e\n"%(p1, m[11]))
	fo.writelines("p2=%e    m=%e\n"%(p2, m[12]))



