import numpy as np

def gaussian_elimination(A,b):
	n = len(b)
	print(b)
	# 前進消去法
	for i in range(1,n):
		pivot = A[i,i]	# pivot : 対角成分
		A[i] = A[i] / pivot #pivotで係数行列を割り A[i,i]を1にする
		b[i] = b[i] / pivot #変数ベクトルをpivotで割る
		#i行目の定数倍をi+1行目以降から引く
		for j in range(i+1,n):
			p = A[j,i]	# i列のi+1行目行目以降(対角成分から下)の値を格納
			A[j] -= p*A[i] # j列全体をi列の値xj列のi直下成分で引く
			b[j] -= p*b[i]
	#後退代入を行う
	x = np.zeros(n)	#解格納
	for i in range(n-1,0,-1): # nから1行目までの処理
		x[i] = b[i]/A[i,i]
		for j in range(i):
			b[j] -= A[j,i]*x[i]
	return A,x


# データ入力
print("接点数")
NP = int(input())
NP+=1
X = np.zeros(NP)
Y = np.zeros(NP)
print("X座標 , Y座標")
for L in range(1,NP):
	print(L)
	nums = [float(e) for e in input().split()]
	X[L] = nums[0]
	Y[L] = nums[1]
	nums.clear()

print("要素数")
NE = int(input())
NE+=1
II = np.zeros(NE)
JJ = np.zeros(NE)
KK = np.zeros(NE)
print("頂点の番号")
for L in range(1,NE):
	print(L)
	nums = [float(e) for e in input().split()]
	II[L] = nums[0]
	JJ[L] = nums[1]
	KK[L] = nums[2]
	nums.clear()

print("境界条件を与える点数")
NB = int(input())
NB+=1
IB = np.zeros(NB)
FB = np.zeros(NB)

print("接点番号 , 関数値")
for L in range(1,NB):
	nums = [float(e) for e in input().split()]
	IB[L] = nums[0]
	FB[L] = nums[1]
	nums.clear()



# 係数行列の作成

N = int(NP)
A = np.zeros((N,N))
U = np.zeros((N))
F = np.zeros((N))
LL = np.zeros((3+1))
EK = np.zeros((3+1,3+1))
B = np.zeros((2+1,3+1))

for L in range(NE):
	I = int(II[L])
	J = int(JJ[L])
	K = int(KK[L])
	print("(I,J,K)\n" ,"(", I , J , K,")")
	S = (X[J] - X[I]) * (Y[K] - Y[I]) - (X[K] - X[I]) * (Y[J] - Y[I])
	G = abs(S)
	B[1,1] = (Y[J] - Y[K])/S
	B[1,2] = (Y[K] - Y[I])/S
	B[1,3] = (Y[I] - Y[J])/S
	B[2,1] = (X[K] - X[J])/S
	B[2,2] = (X[I] - X[K])/S
	B[2,3] = (X[J] - X[I])/S

	for IE in range(1,3+1):
		for JE in range(1,3+1):
			EK[IE,JE] = G*(B[1,IE]*B[1,JE] + B[2,IE]*B[2,JE])
	LL[1]=I
	LL[2]=J
	LL[3]=K

	for IE in range(1,3+1):
		IT = int(LL[IE])
		for JE in range(1,3+1):
			JT = int(LL[JE])
			A[IT,JT]=A[IT,JT]+EK[IE,JE]

for L in range(1,NB):
	K =int(IB[L])
	for J in range(1,N):
		A[K,J]=0
	A[K,K]=1
	F[K]=FB[L]

# 連立方程式を解く

result = gaussian_elimination(A,F)
A = result[0]
U = result[1]

print("計算結果")
print("節点番号 , 関数値")
for L in range(1,NP):
	print("  ",L,"   ,   ",U[L])
