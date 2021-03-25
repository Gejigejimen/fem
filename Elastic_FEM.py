import numpy as np


# 接点座標
Points=np.array([
				[0,0],
				[1,0],
				[0,1],
				[1,1]
				])
print(Points[1][1])
# 接点数
NP = len(Points)
X = np.zeros(NP)
Y = np.zeros(NP)
print("X座標 , Y座標")
for L in range(NP):
	X[L]=Points[L][0]
	Y[L]=Points[L][1]
	print("Point",L,"[",X[L],Y[L],"]")

# 要素数
NE = 2
II = np.zeros(NE,dtype=np.int8)
JJ = np.zeros(NE,dtype=np.int8)
KK = np.zeros(NE,dtype=np.int8)
# 各要素の接点番号格納
Elements_points=np.array([
						[0,3,2],
						[0,1,3]
						],dtype=np.int8)

# x軸方向BCを与える節点番号
BC_points_x = np.array([
						[0,0],
						[2,0]
						])
IBx = np.zeros(len(BC_points_x),dtype=np.int8)
FBx = np.zeros(len(BC_points_x))
for L in range(len(BC_points_x)):
	IBx[L]=BC_points_x[L][0]
	FBx[L]=BC_points_x[L][1]
# y軸方向BCを与える接点番号
BC_points_y = np.array([
						[0,0],
						[1,0]
						])
IBy = np.zeros(len(BC_points_y),dtype=np.int8)
FBy = np.zeros(len(BC_points_y))
for L in range(len(BC_points_y)):
	IBy[L]=BC_points_y[L][0]
	FBy[L]=BC_points_y[L][1]


print("要素番号  ,  節点番号")
for L in range(NE):
	II[L]=int(Elements_points[L][0])
	JJ[L]=int(Elements_points[L][1])
	KK[L]=int(Elements_points[L][2])
	print("Element",L,"[",II[L],JJ[L],KK[L],"]")

Young = 1.92 * (10**4)#ヤング率
nu = 0.2#ポアソン比

# Dマトリクス作成
D = np.zeros((3,3))
D_coef=Young/(1-nu**2)
D=([
	[1,nu,0],
	[nu,1,0],
	[0,0,(1-nu)/2]
	])
D=np.dot(D_coef,D)
print("Dmatrix\n",D)


# 点に作用する荷重 = 荷重/作用する節点数
sigmaAll = 100 #面に作用する荷重
sigmaPoints = 2 #面を構成する節点の数
sigma = sigmaAll/sigmaPoints #一点に作用する荷重
# 0=荷重あり、1=荷重無し
force_elnum = np.array([1,0])
# 1点に作用する荷重ベクトル
points_force = np.array([sigma,0])

GK=np.zeros((2*NP,2*NP))
LL=np.zeros((3))
EF=np.zeros((6))
GF=np.zeros((2*NP))
Right_vec = np.zeros((2*NP))
El_B = np.zeros((NE,3,6))
#接点ごとに計算する
for L in range(NE):
	I = int(II[L])
	J = int(JJ[L])
	K = int(KK[L])
	print("(I,J,K)\n" ,"(", I , J , K,")")
	# 行列式の値S	
	S = (X[J] - X[I]) * (Y[K] - Y[I]) - (X[K] - X[I]) * (Y[J] - Y[I])
	# 要素三角形の面積G
	G = abs(S)/2
	print("要素面積",G)
	# Bマトリクス作成
	B = np.dot(1/(2*G),np.array([
				[Y[J]-Y[K],0,Y[K]-Y[I],0,Y[I]-Y[J],0],
				[0,X[K]-X[J],0,X[I]-X[K],0,X[J]-X[I]],
				[X[K]-X[J],Y[J]-Y[K],X[I]-X[K],Y[K]-Y[I],X[J]-X[I],Y[I]-Y[J]]
				]))
	El_B[L]=B
	print("Bmatrix\n",B)
	# 要素Kマトリクス作成
	EK = np.dot(G,np.dot(B.T,np.dot(D,B)))
	print("要素剛性マトリクスEK\n",EK)
	# 要素マトリクスEKの行番号と全体係数行列GKの行番号
	LL[0]=I
	LL[1]=J
	LL[2]=K
	# 全体行列への足し合わせ
	for IE in range(3):
		IT = int(LL[IE])
		for JE in range(3):
			JT = int(LL[JE])
			GK[2*IT,2*JT]+=EK[2*IE,2*JE]
			GK[2*IT,2*JT+1]+=EK[2*IE,2*JE+1]
			GK[2*IT+1,2*JT]+=EK[2*IE+1,2*JE]
			GK[2*IT+1,2*JT+1]+=EK[2*IE+1,2*JE+1]

	print("全体剛性マトリクス\n",GK)
	
	#要素荷重ベクトル作成
	if(force_elnum[L]==0):
		for IE in range(3):	
			EF[2*IE] = points_force[0]
			EF[2*IE+1]=points_force[1]
	print("要素荷重ベクトル\n",EF)
	#全体荷重ベクトルへの足し合わせ
	for IE in range(3):
		IT = int(LL[IE])
		GF[2*IT]+=EF[2*IE]
		GF[2*IT+1]+=EF[2*IE+1]
	print("全体荷重ベクトル\n",GF)
	Right_vec = GF

BCV = np.zeros((2*NP))
#境界条件処理
# x軸方向
for L in range(len(BC_points_x)):
	K = IBx[L]
	for J in range(2*NP):
		GK[2*K,J]=0
	GK[2*K,2*K]=1
	Right_vec[2*K]=FBx[L]
# y軸方向
for L in range(len(BC_points_y)):
	K = IBy[L]
	for J in range(2*NP):
		GK[2*K+1,J]=0
	GK[2*K+1,2*K+1]=1
	Right_vec[2*K+1]=FBy[L]
print("境界条件による右辺ベクトル\n",BCV)
print("最終的な全体剛性マトリクス\n",GK)

print("最終的な右辺ベクトル\n",Right_vec)

def gaussian_elimination(A, b):
    n = len(b)
 
    # 前進消去を行う
    for i in range(n):
        pivot = A[i, i]                 # 対角成分をpivotに代入
        A[i] = A[i] / pivot             # pivotで係数行列を割り、A[i,i]を1にする
        b[i] = b[i] / pivot             # 定数ベクトルもpivotで割り同値変形する
 
        # i行目の定数倍をi+1行目以降から引くループ
        for j in range(i+1, n):
            p = A[j, i]                 # i+1行目以降i列の数値を格納
            A[j] -= p * A[i]            # 係数行列のi+1行目からi行目の定数倍を引く
            b[j] -= p * b[i]            # 定数ベクトルのi+1行目からi行目の定数倍を引く
 
    # 後退代入を行う
    x = np.zeros(n)                     # 解の入れ物を用意
    for i in reversed(range(n)):        # 最終行から後退処理する
        x[i] = b[i] / A[i, i]           # 解を求める
        for j in range(i):
            b[j] -= A[j, i] * x[i]      # 解が求まった列分bの値を上から更新する
    return x

result = gaussian_elimination(GK,Right_vec)
U = np.zeros((2*NP))
U = result

print("解\n",U)

for L in range(NE):
	I = int(II[L])
	J = int(JJ[L])
	K = int(KK[L])
	LL[0]=I
	LL[1]=J
	LL[2]=K
	print("要素番号:",L,"節点番号(I,J,K):" ,"(", I , J , K,")")
	u_vec=np.zeros((6))
	for IE in range(3):
		IT=int(LL[IE])
		u_vec[2*IE]=U[2*IT]
		u_vec[2*IE+1]=U[2*IT+1]
	print(u_vec)
	ep_vec=np.dot(El_B[L],u_vec)
	print("ep_vec\n",ep_vec)
	sig_vec=np.dot(D,ep_vec)
	print("sig_vec\n",sig_vec)
