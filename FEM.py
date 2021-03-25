import numpy as np

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

# 接点座標
Points=np.array([
				[0,0],
				[0.5,0],
				[1,0],
				[0,0.5],
				[0.5,0.5],
				[1,0.5],
				[0,1],
				[0.5,1],
				[1,1]
				])
print(Points[1][1])
# 接点数
NP = len(Points)
X = np.zeros(NP)
Y = np.zeros(NP)
print("X座標 , Y座標")
for L in range(NP):
	print(L)
	X[L]=Points[L][0]
	Y[L]=Points[L][1]
	print("[",X[L],Y[L],"]")

# 要素数
NE = 8
II = np.zeros(NE,dtype=np.int8)
JJ = np.zeros(NE,dtype=np.int8)
KK = np.zeros(NE,dtype=np.int8)
# 各要素の接点番号格納
Elements_points=np.array([
						[0,1,4],
						[1,2,5],
						[0,3,4],
						[1,4,5],
						[3,7,4],
						[4,5,8],
						[3,6,7],
						[4,7,8]
						],dtype=np.int8)

print("要素番号  ,  接点番号")
for L in range(NE):
	print(L)
	II[L]=int(Elements_points[L][0])
	JJ[L]=int(Elements_points[L][1])
	KK[L]=int(Elements_points[L][2])
	print("[",II[L],JJ[L],KK[L],"]")

# 境界条件を与える接点数
NB = 5
IB = np.zeros(NB,dtype=np.int8)
FB = np.zeros(NB)
# 各接点の境界条件 [接点番号,関数値]
Boundary_points = np.array([
							[2,0],
							[5,0],
							[8,0],
							[7,0],
							[6,0]
							])
print("接点番号 , 関数値")
for L in range(NB):
	print(L)
	IB[L]=Boundary_points[L][0]
	FB[L]=Boundary_points[L][1]
	print("[",IB[L],FB[L],"]")


# 係数行列の作成

N = int(NP)
A = np.zeros((N,N))
U = np.zeros((N))
F = np.zeros((N))
# 各要素の接点番号配列
LL = np.zeros((3))
# 要素ごとの剛性マトリクス
EK = np.zeros((3,3))
# 基底の微分を格納した行列
B = np.zeros((2,3))


#接点ごとに計算する
for L in range(NE):
	I = int(II[L])
	J = int(JJ[L])
	K = int(KK[L])
	print("(I,J,K)\n" ,"(", I , J , K,")")
	print("(X[I],X[J],X[K])\n" ,"(", X[I] , X[J] , X[K],")")
	# 行列式の値S	
	S = (X[J] - X[I]) * (Y[K] - Y[I]) - (X[K] - X[I]) * (Y[J] - Y[I])
	# 要素三角形の面積G
	G = abs(S)/2
	B[0,0] = (Y[J] - Y[K])/S
	B[0,1] = (Y[K] - Y[I])/S
	B[0,2] = (Y[I] - Y[J])/S
	B[1,0] = (X[K] - X[J])/S
	B[1,1] = (X[I] - X[K])/S
	B[1,2] = (X[J] - X[I])/S

	# 要素マトリクスの作成
	for IE in range(3):
		for JE in range(3):
			EK[IE,JE] = G*(B[0,IE]*B[0,JE] + B[1,IE]*B[1,JE])
	# 要素マトリクスEKの行番号と全体係数行列Aの行番号
	LL[0]=I
	LL[1]=J
	LL[2]=K
	# 全体行列への足し合わせ
	for IE in range(3):
		IT = int(LL[IE])
		for JE in range(3):
			JT = int(LL[JE])
			A[IT,JT]=A[IT,JT]+EK[IE,JE]

	# 右辺ベクトルがあるならここで作成する
	EF = np.array([abs(S)/6, abs(S)/6, abs(S)/6])#要素右辺ベクトル
	for IE in range(3):
		IT = int(LL[IE])
		F[IT] = F[IT]+EF[IE]

# 境界条件が設定された箇所は全体ベクトルから排除する
for L in range(NB):
	K =int(IB[L])
	for J in range(N):
		A[K,J]=0
	A[K,K]=1
	F[K]=FB[L]

# 連立方程式を解く
print(A)
print(F)
print(U)

result = gaussian_elimination(A,F)
U = result

print("計算結果")
print("節点番号 , 関数値")
for L in range(NP):
	print("  ",L,"   ,   ",U[L])
