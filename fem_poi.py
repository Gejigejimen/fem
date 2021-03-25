import numpy as np

import numpy as np
 
def gaussian_elimination(A, b):
    try:
        n = len(b)
        # 前進消去を行う
        for i in range(n):
            # pivot選択機能---------------------------------------------------------------------------------
            order = np.argmax(np.abs(A[i:, i]))  # pivot列i行目以上で最も絶対値の高い行を検索
            temp_A = A[i + order].copy()         # Aについて、最も絶対値の高い行成分を抽出し、別メモリで一時保持する
            temp_b = b[i + order].copy()         # bについて、最も絶対値の高い行成分を抽出し、別メモリで一時保持する
            A[i + order] = A[i]                  # Aについて、一時保持した行にi行目成分を置換する
            A[i] = temp_A                        # Aについて、保持しておいた行をi行目に置換する
            b[i + order] = b[i]                  # bについて、一時保持した行にi行目成分を置換する
            b[i] = temp_b                        # bについて、保持しておいた行をi行目に置換する
            # ----------------------------------------------------------------------------------------------
            pivot = A[i, i]                      # 対角成分をpivotに代入
 
            # pivotが小さい時はエラーとする
            if np.abs(pivot) < 1e-6:
                print('pivot=', pivot)
                raise ZeroDivisionError
 
            A[i] = A[i] / pivot                  # pivotで係数行列を割り、A[i,i]を1にする
            b[i] = b[i] / pivot                  # 定数ベクトルもpivotで割り同値変形する
 
            # i行目の定数倍をi+1行目以降から引くループ
            for j in range(i+1, n):
                p = A[j, i]                      # i+1行目以降i列の数値を格納
                A[j] -= p * A[i]                 # 係数行列のi+1行目からi行目の定数倍を引く
                b[j] -= p * b[i]                 # 定数ベクトルのi+1行目からi行目の定数倍を引く
 
        # 後退代入を行う
        x = np.zeros(n)                          # 解の入れ物を用意
        for i in reversed(range(n)):             # 最終行から後退処理する
            x[i] = b[i] / A[i, i]                # 解を求める
            for j in range(i):
                b[j] -= A[j, i] * x[i]           # 解が求まった列分bの値を上から更新する
        print('Normal termination')
    except:
        print('Error termination : pivot too small.')
        x = np.nan
    return A,x,b

A=	np.array([
	[2,-1,-1,0],
	[-1,2,0,-1],
	[-1,0,2,-1],
	[0,-1,-1,2]
	])

b=np.array(
	[2,1,1,2]
	)


result = gaussian_elimination(A,b)
A = result[0]
x = result[1]/(24/2)
print("2要素結果")
print("全体マトリクス")
print(A)
print("解")
print(x[1])

A=np.array([
	[2,-1,0,-1,0,0,0,0,0],
	[-1,4,-1,0,-2,0,0,0,0],
	[0,-1,2,0,0,-1,0,0,0],
	[-1,0,0,4,-2,0,-1,0,0],
	[0,-2,0,-2,8,-2,0,-2,0],
	[0,0,-1,0,-2,4,0,0,-1],
	[0,0,0,-1,0,0,2,-1,0],
	[0,0,0,0,-2,0,-1,4,-1],
	[0,0,0,0,0,-1,0,1,2]
])
b=np.array(
	[2,3,1,3,6,3,1,3,2]
)
result = gaussian_elimination(A,b)
A = result[0]
x = result[1]/(1599/125)
print("8要素結果")
print("全体マトリクス")
print(A)
print("解")
print(x[1])
