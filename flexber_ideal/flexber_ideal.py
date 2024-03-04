import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

'''
import sympy as sp
from scipy.special import *
'''
##--変数定義．大文字のみ定義すればいい--

L = 2.54
E = 70*10**9
D = 4.2*10**(-3) #円柱断面形状を想定している．
sigma_y = 210 * 10**6
H,B = 4.2*10**(-3),4.2*10**(-3) #Hは高さ，Bは幅
I = B*H**3/12
Zp = B * H**2/6
#I = np.pi*(D**4)/64
#Zp = D**3/6 #全塑性モーメントの計算に使う．円柱断面形状です

ei = E * I

##--先端たわみ荷重からpの算出--
phi_0_data = np.loadtxt("phi_0.csv" ,delimiter=",", encoding="shift-jis") #先端たわみ角のimport
phi_0_data = phi_0_data * np.pi / 180
nn = len(phi_0_data)
p = np.zeros(nn)
for i in range(nn):
    p[i] = ((1 + np.sin(phi_0_data[i]))*0.5)**0.5
print("進捗度10%|pの導出完了")

##--phi_1の算出---
phi_1_data = np.zeros(nn)
for i in range(nn):
    phi_1_data[i] = np.arcsin((np.sqrt(2)*p[i])**(-1))
print("進捗度20%|phi_1 の導出完了")


##--関数の準備F(phi,k)------------------
def Finte(a1,a2,p):
    n = 500
    d_theta = (a2 - a1) / n
    theta_before = a1
    Fanswer = 0
    
    for i in range(n):
        theta_after = theta_before + d_theta
        Fanswer = Fanswer  + (((1 - (p**2)*(np.sin(theta_before))**2)**(-0.5)) + ((1 - (p**2)*(np.sin(theta_after))**2)**(-0.5)) )*d_theta*0.5
        theta_before = theta_after
    return Fanswer

##--関数の準備E(phi,k)------------------
def Einte(a1,a2,p):
    n = 500
    d_theta = (a2 - a1) / n
    theta_before = a1
    Eanswer = 0
    
    for i in range(n):
        theta_after = theta_before + d_theta
        Eanswer = Eanswer + (((1 - (p**2)*(np.sin(theta_before))**2)**(0.5)) + ((1 - (p**2)*(np.sin(theta_after))**2)**(0.5)))*d_theta*0.5
        theta_before = theta_after
    return Eanswer

'''
##--関数の準備（第1種不完全積分F）
def Fintegrate(x,p):
    a = sp.elliptic_f(x, p**2)
    return a

##--関数の準備（第2種不完全積分E）
def Eintegrate(x,p):
    a = sp.elliptic_e(x, p**2)
    return a

##--関数の準備（第１種完全積分K）
def Kintegrate(p):
    a = ellipk(p**2)
    return a
理由はまったくわからないが，残念ながら使用不可．
'''


##--qの導出
q = np.zeros(nn)
for i in range(nn):
    q[i] = Finte(phi_1_data[i], np.pi*0.5, p[i])
    #q[i] = Kintegrate(p[i])-Fintegrate(phi_1_data[i],p[i])
print("進捗度30%|qの導出完了")

##--Pの導出
P_load = np.zeros(nn)
for i in range(nn):
    P_load[i] = (q[i]**2 * ei) / (L**2)
print("進捗度50%|P_loadの導出完了")

##--Δ水平方向変位の導出の導出
delta_h = np.zeros(nn)
for i in range(nn):
    delta_h[i] = L - ((2 * ei * ((2 * p[i]**2) -1))/P_load[i])**0.5
print("進捗度70%|水平方向変位の導出完了")

##--δ垂直方向変位
delta_v = np.zeros(nn)
for i in range(nn):
    a = (ei/P_load[i])**0.5    
    b = Finte(phi_1_data[i],  np.pi*0.5,  p[i])
    c = Einte(phi_1_data[i],  np.pi*0.5,  p[i])
    delta_v[i] = a * (b - 2 * c)    
    
    #delta_v[i] = a * (Kintegrate(p[i])- Fintegrate(phi_1_data[i],p[i]) - 2 * Eintegrate(np.pi* 0.5,p[i]) - 2 * Eintegrate(phi_1_data[i],p[i]) )
print("進捗度90%|垂直方向変位の導出完了")

##--微小変形理論に基づく垂直方向変位--
delta_v2 = np.zeros(nn)
for i in range(nn):
    delta_v2[i] = (P_load[i] * L**3)/(3 * ei)

##--弾塑性判定についての計算--
Mp = Zp*sigma_y
safety_ratio = np.zeros(nn)
Mx= np.zeros(nn)
for i in range(nn):
    Mx[i] = P_load[i] * (L - delta_h[i])
    safety_ratio[i] = Mp / Mx[i]
    if safety_ratio[i] < 1 and safety_ratio[i-1] > 1:
        P_safety_ratio = P_load[i]

##--サブ計算
epsron_x = np.zeros(nn)
epsron_y = np.zeros(nn)
final_x = np.zeros(nn)
final_y = np.zeros(nn)
nondimention_P = np.zeros(nn)
epsron_y2 = np.zeros(nn)
for i in range(nn):
    epsron_x[i] = (L - delta_h[i])/L
    epsron_y[i] = delta_v[i] / L
    epsron_y2[i] = delta_v2[i]/L
    final_x[i] = L - delta_h[i]
    final_y[i] = -delta_v[i]
    nondimention_P[i] = P_load[i] * L**2/(ei)


##--データ整理・出力
datasheet = np.ndarray((nn,8))
phi_0_data = np.loadtxt("phi_0.csv" ,delimiter=",", encoding="shift-jis")
for i in range(nn):
    datasheet[i,0] = phi_0_data[i]
    datasheet[i,1] = p[i]
    datasheet[i,2] = q[i]
    datasheet[i,3] = P_load[i]
    datasheet[i,4] = delta_h[i]
    datasheet[i,5] = delta_v[i]
    datasheet[i,6] = delta_v2[i]
    datasheet[i,7] = safety_ratio[i]

def fig():
    ##---荷重による変位量を3種類プロット
    fig, ax1 = plt.subplots()
    x = datasheet[:,3]
    y1 = datasheet[:,4]
    y2 = datasheet[:,5]
    y3 = datasheet[:,6]

    ax1.plot(x, y1, 'o', color='r', markersize=4)
    ax1.plot(x, y2, '*', color='b', markersize=4)
    #ax1.plot(x, y3, '+', color='k', markersize=4)
    ax1.axvline(x = P_safety_ratio)

    #ax1.legend(["H_Displacement","V_deflection","V_small_deflection","弾塑性の境界"],prop = {"family" : "MS Gothic"})
    ax1.legend(["H_Displacement","V_deflection","弾塑性の境界"],prop = {"family" : "MS Gothic"})


    ax1.set_title("荷重による変位量の推移 | L = " + str(L) + "[m]",fontname = 'MS Gothic')
    ax1.set_xlabel("Load[N] | P")
    ax1.set_ylabel("Displacement[m]")
    print("進捗度100%|計算完了")
    print("|phi_0_data(deg)|  p  |  q  |P_load[N]|delta_H|delta_V|")
    plt.savefig("荷重による変位量の推移")


    ##---安全率をプロット
    fig, ax2 = plt.subplots()
    x = datasheet[:,3]
    y1 = datasheet[:,7]
    ax2.axhline(y=1)
    ax2.plot(x, y1, 'o', color='r', markersize=4)
    ax2.set_title("荷重による安全率の推移",fontname = 'MS Gothic')
    ax2.set_title("安全率の推移",fontname = 'MS Gothic')
    ax2.set_xlabel("Load[N] | P")
    ax2.set_ylabel("safety_ratio[-]")
    plt.savefig("荷重による安全率の推移")

    ##---ひずみのプロット
    fig, ax3 = plt.subplots()
    y = nondimention_P
    x1 = epsron_x
    x2 = epsron_y
    x3 = epsron_y2
    ax3.plot(x1, y, 'o', color='r', markersize=4)
    ax3.plot(x2, y, '*', color='b', markersize=4)
    ax3.plot(x3, y, '+', color='k', markersize=4)
    ax3.set_xlim(0.1,1)
    ax3.set_ylim(0,10)

    ax3.set_title("ひずみ量の比較",fontname = 'MS Gothic')
    ax3.legend(["(L - Δ)/L","δ/L","small_deflection"])
    ax3.set_xlabel("Diflection[-]")
    ax3.set_ylabel("PL**2/EI|[-]")
    plt.savefig("ひずみ量の比較")

    ##---自由端の移動
    fig, ax4 = plt.subplots()
    x1 = final_x
    y1 = final_y

    ax4.plot(x1, y1, 'o', color='r', markersize=4)
    ax4.set_title("自由端の移動 | L = " + str(L) + "[m]",fontname = 'MS Gothic')
    ax4.set_xlabel("x [m]")
    ax4.set_ylabel("y [m]")
    plt.savefig("自由端の移動")

    print("進捗度100%|計算完了")
    print("画像のプロットを行います")

datasheet_pd = pd.DataFrame(datasheet,columns=["phi0","p","q","P[N]","delta_H","delta_V","small_delta_v","safety_ratio"])
datasheet_pd.to_csv('flexber_ideal_cluclate_data.csv')
fig()
plt.show()
