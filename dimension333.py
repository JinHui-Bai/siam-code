# coding=utf-8

import math
import random
import matplotlib.pyplot as plt
from decimal import *
getcontext().prec= 40
class algorithm():
	def __init__(self,s,theta,gamma,ST,order,fstar):
		self.s = s
		self.theta = theta
		self.gamma = gamma
		self.Z_Ln = [[[Decimal(1)]]]
		self.Z_Lnt = [[[Decimal(1)]]]
		self.ZnX = []
		self.fn = [[[Decimal(0)]]]
		self.fnbar = [[[Decimal(0)]]]
		#self.fstar = fstar
		self.fstar = [[[1/Decimal((j+1)**1.5) for ll in range(k+1)] for k in range(j+1)] for j in range(order+1)]
		self.n = 1
		self.Ln = 0
		self.c_2Lnf = Decimal(1)
		self.ST = ST
		#阶乘常数(2n-1)!!
		self.factoral = 1
		self.coefficient = [[[Decimal(1)]]]
		self.Su = []
	
	def factorial(self):
		S = 1
		for i in range(self.Ln):
			S *= (i+1)
		self.factoral = S
		return S
	
	def update_Z_Ln(self):
		#Ln必须是更新过的Ln
		self.Z_Ln = []
		self.factorial()
		self.c_2Lnf = (2*Decimal(self.Ln)-1)*self.c_2Lnf/Decimal(self.factoral)
		#更新系数矩阵
		Alpha = self.coefficient[-1][:]
		Alpha.append([Decimal(0) for j in range(self.Ln+1)])
		Beta = self.coefficient[-1][:]
		Beta = [ a+[Decimal(0)] for a in Beta]
		Beta.insert(0,[Decimal(0)])
		Gamma = self.coefficient[-1][:]
		Gamma = [ [Decimal(0)]+a for a in Gamma]
		Gamma.insert(0,[Decimal(0)])
		self.coefficient.append([[ a+b+c for a,b,c in zip(alpha,beta,gamma)] for alpha,beta,gamma in zip(Alpha,Beta,Gamma)])
		if self.Ln%2 == 0:
			C_2jf = [Decimal(1)]
			c = Decimal(-1-2*self.Ln)
			C_lambda = [Decimal(1)]
			factorial = Decimal(1)
			for j in range(1,self.Ln//2+1):
				C_2jf.append(2*j*C_2jf[-1])
				c = c+2
				C_lambda.append(c*C_lambda[-1])
			C_2jf.reverse()
			C_lambda.reverse()
			for j in range(self.Ln+1):
				if j%2 == 1:
					self.Z_Ln.append([[ Decimal(0) for ll in range(k+1)] for k in range(j+1)])
					factorial = factorial*j
				else:
					if j > 0:
						factorial = factorial*j
					self.Z_Ln.append([[self.c_2Lnf*a/(C_2jf[j//2]*C_lambda[j//2]*factorial) for a in A] for A in self.coefficient[j]])
		else:
			C_2jf = [Decimal(1)]
			c = Decimal(-1-2*self.Ln)
			C_lambda = [Decimal(1)]
			factorial = Decimal(1)
			for j in range(1,(self.Ln-1)//2+1):
				C_2jf.append(2*j*C_2jf[-1])
				c = c+2
				C_lambda.append(c*C_lambda[-1])
			C_2jf.reverse()
			C_lambda.reverse()
			for j in range(self.Ln+1):
				if j%2 == 0:
					if j > 0:
						factorial = factorial*j
					self.Z_Ln.append([[ Decimal(0) for ll in range(k+1)] for k in range(j+1)])
				else:
					factorial = factorial*j
					self.Z_Ln.append([[self.c_2Lnf*a/(C_2jf[(j-1)//2]*C_lambda[(j-1)//2]*factorial) for a in A] for A in self.coefficient[j]])
					
	def update_Z_Lnt(self):
		self.update_Z_Ln()
		omega_Ln = (Decimal(self.Ln+1))**(-4*self.s)
		self.Z_Lnt.append([[Decimal(0) for k in range(j+1)] for j in range(self.Ln+1)])
		self.Z_Lnt = [[[omega_Ln*aLn + aLnt for aLn,aLnt in zip(alphaLn,alphaLnt)] for alphaLn,alphaLnt in zip(AlphaLn,AlphaLnt)] for AlphaLn,AlphaLnt in zip(self.Z_Ln,self.Z_Lnt)]
		
	def Z_LntX_cal(self,x1,x2,x3):
		self.ZnX = [self.Z_Lnt[0][:]]
		Line = [[Decimal(1)]]
		L = len(self.Z_Lnt)
		for j in range(1,L):
			x2x3i = Line[-1][:]
			x3i = x2x3i[-1]
			x2x3i = [ a*x2 for a in x2x3i]
			x2x3i.append(x3i*x3)
			Line = [[x1*a for a in A] for A in Line]
			Line.append(x2x3i[:])
			self.ZnX.append([[ a*b for a,b in zip(alpha,beta)] for alpha,beta in zip(Line,self.Z_Lnt[j])])
			
	def f_n_cal(self,x1,x2,x3):
		S = self.fn[0][0][0]
		Line = [[Decimal(1)]]
		L = len(self.fn)
		for j in range(1,L):
			x2x3i = Line[-1][:]
			x3i = x2x3i[-1]
			x2x3i = [ x2*a for a in x2x3i]
			x2x3i.append(x3i*x3)
			Line = [[x1*a for a in A] for A in Line]
			Line.append(x2x3i[:])
			S += sum([sum([ a*b for a,b in zip(fnline,xline)]) for fnline,xline in zip(self.fn[j],Line)])
		return S
				
	def f_nbar_cal(self,x1,x2,x3):
		S = self.fnbar[0][0][0]
		Line = [[Decimal(1)]]
		L = len(self.fnbar)
		for j in range(1,L):
			x2x3i = Line[-1][:]
			x3i = x2x3i[-1]
			x2x3i = [ x2*a for a in x2x3i]
			x2x3i.append(x3i*x3)
			Line = [[x1*a for a in A] for A in Line]
			Line.append(x2x3i[:])
			S += sum([sum([ a*b for a,b in zip(fnbarline,xline)]) for fnbarline,xline in zip(self.fnbar[j],Line)])
		return S
		
	def algorithm(self,X,Y):
		coe = self.gamma*(Y - self.f_n_cal(X[0],X[1],X[2]))
		if self.n**self.theta >= (self.Ln+1)**2:
			self.Ln = self.Ln + 1
			self.update_Z_Lnt()
			self.fn.append([ [ Decimal(0) for k in range(j+1)] for j in range(self.Ln+1)])
			self.fnbar.append([ [ Decimal(0) for k in range(j+1)] for j in range(self.Ln+1)])
		self.Z_LntX_cal(X[0],X[1],X[2])
		for j in range(self.Ln+1):
			self.fn[j] = [[ coe*azLnX + afn for azLnX,afn in zip(ZLnX,fnline)] for ZLnX,fnline in zip(self.ZnX[j],self.fn[j])]
			self.fnbar[j] = [[ (Decimal(self.n)*afnbar+afn)/(Decimal(self.n)+1) for afnbar,afn in zip(fnbarl,fnl)] for fnbarl,fnl in zip(self.fnbar[j],self.fn[j])]
		
	def fstar_cal(self,x1,x2,x3):
		#这里Decimal没有加全
		S = Decimal(self.fstar[0][0][0])
		Line = [[Decimal(1)]]
		L = len(self.fstar)
		for j in range(1,L):
			x2x3i = Line[-1][:]
			x3i = x2x3i[-1]
			x2x3i = [ x2*a for a in x2x3i]
			x2x3i.append(x3i*x3)
			Line.append(x2x3i)
			S += sum([sum([ a*b for a,b in zip(fsline,xline)]) for fsline,xline in zip(self.fstar[j],Line)])
		return S
	
	def construct_data(self):
		phi = random.uniform(0,2*math.pi)
		theta = math.acos(random.uniform(-1, 1))
		X = [Decimal(math.sin(theta)*math.cos(phi)),Decimal(math.sin(theta)*math.sin(phi)),Decimal(math.cos(theta))]
		Y = self.fstar_cal(X[0],X[1],X[2])
		Z = Y + Decimal(random.uniform(-0.2,0.2))
		return (X,Z,Y)
	
	def updating(self,n):
		for i in range(1,n+1):
			self.n = i
			Data = self.construct_data()
			self.algorithm(Data[0],Data[1])
			if self.n%self.ST == 0:
				print('ab',self.Ln,self.n)
			T = 10000
			if self.n%self.ST == 0:
				S = Decimal(0)
				for j in range(T):
					Z = self.construct_data()
					S +=(self.f_nbar_cal(Z[0][0],Z[0][1],Z[0][2])-Z[2])**2
				self.Su.append(math.log10(S/T))
				print(S/T)
		
s = Decimal(10)/10
theta = 1/(2*s+1)
gamma = Decimal(2)/10
ST = 800
fstarorder = 10

ALG = algorithm(Decimal(1),theta,gamma,ST,fstarorder,0)
'''
for j in range(fstarorder):
	ALG.update_Z_Lnt()
ALG.Z_LntX_cal(1/math.sqrt(3),1/math.sqrt(3),1/math.sqrt(3))
fstar = ALG.ZnX
'''
fstar = 0

Algorithm = algorithm(s,theta,gamma,ST,fstarorder,fstar)
Algorithm.updating(100000)
S = Algorithm.Su
X = [ math.log10(i*ST+1) for i in range(len(S))]
s = float(s)
A = len(S)
Y = [ float(S[-1]) -(2*s/(2*s+1))*(i-(math.log10(A*ST))) for i in X]
plt.plot(X,S)
plt.plot(X[1:],Y[1:])
plt.show()

		
		
		
		
		
		
		
		
