#opt algorithms by lwj
import numpy as np
from numpy import linalg
import math
from math import sqrt,pow

def func(x):
	G = NULL 
	b = NULL
	c = NULL 

	return 1/2*x@G@x + b@x + c

def grad(x):
	return NULL

def Hess(x):
	return NULL 

def phi(alpha,x,d):
	return func(x + alpha*d)

def grad_phi(alpha,x,d):
	return NULL 

def Hess_phi(alpha,x,d):
	return NULL

# directions

def get_decent_direction(x):
	d = -1*grad(x)
	return d

def get_newton_direction(x):
	return -1*linalg.inv(Hess(x)) @ grad(x)

def get_conjugate_direction(x,x_old,d_old):

	def get_crowder_beta(g_new,g,d):
	return g_new @(g_new-g)/( d@(g_new - g))

	def get_fletcher_beta(g_new,g,d):
		return g_new@g_new / (g@g)

	def get_polak_beta(g_new,g,d):
		return g_new@(g_new - g)/(g@g)

	def get_dixon_beta(g_new,g,d):
		return -1*g_new@g_new/(d@g)

	g,g_old = grad(x),grad(x_old)
	beta = get_crowder_beta(g,g_old,d_old)

	return -1*g + beta*d_old

# alphas 
def get_interval(x,d):
	alpha0 = 0
	h0 = 0.01

	h = h0
	a = alpha0
	k = 0
	while k<1000:
		b = a + h
		if phi(b,x,d) < phi(a,x,d):
			a = b 
			h = 2*h
			k = k+1
		else:
			if k = 0:
				h = -1*h
				k = k+1

			else:
				return a,b

def get_0618alpha(x,d):
	a,b = get_interval(x,d)
	tol = 0.01

	while True:
		lam = a + 0.382*(b-a)
		mu = a + 0.618*(b-a)

		if phi(lam,x,d) < phi(mu,x,d):
			a = a
			b = mu
		else:
			a = lam
			b = b 


		if abs(a-b) < tol:
			return (a+b)/2


def get_bisection_alpha(x,d):
	a,b = get_interval(x,d)
	tol = 0.01

	while True:
		c = (a+b)/2

		if grad_phi(c,x,d) > 0:
			b = c
			a = a
		else:
			a = c
			b = b

		if abs(a-b) < tol:
			return (a+b)/2


def get_1pt2rd_alpha(x,d):
	a,b = get_interval(x,d)

	k = 0
	while abs(grad_phi(a,x,d)) > tol:
		a = a - grad_phi(a,x,d)/Hess_phi(a,x,d)
		k = k + 1
	return a

def get_2pt2rd_alpha(x,d):
	a,b = get_interval(x,d)

	k = 0
	while abs(grad_phi(a,x,d)) > tol:
		tmp = a
		a = a - 0.5*(a-b)*grad_phi(a,x,d)/(grad_phi(a,x,d)-(phi(a,x,d)-phi(b,x,d))/(a-b))
		b = tmp
		k = k + 1

	return a

def get_3pt2rd_alpha(x,d):
	a,b = get_interval(x,d)
	c = (a+b)/2

	k = 0
	while abs(grad_phi(a,x,d)) > tol:
		tmp1 = a
		tmp2 = b
		a = 0.5*(a+b) + 0.5*(phi(a,x,d)-phi(b,x,d))*(b-c)*(c-a)/( (b-c)*phi(a,x,d)+(c-a)*phi(b,x,d)+(a-b)*phi(c,d,x)  )
		b = tmp1
		c = tmp2 
		k = k + 1

	return a 

# check condition
def check_reduction(x_new,x,d):
	if func(x)-func(x_new) > 0:
		return True
	else:
		return False

def check_AG_condition(x_new,x,d):
	rho = 0.1

	if func(x)-func(x_new) > -1*rho*grad(x)@d and func(x)-func(x_new) < -1*(1-rho)*grad(x)@d:
		return True
	else:
		return False

def check_WP_condition(x_new,x,d):
	rho = 0.1
	sigma = 0.2

	if func(x)-func(x_new) > -1*rho*grad(x)@d and grad(x_new)@d > sigma*grad(x)@d:
		return True
	else:
		return False 

# algs
def get_newton_minimal():
	dim = 10

	tol = 1e-4

	k = 0
	x = np.random.randn(dim)
	while sqrt(grad(x)@grad(x)) > tol:
		d = get_newton_direction(x)
		alpha = get_bisection_alpha(x,d)
		x = x + alpha*d
		k = k + 1

	return x

def get_gradient_decent_minimal():
	dim = 10

	tol = 1e-4

	k = 0
	x = np.random.randn(dim)
	while  sqrt(grad(x)@grad(x)) > tol:
		d = get_decent_direction(x)
		alpha = get_bisection_alpha(x,d)
		x = x + alpha*d
		k = k + 1

	return x

def get_conjugate_minimal():
	dim = 10

	tol = 1e-4

	k = 0
	x_old = np.random.randn(dim)
	d_old = -1*grad(x)
	alpha = get_bisection_alpha(x,d)
	x = x + alpha*d

	while  sqrt(grad(x)@grad(x)) > tol:
		d = get_conjugate_direction(x,x_old,d_old)
		alpha = get_bisection_alpha(x,d)
		x = x + alpha*d
		k = k + 1

	return x



		



 
