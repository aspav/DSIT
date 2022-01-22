import numpy as np
from numpy.linalg import inv, norm
import matplotlib.pyplot as plt


# Define the polynomial model
def polynomial_model(M, x):
    polynomial_terms = []
    for i in range(M+1):
        polynomial_terms.append(np.array(x)**i)
    return polynomial_terms


# Define Sine function
def sinusoidal(x):
    return np.sin(2 * np.pi * x)


# Function for generating training data , where y_train = y_true + gaussian noise
def generate_toy_data(sample_size, std):
    x = np.linspace(0, 1, sample_size)
    y = sinusoidal(x) + np.random.normal(scale=std, size=x.shape)
    return x, y

#  Function that predicts m and s^2
def bayesian_prediction(a, b, Phi_test, Phi_train, y_train):
    Sigma = np.linalg.inv(a * np.eye(M + 1) + b * Phi_train.T @ Phi_train)
    m = b * Phi_test @ Sigma @ Phi_train.T @ y_train
    sigma = 1 / b + Phi_test @ Sigma @ Phi_test.T
    var = np.diag(sigma)

    return m, var


# Number of training data points
N = 10

# Degree of polynomial for φ(x)
M = 9

# Gaussian noise
m_noise = 0
beta = 11.1
s_noise = 1/beta

# Prior parameter
alpha = 0.0005


# Generate train data
x_train, y_train = generate_toy_data(sample_size=N, std=s_noise)

# Generate Φ(x)
Phi_train = []
for i in range(N):
    Phi_train.append(polynomial_model(M, x_train[i]))
Phi_train = np.array(Phi_train)

# Generate test data
x_test = np.linspace(0, 1, 1000)
Phi_test = np.array([polynomial_model(M=M, x=x) for x in x_test])

# Calculate the true and the predicted value
t = sinusoidal(x_test)
m, var = bayesian_prediction(a=alpha, b=beta, Phi_test=Phi_test, Phi_train=Phi_train, y_train=y_train)


plt.figure()
plt.scatter(x_train, y_train, facecolor="none", edgecolor="mediumseagreen", s=50, label="Training data")
plt.plot(x_test, t, color='steelblue', label='True Model')
plt.plot(x_test, m, color='firebrick', label='Bayesian Linear Regression')
plt.fill_between(x_test, m - var, m + var, alpha=0.2, color='firebrick', label='Variance')
plt.xlabel('$x$')
plt.xlim(-0.02, 1.01)
plt.ylabel('$y$')
plt.title('Bayesian Linear Regression vs True Model')
plt.legend()
plt.savefig('./output/ex5/bayesian_approach.pdf')



