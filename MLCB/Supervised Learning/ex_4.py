import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt


# Define the polynomial model
def polynomial_model(M, x):
    polynomial_terms = []
    for i in range(M+1):
        polynomial_terms.append(np.array(x)**i)
    return polynomial_terms


# Sine function
def sinusoidal(x):
    return np.sin(2 * np.pi * x)


# Create data  by adding gaussian noise
def generate_toy_data(sample_size, std):
    x = np.linspace(0, 1, sample_size)
    y = sinusoidal(x) + np.random.normal(scale=std, size=x.shape)
    return x, y


# Number of training data points
N = [10, 100]

# Degree of polynomial
M = [2 ,3 ,4 ,5, 9]

# Gaussian noise
m_noise = 0
s_noise = 1

# Test set Generation (random points over [0,1])
T = 1000

x_test = x = np.linspace(0, 1, T)
t = sinusoidal(x_test)


# Least squares method for different degree polynomials
# Vectorized Solution (coefficients = w = (X^T*X)^-1) * X^T * y) and predicted values y = X * w
with open('./output/ex4/Coefficients_RMSE.txt', 'w') as file:
    file.write("Sample Size \t degree of polynomial \t coefficients \t RMSE \n")
    for n in N:

        RMSE_list =[]

        x_train, y_train = generate_toy_data(sample_size=n, std=np.sqrt(s_noise))

        for m in M:
            X = np.array([polynomial_model(m, x) for x in x_train])
            w = (np.linalg.inv(X.T @ X) @ X.T) @ y_train
            print("For a "+ str(m) + " degree polynomial, the coefficients are: " + str(w))

            # Vectorized y_pred and y calculation
            y_pred = np.array([polynomial_model(m, x) for x in x_test]) @ w

            # Root Mean Squared Error (RMSE), RMSE = MSE^(1/2), MSE = 1/N * sum{(y_pred- t)^2}
            RMSE = np.sqrt((1 / m) * sum((y_pred - t) ** 2))
            RMSE_list.append(RMSE)
            print("Root Mean Squared Error: " + str(RMSE))

            file.write(str(n) + '\t' + str(m) + '\t' + str(w) + '\t' + str(RMSE) + '\n')

            plt.figure()
            plt.scatter(x_train, y_train, facecolor="none", edgecolor="mediumseagreen", s=50, label="Training data")
            plt.plot(x_test, t, color='steelblue', label='True Model')
            plt.plot(x_test, y_pred, color='firebrick', label='L-S Method' + ' (M=' + str(m) + ')')
            plt.xlabel('$x$')
            plt.xlim(-0.02, 1.01)
            plt.ylabel('$y$')
            if m == 2:
                # plt.title(str(m) + '$^{nd}$ Degree Polynomial')
                plt.title('Least Squares Method vs True Model \n ' + str(m) + '$^{nd}$ Degree Polynomial')
            elif m == 3:
                # plt.title(str(m) + '$^{rd}$ Degree Polynomial')
                plt.title('Least Squares Method vs True Model \n ' + str(m) + '$^{rd}$ Degree Polynomial')
            else:
                # plt.title(str(m) + '$^{th}$ Degree Polynomial')
                plt.title('Least Squares Method vs True Model \n ' + str(m) + '$^{th}$ Degree Polynomial')
            plt.legend()
            plt.savefig('./output/ex4/' + str(n) + '-' + str(m) + '-degree.pdf')

        plt.figure()
        plt.scatter(M, RMSE_list, facecolor="none", edgecolor="firebrick", s=30, label='values of RMSE')
        plt.plot(M, RMSE_list, color='mediumseagreen', label='RMSE')
        plt.title('Root Mean Square Error vs Degree of Polynomial \n $N=$' + str(n))
        plt.xlabel('Degree of Polynomial ($M$)')
        plt.xlim(1.9, 9.1)
        plt.ylabel('RMSE')
        plt.legend()
        plt.savefig('./output/ex4/' + str(n) + '-RMSE.pdf')
