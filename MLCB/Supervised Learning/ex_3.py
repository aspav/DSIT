import numpy as np
import matplotlib.pyplot as plt


def gaussian(x, mu, sigma):
    # return 1 / (np.sqrt(2 * np.pi) * sigma) * np.exp(-(1 / 2) * ((x - mu) / sigma) ** 2)
    return np.exp(-(x - mu) ** 2 / (2 * sigma)) / np.sqrt(2 * np.pi * sigma)


def draw_gaussian(mu, sigma):
    x = np.linspace(mu-3*np.sqrt(sigma), mu+3*np.sqrt(sigma), 1000)
    y = gaussian(x, mu, sigma)
    return x, y


# Calculate average x
def avg(x):
    avg_x = 1 / len(x) * sum(x)
    return avg_x


# Calculate the posterior distribution of mean and variance p(μ|Χ)
def posterior(N, x, mu_0, sigma_0, sigma):
    mu_N = (N * sigma_0 * avg(x) + sigma * mu_0) / (N * sigma_0 + sigma)
    sigma_N = (sigma * sigma_0) / (N * sigma_0 + sigma)
    return mu_N, sigma_N


mu_true = 7
sigma_true = 16
N = [1, 5, 10, 20, 50, 100, 1000]

mu_prior = 0
sigma_prior = 4

# x_train = np.random.normal(loc=mu_true, scale=np.sqrt(sigma_true), size=1)
#
# print(x_train)
with open('./output/ex3/posterior.txt', 'w') as file:
    file.write("Number of observations, posterior distribution's mean, posterior distribution's variance \n")
    for size in N:
        x_train = np.random.normal(loc=mu_true, scale=np.sqrt(sigma_true), size=size)
        mu_post, sigma_post = posterior(N=size, x=x_train, mu_0=mu_prior, sigma_0=sigma_prior, sigma=sigma_true)
        file.write(str(size) + ', ' + str(mu_post) + ', ' + str(sigma_post) + ', ' + '\n')

        # Figures
        plt.figure()
        plt.scatter(x_train, gaussian(x=x_train, mu=mu_true, sigma=sigma_true), facecolor="none", edgecolor="mediumseagreen", s=50, label="Training data")
        plt.plot(draw_gaussian(mu=mu_true, sigma=sigma_true)[0], draw_gaussian(mu=mu_true, sigma=sigma_true)[1], color='mediumseagreen', label='Distribution Generating the data')
        plt.plot(draw_gaussian(mu=mu_prior, sigma=sigma_prior)[0], draw_gaussian(mu=mu_prior, sigma=sigma_prior)[1], color='steelblue', label='Prior Distribution')
        plt.plot(draw_gaussian(mu=mu_post, sigma=sigma_post)[0], draw_gaussian(mu=mu_post, sigma=sigma_post)[1], color='firebrick', label='Estimated Posterior Distribution')
        plt.xlabel('$x$')
        plt.ylabel('$y$')
        plt.ylim(0)
        plt.title('$N=$' + str(size))
        plt.legend()
        plt.savefig('./output/ex3/' + str(size) + 'posterior.pdf')