% DETECTION & ESTIMATION THEORY / ASSIGNMENT-EXERCISE 3

% In this exercise, we 've got a data set which definetely derives from
% either an exponential or a Rayleigh distribution. We will execute two
% methods, in order to detect where the data come from, i.e to estimate the
% distribution (either exponential or Rayleigh) that the data set derives
% from.

% Load the data from 'data.mat' file.
DATA=load('data.mat');
data=DATA.data;

% First Step: Histogram Method.
% Given the data set, we will plot the histogram of the data as well as the
% Probability Density Functions of both exponential distribution (for an
% arbitrary parameter λ (lambda) -it doesn't matter, since we want to compare
% only the shape of the histogram with the shape of the pdf of exponential
% distribution) and Rayleigh distribution (for an arbitrary parameter σ
% (sigma) -it doesn't matter, since we want to compare only the shape of
% the histogram with the shape of the pdf of Rayleigh distribution).

% Histogram of data
figno=1;
figure(figno)
clf
histogram(data)
xlabel('Value')
ylabel('Frequency of Occurence')
title('Histogram of Data')
figno=figno+1;

% Plot the Probability Density Functions (PDFs) of both exponential and
% Rayleigh distribution.
% We seek for similarities in the shapes of histogram and some of these two
% distributions. As a result, it doeasn't matter what parameters λ and σ
% will be chosen for exponential and Rayleigh distribution in order to make
% the graphs of their pdfs.

% Define range of x values
x=0:0.1:10;
lambda_values=[0.5, 1, 1.5];
sigma_values=[0.5, 1, 2, 3, 4];

% PDF of exponential distribution, for three different values of the λ
% parameter of exponential distribution.
% In genereal: PDF of exponetial with parameter λ: p(x;λ)=λ*exp(-λx), for
% x>=0.

figure(figno)
hold on
for i=1:length(lambda_values)
    pdf_exp_values=lambda_values(i)*exp(-lambda_values(i)*x(2:length(x)));
    plot(x(2:length(x)),pdf_exp_values,'LineWidth',2,'DisplayName',['\lambda = ',num2str(lambda_values(i))])
end
xlabel('x')
ylabel('p(x) - Probability Density')
title('Exponential Distribution - Probability Density Function')
legend()
figno=figno+1;

% PDF of Rayleigh distribution, for five different values of the σ
% parameter of Rayleigh distribution.
% In genereal: PDF of Rayleigh with parameter σ: 
% p(x;σ)= [x / (σ^2)] * exp( [-(x^2) / (2*(σ^2))] ), for x>=0.
figure(figno)
hold on
for i=1:length(sigma_values)
    pdf_Rayleigh_values = ( x / (sigma_values(i).^2) ).*exp( -x.^2 / (2 * sigma_values(i).^2) );
    plot(x,pdf_Rayleigh_values,'LineWidth',2,'DisplayName',['\sigma = ',num2str(sigma_values(i))])
end
xlabel('x')
ylabel('p(x) - Probability Density')
title('Rayleigh Distribution - Probability Density Function')
legend()
figno=figno+1;

% Having the above three graphs (Histogram, pdf of exponential, pdf of
% Rayleigh), we can make an initial estimation of what distributiion the
% data derive from.

% Second Step: Empirical Cumulative Density Function Method, using MLE
% estimation Theory to estimate the parameters of both possible exponential
% and possible Rayleigh distribution that data may derive from.
% In this analysis, we use the collected data in order to result in
% Empirical Cumulative Density Function (ECDF). The ECDF is a function that
% represents the probability an observation in this data set to get smaller
% value than a specific value x in the data set, based on the collected
% data, i.e the data set itself. To calculate the ECDF given a data set,
% the data must first be sorted in ascending order. Seconldy, the ECDF is
% calculated by setting the ECDF value i/n to the i-th element of the
% vector that contains the sorted data, where n is the length of data set.
% If the data in the data set derived from some specific distibution, then
% the ECDF will have similarities in shape and maybe in values with the CDF
% of this distribution.
% So, since we know that the data come from either an exponential or a
% Rayleigh distribution, we will execute the Empirical Cumulative Density
% Function Method, as it will be described below.
% First of all we are going to compute the ECDF given the data. 
sorted_data=sort(data);
ecdf=(1:length(data)) / length(data);
figure(figno)
plot(sorted_data,ecdf,'b-','LineWidth',2)
xlabel('x (Data)')
ylabel('P(X<x) - Cumulative Probability')
title('Empirical Cumulative Density Function')
figno=figno+1;

% Then we will make two assumptions, one at a time.

% At first, we assume that the data really derive from an exponential
% distribution. The exponential distribution is defined by its parameter λ.
% Since, we don't know the actual value of the parameter λ of the possible
% exponential distribution that the data set may come from, we have to
% estimate it. We use MLE in order to estimate the unknown parameter λ that
% defines the possible exponential distribution. We use MLE, as it finds
% parameter estimates that are statistically optimal for the observed data.
% In order to execute the MLE estimator for the parameter λ, we assume that
% the data are iid. The assumption that the data are independant can be
% valid as we don't have any other evidence that shows a correlation
% between them. So, if we assume iid data, we can easily extract the joint
% pdf of the data just considering that the PDF of exponential with parameter
% λ is: p(x[i];λ)=λ*exp(-λx[i]), for x[i]>=0, where x[i] is the i-th element
% of data vector each time. The Estimation Theory gives the MLE for λ as
% the value of λ that maximizes the likelihood function p(x;λ) where
% x=[x[1], ..., x[N]] the data vector which is fixed. 
% So, MLE : λ_MLE = N / sum(x[i]) = 1/mean(x), where N is the length of
% data set, i.e data vector.
% We will use the value of MLE as the parameter λ in order to define the
% exponential distribution that is our estimation for the distribution that
% the data come from. To detect if that is a good estimation we will
% present the graphs of ecdf and cdf of the exponential distribution
% defined by parameter λ=λ_MLE together as well as a graph that presents
% the absolute error of ecdf from cdf above.

lambda_MLE = 1 / mean(data);
disp('-------- Assumption No1: Dataset comes from an exponential distribution --------')
disp('Estimate the parameter λ of the exponential distribution with MLE.')
disp(['The number of collected data in the dataset is: N=',num2str(length(data))])
disp(['The MLE for the parameter λ considering the collected data: λ_MLE=',num2str(lambda_MLE)])
fprintf('\n')

% Compute theoritical CDF of exponential distribution with parameter
% lambda_MLE.
theoritical_cdf_exp = 1 - exp(-lambda_MLE*sorted_data);
absolute_error_exp = abs(theoritical_cdf_exp-ecdf);

figure(figno)
tiledlayout(1,2)
nexttile
plot(sorted_data,ecdf,'b-','LineWidth',2,'DisplayName','Empirical CDF')
hold on
plot(sorted_data,theoritical_cdf_exp,'r-','LineWidth',2,'DisplayName',['Exponential CDF (\lambda =',num2str(lambda_MLE),')'])
xlabel('x (Data)')
ylabel('P(X<x) - Cumulative Probability')
title(['Empirical Cumulative Density Function vs Theoritical CDF of exp(\lambda =',num2str(lambda_MLE),')'])
legend('Location','Best');

nexttile
plot(sorted_data,absolute_error_exp,'r-','LineWidth',2)
xlabel('x')
ylabel('Absolute Error')
title(['Absolute Error of ECDF from Theoritical Exponential CDF with \lambda =',num2str(lambda_MLE)])

figno=figno+1;

% Seconldy, we assume that the data really derive from a Rayleigh
% distribution. The Rayleigh distribution is defined by its parameter σ.
% Since, we don't know the actual value of the parameter σ of the possible
% Rayleigh distribution that the data set may come from, we have to
% estimate it. We use MLE in order to estimate the unknown parameter σ that
% defines the possible Rayleigh distribution. We use MLE, as it finds
% parameter estimates that are statistically optimal for the observed data.
% In order to execute the MLE estimator for the parameter σ, we assume that
% the data are iid. The assumption that the data are independant can be
% valid as we don't have any other evidence that shows a correlation
% between them. So, if we assume iid data, we can easily extract the joint
% pdf of the data just considering that the PDF of Rayleigh with parameter
% σ is: p(x[i];σ)= [x[i] / (σ^2)] * exp( [-(x[i]^2) / (2*(σ^2))] ),
% for x[i]>=0 where x[i] is the i-th element of data vector each time. The
% Estimation Theory gives the MLE for σ as the value of σ that maximizes
% the likelihood function p(x;σ) where x=[x[1], ..., x[N]] the data vector
% which is fixed. 
% So, MLE : σ_MLE = sqrt( sum(x[i]^2) / N ), where N is the length of data
% set, i.e data vector.
% We will use the value of MLE as the parameter σ in order to define the
% Rayleigh distribution that is our estimation for the distribution that
% the data come from. To detect if that is a good estimation we will
% present the graphs of ecdf and cdf of the Rayleigh distribution
% defined by parameter σ=σ_MLE together as well as a graph that presents
% the absolute error of ecdf from cdf above.

squared_data=data.^2;
sigma_MLE=sqrt(mean(squared_data));

disp('-------- Assumption No2: Dataset comes from a Rayleigh distribution --------')
disp('Estimate the parameter σ of the Rayleigh distribution with MLE.')
disp(['The number of collected data in the dataset is: N=',num2str(length(data))])
disp(['The MLE for the parameter σ considering the collected data: σ_MLE=',num2str(sigma_MLE)])

% Compute theoritical CDF of exponential distribution with parameter
% lambda_MLE.
theoritical_cdf_Rayleigh = 1 - exp(-sorted_data.^2/(2*sigma_MLE.^2));
absolute_error_Rayleigh = abs(theoritical_cdf_Rayleigh-ecdf);

figure(figno)
tiledlayout(1,2)
nexttile
plot(sorted_data,ecdf,'b-','LineWidth',2,'DisplayName','Empirical CDF')
hold on
plot(sorted_data,theoritical_cdf_Rayleigh,'r-','LineWidth',2,'DisplayName',['Rayleigh CDF (\sigma =',num2str(sigma_MLE),')'])
xlabel('x (Data)')
ylabel('P(X<x) - Cumulative Probability')
title(['Empirical Cumulative Density Function vs Theoritical CDF of Rayleigh(\sigma =',num2str(sigma_MLE),')'])
legend('Location','Best');

nexttile
plot(sorted_data,absolute_error_Rayleigh,'r-','LineWidth',2)
xlabel('x')
ylabel('Absolute Error')
title(['Absolute Error of ECDF from Theoritical Exponential CDF with \sigma =',num2str(sigma_MLE)])

figno=figno+1;

% Comparing the results from two assumptions we can make our final
% estimation of what destribution, between exponential and Rayleigh, the 
% data derived from. We will have resulted in an estimation for the
% parameter that defines the distribution, i.e an estimaton for the exact
% exponential or Rayleigh distribution, as well.

% Moreover, plot again the histogram of the data as well as the pdfs of both
% the exponential distribution defined by the parameter λ which is estimated
% with λ_MLE and the Rayleigh distribution defined by the parameter σ which
% is estimated with σ_MLE. Mark what of the two pdf graphs has stronger
% resemblance with histogram.
figure(figno)
tiledlayout(2,2)

nexttile([1 2])
histogram(data)
xlabel('Value')
ylabel('Frequency of Occurence')
title('Histogram of Data')

x=min(data):0.01:max(data);
nexttile
pdf_exp_values=lambda_MLE*exp(-lambda_MLE*x);
plot(x,pdf_exp_values,'LineWidth',2)
xlabel('x')
ylabel('p(x) - Probability Density')
title(['Probability Density Function - Exponential Distribution with \lambda = ',num2str(lambda_MLE)])

nexttile
pdf_Rayleigh_values = ( x / (sigma_MLE.^2) ).*exp( -x.^2 / (2 * sigma_MLE.^2) );
plot(x,pdf_Rayleigh_values,'LineWidth',2)
xlabel('x')
ylabel('p(x) - Probability Density')
title(['Probability Density Function - Rayleigh Distribution with \sigma = ',num2str(sigma_MLE)])

figno=figno+1;