% DETECTION & ESTIMATION THEORY / ASSIGNMENT-EXERCISE 2

% Assume the measurement system x[n]=h*θ + w[n] , n=0,1,...,N-1 where w[n]
% is gaussian noise with zero mean and variance σ^2=4 and h is a known
% parameter of the measurement system. Let h=0.5 for all the following
% cases. The parameter θ is an unknown parameter to be estimated. We will
% consider different cases for this parameter θ and the measurement system
% and we will try to find out estimator for it in each case.

h=0.5;
% Declare the characteristics of the gaussian noise w of the measurement
% system(noise mean and noise variance).
mean_w=0;
var_w=4;

% 1) In this case, assume that the unknown parameter θ to be estimated is a
% gaussian random variable that follows normal distribution with mean μ=4
% and variance σ0^2=1. In fact, the parameter θ to be estimated is a
% realization of the random variable Θ with the above characteristics.
% In this analysis we will try to find the MMSE estimator for the θ using
% the Bayesian Theory and Approach.
% The measurement system is given by: x[n]=h*θ + w[n] , n=0,1,...,N-1 where
% N is the number of observation we take from the system. Let N=50
% observations. Considering the form of the measurement system, this can be
% presented by general bayesian linear model. Since both θ and w are
% gaussian and assuming that the noise w is independent of θ, we have the
% general linear gaussian model. So, the posterion pdf p(θ|x), where
% x=transpose( [x[0], x[1], ... x[N-1]] ), is also gaussian with mean
% E[θ|x]= mean(θ) + C(Θ)*transpose(H)*inverse(H*C(θ)*transpose(H)+C(w))*
% (x-H*mean(θ))
% or E[θ|x]=mean(θ) + inverse(transpose(H)*inverse(C(w))*H + inverse(C(θ)))*
% transose(H)*inverse(C(w))*(x-H*mean(θ)) where H is the measurement matrix
% and so, in this case, is a vector of length N with all its elements to be
% equal to h, C(θ) is the covariance matrix of θ in general and so, in this
% case, the variance of θ, C(w) is the covariance matrix of w.
% Note that Bmse is the covariance matrix of the posterior pdf.
N=50;

% Generate the noise for the measurement system.
w=mean_w + sqrt(var_w)*randn(N,1);
Cov_w=var_w*eye(N);

% Declare theta as the realization of the random variable θ that follows
% normal distribution with mean μ=4 and variance σ0^2=1. So, theta
% represents the unknown parameter to be estimated by MMSE.
mean_theta1=4;
var_theta1=1;
theta1=mean_theta1 + sqrt(var_theta1)*randn;

% Generate synthetic observations taken from the measurement system.
x=h*theta1 + w;

% Use general gaussian linear model analysis to declare the MMSE estimator
% for theta, since we know that MMSE(θ)=E[θ|x].
% Declare the observation matrix H.
H=h*ones(N,1);

% Compute MMSE and Bmse
MMSE_theta1 = mean_theta1 + var_theta1*transpose(H)*inv(H*var_theta1*transpose(H) + Cov_w)*(x-H*mean_theta1);
Bmse_1 = inv(1/var_theta1 + transpose(H)*inv(Cov_w)*H);

% Display results.
disp("-------               MMSE for θ / Results               -------")
disp("------- Measurement System: x[n]=h*θ+w[n], n=0,1,...,N-1 -------")
fprintf('\n')
disp(['Known parameter h=',num2str(h) ,' . Number of observations taken: N=',num2str(N)])
disp(['Take θ from distribution N(4,1). The actual value of the parameter θ: θ=', num2str(theta1)])
disp(['The MMSE estimator using Bayessian Approach for General Gaussian Linear Model: θ_MMSE=',num2str(MMSE_theta1)])
disp(['The Bmse: Bmse=',num2str(Bmse_1)])

% Do the same for 5 different experiments in order to see how the results
% may change.
figno=1;
figure(figno)
hold on
for i=1:5
    theta_true= mean_theta1 +sqrt(var_theta1)*randn;
    plot(i,theta_true,'.b','MarkerSize',15,'DisplayName','True value of θ')
    w=sqrt(var_w)*randn(N,1);
    x = h*theta_true + w;
    MMSE_theta = mean_theta1 + var_theta1*transpose(H)*inv(H*var_theta1*transpose(H) + Cov_w)*(x-H*mean_theta1);
    Bmse = inv(1/var_theta1 + transpose(H)*inv(Cov_w)*H);
    
    plot(i,MMSE_theta,'.r','MarkerSize',15,'DisplayName','MMSE of θ')

    % Display Results
    disp(['---- Experiment number ',num2str(i),' ----'])
    disp(['Actual value of the parameter θ: θ=', num2str(theta_true)])
    disp(['The MMSE estimator : θ_MMSE=',num2str(MMSE_theta)])
    disp(['The Bmse: Bmse=',num2str(Bmse)])
end
legend('True value of θ','MMSE of θ','Location','Best')
title('MMSE of parameter θ (when derived from gaussian distribution)')
grid on;
figno=figno+1;

%-------------------------------------------------------------------------
% 2) In this case, θ is considered as the unknown parameter we want to estimate
% but it is now deterministic. We have the same model for the measurement
% system, i.e the general linear model. We want to find out the MLE for the
% unknown parameter θ. The MLE is the value of θ that maximizes the
% likelihood function p(x|θ) for fixed x. Since we have general linear
% model, we can use the expression derived from MLE Approach for GLM:
% θ_MLE=inv(transpose(H)*inv(C(w))*H)*transpose(H)*inv(C(w))*x, where C(w)
% and H as before. 
% The covariance matrix for the MLE (since MLE is a random variable now that
% follows normal distribution) is given by: C(θ_MLE)=inv(transpose(H)*C(w)*H).
% In this case, where θ is scalar, the above expression gives the variance
% of the MLE.
% Again we consider the number of observations taken from the system to be
% equal to 50.
N=50;
H=h*ones(N,1);

% Generate the gaussian noise vector.
w=mean_w + sqrt(var_w)*randn(N,1);
Cov_w=var_w*eye(N);

% θ is now deterministic but considered as unknown and we want to estimate
% it. We choose an arbitrary value as the value of θ not taken from any
% distribution. Assume θ=4, in order to be equal to μ.
theta2=4;

% Generate synthetic observations taken from the measurement system.
x = h*theta2 + w;

% Find MLE and the variance of the MLE.
var_MLE_theta2=inv(transpose(H)*inv(Cov_w)*H);
MLE_theta2 = var_MLE_theta2*transpose(H)*inv(Cov_w)*x;

% Display results.
fprintf('\n')
disp("-------               MLE for θ / Results               -------")
disp("------- Measurement System: x[n]=h*θ+w[n], n=0,1,...,N-1 -------")
fprintf('\n')
disp(['Known parameter h=',num2str(h) ,' . Number of observations taken: N=',num2str(N)])
disp(['The unknown parameter θ is deterministic. So the value of θ is setted arbitrarily: θ=', num2str(theta2)])
disp(['The MLE estimator using the MLE Approach for General Linear Model: θ_MLE=',num2str(MLE_theta2)])
disp(['The variance of the MLE for estimation of θ: Var(θ_MLE)=',num2str(var_MLE_theta2)])

% Do the same for 5 different experiments in order to see how the results
% may change.
figure(figno)
hold on
for i=1:5
    theta_true=4; % Arbitrary choice of θ. We want it to be deterministic but unknown.
    plot(i,theta_true,'.b','MarkerSize',15,'DisplayName','True value of θ')
    w=sqrt(var_w)*randn(N,1);
    x = h*theta_true + w;
    var_MLE_theta=inv(transpose(H)*inv(Cov_w)*H);
    MLE_theta = var_MLE_theta*transpose(H)*inv(Cov_w)*x;

    plot(i,MLE_theta,'.r','MarkerSize',15,'DisplayName','MLE of θ')

    % Display Results
    disp(['---- Experiment number ',num2str(i),' ----'])
    disp(['Actual value of the parameter θ: θ=', num2str(theta_true)])
    disp(['The MLE estimator : θ_MLE=',num2str(MLE_theta)])
    disp(['The variance of the MLE: Var(θ_MLE)=',num2str(var_MLE_theta)])
end
legend('True value of θ','MLE of θ','Location','Best')
title('MLE of deteministic (but unknown) parameter θ')
grid on;
figno=figno+1;

%-------------------------------------------------------------------------
% 3) In the third step, we want to examine for both the two estimators
% (MMSE and MLE), we found above, the convergence of the Mean Squared Error
% (MSE) as a function of the number of observations taken per experiment.
% So, we will make the graph of the MSE as a function of the number of
% observations taken per experiment for 500 measurement experiments and a
% maximum number of observations per experiment 50, executing the whole
% process for both MMSE and MLE estimators for the unknown parameter θ of
% the measurement system.
% In order to examine the convergence of the MSE of some of the two
% estimators, we will perform 500 measurment experiments varying the number
% of observations taken per experiment from 1 to 50 while executing them.
% An experiment is actually the process of collecting a number of data x[n]
% from the measurement system x[n]=h*θ + w[n] , n=0,1,...,N-1 where N is
% exactly this number of observations taken per experiment that is going to
% vary. 
% We want to make the graph of MSE for each estimator as a function of the
% number of observations taken per experiment varying the number of
% observations taken per experiment from 1 to 50. So, we need to distribute
% the 500 experiments across these different values of the number of
% observations taken per experiment, i.e the N. Since we have a total
% number of experiments 500 and a maximum number of observations taken per
% experiment 50, we can evenly distribute these experiments so that for
% each value of N from 1 to 50 we perform 500/50=10 experiments and take
% the Mean Squared Error of the estimation of θ that refers to a specific
% value of the number of obsrvations taken per experiment (N) each time.
% The Mean Squared Error of the estimation of the unknown parameter θ, that
% refers to the estimation of θ while taking N observations per experiment,
% is given by averaging over all the 'Squared Errors' of the estimation of
% θ that derive from the 10 different experiments which refer to this
% specific number N of observations taken in each of them.

% a) In this case, θ is considered to be gaussian with mean μ=4 and
% variance σ0^2=1. The measurement system is a general linear gaussian
% model: x[n]=h*θ + w[n] , n=0,1,...,N-1 , so the MMSE for the unknown θ is
% given by (as we declare in the first step of the exercise):
% E[θ|x]= mean(θ) + C(Θ)*transpose(H)*inverse(H*C(θ)*transpose(H)+C(w))*
% (x-H*mean(θ))
% or E[θ|x]=mean(θ) + inverse(transpose(H)*inverse(C(w))*H + inverse(C(θ)))*
% transose(H)*inverse(C(w))*(x-H*mean(θ)) where H is the measurement matrix
% and so, in this case, is a vector of length N with all its elements to be
% equal to h, C(θ) is the covariance matrix of θ in general and so, in this
% case, the variance of θ, C(w) is the covariance matrix of w.

% Declare the mean and the variance of the θ.
mean_theta3a=4;
var_theta3a=1;

N_total_experiments=500; % Total number of experiments to be made.
n_observations_max=50; % Maximum number of observations per experiment.
 
% Number of experiments to be made for each value of number of observations
% taken per experiment.
N_experiments=N_total_experiments/n_observations_max;

% Define MSE as a vector of length n_observations_max.
MSE=zeros(n_observations_max,1);

for n_obs=1:n_observations_max
    MSE_exp_per_n_obs=zeros(N_experiments,1);

    % Define the observation matrix H.
    H=h*ones(n_obs,1);

    for n_exp=1:N_experiments
        theta3a=mean_theta3a + sqrt(var_theta3a)*randn; % Reliazation of θ.
        
        % Generate the gaussian noise vector.
        w=mean_w + sqrt(var_w)*randn(n_obs,1);
        Cov_w=var_w*eye(n_obs);
        
        x = h*theta3a + w; % Generate synthetic data

        % Compute MMSE and the Squared Error for the estimation of the θ.
        MMSE_theta3a = mean_theta3a + var_theta3a*transpose(H)*inv(H*var_theta3a*transpose(H) + Cov_w)*(x-H*mean_theta3a);
        MSE_exp_per_n_obs(n_exp) = (theta3a-MMSE_theta3a).^2;
    end

    % Average over the 'Squared Errors' computed from the experiments for
    % this value of number of observations to take the MSE for this number.
    MSE(n_obs) = mean(MSE_exp_per_n_obs);

end

figure(figno)
clf
plot(1:n_observations_max, MSE, '.b','MarkerSize',15)
ylabel('Mean Squared Error (MSE)');
xlabel('Number of observations');
title('Convergence of MSE for MMSE Estimator of unknown parameter θ (θ taken from Normal Distribution (4,1) )');
grid on;
figno=figno+1;


% b) In this case, θ is considered to be deterministic but unknown. Since,
% we have general linear model: x[n]=h*θ + w[n] , n=0,1,...,N-1 , the MLE
% is given by: θ_MLE=inv(transpose(H)*inv(C(w))*H)*transpose(H)*inv(C(w))*x.

% θ is now deterministic but considered as unknown and we want to estimate
% it. We choose an arbitrary value as the value of θ not taken from any
% distribution. Assume θ=4.
theta3b=4;

% Generate synthetic observations taken from the measurement system.
x = h*theta3b + w;

N_total_experiments=500; % Total number of experiments to be made.
n_observations_max=50; % Maximum number of observations per experiment.
 
% Number of experiments to be made for each value of number of observations
% taken per experiment.
N_experiments=N_total_experiments/n_observations_max;

% Define MSE as a vector of length n_observations_max.
MSE=zeros(n_observations_max,1);

for n_obs=1:n_observations_max
    MSE_exp_per_n_obs=zeros(N_experiments,1);

    % Define the observation matrix H.
    H=h*ones(n_obs,1);

    for n_exp=1:N_experiments
        % θ is deterministic and is defined with the value 2. However, it
        % is considered as unknown parameter to be estimated.
        
        % Generate the gaussian noise vector.
        w=mean_w + sqrt(var_w)*randn(n_obs,1);
        Cov_w=var_w*eye(n_obs);
        
        x = h*theta3b + w; % Generate synthetic data

        % Compute MLE and the Squared Error for the estimation of the θ.
        var_MLE_theta3b=inv(transpose(H)*inv(Cov_w)*H);
        MLE_theta3b = var_MLE_theta3b*transpose(H)*inv(Cov_w)*x;
        MSE_exp_per_n_obs(n_exp) = (theta3b-MLE_theta3b).^2;
    end

    % Average over the 'Squared Errors' computed from the experiments for
    % this value of number of observations to take the MSE for this number.
    MSE(n_obs) = mean(MSE_exp_per_n_obs);

end

figure(figno)
clf
plot(1:n_observations_max, MSE, '.b','MarkerSize',15)
ylabel('Mean Squared Error (MSE)');
xlabel('Number of observations');
title('Convergence of MSE for MLE Estimator of unknown parameter θ (θ is determinitic but unknown parameter)');
grid on;
figno=figno+1;

%-------------------------------------------------------------------------
% 4) In this step, we will execute the Bayesian Approach to find the MMSE
% estimator for the unknown parameter θ considering the observations x[n]
% taken from the measurement system x[n]=h*θ + w[n] , n=0,1,...,N-1 , but
% the choice of prior of the unknown parameter θ is now "wrong". We will
% examine the consequences of a "wrong" choice of the prior distribution of
% θ to the estimation of the parameter usisng the MMSE estimator.
% To model this "wrong" choice of prior pdf for the unknown parameter θ to
% be estimated, we will keep the same model for the MMSE estimator we
% found in the first step of the exercise and we will derive the value of
% the unknown parameter θ from a) uniform distribution and b) exponential
% distribution. The model we found for the MMSE estimator of θ refers to
% the MMSE when the prior pdf of θ is gaussian and so we have General
% Linear Gaussian Model for the measurement system. This MMSE model is the
% same with the model for LMMSE when we have General Linear Model. However,
% the LMMSE is an optimal Bayesian estimator given in closed form but it's
% not exactly the MMSE for all the general cases. In fact, LMMSE is the
% MMSE only when we have jointly gaussian assumtpions for the unknown
% parameter θ and the data x. So, since we will assume a) uniform and b)
% exponential distribution as the prior of the unkown parameter θ, if we
% use the model of MMSE we found before, we will result in the model of
% LMMSE for the unknown parameter θ considering the data.
% In order to investigate the consequences of the "wrong" choices for the
% prior pdf of θ, we are going to make the graph of MSE for each MMSE 
% estimator (i.e in each case a) or b) ) as a function of the number of
% observations taken per experiment varying the number of observations
% taken per experiment from 1 to 50, as we did in the 3rd step.

% a) Assume that θ is a random variable that follow a uniform distribution.
% Let θ is uniform distributed in [A,B]=[0,8], so that θ follows a uniform
% distribution with mean equal to μ. Declare the mean and the variance of θ.
% Reminder: The unknown parametet θ we want to estimate from the data taken
% from the measurement system is a realization of this random variable.
a=0;
b=8;
mean_theta4a=(a+b)/2;
var_theta4a=((b-a).^2) / 12;
% theta4a=unifrnd(a,b)

N_total_experiments=500; % Total number of experiments to be made.
n_observations_max=50; % Maximum number of observations per experiment.
 
% Number of experiments to be made for each value of number of observations
% taken per experiment.
N_experiments=N_total_experiments/n_observations_max;

% Define MSE as a vector of length n_observations_max.
MSE=zeros(n_observations_max,1);

for n_obs=1:n_observations_max
    MSE_exp_per_n_obs=zeros(N_experiments,1);

    % Define the observation matrix H.
    H=h*ones(n_obs,1);

    for n_exp=1:N_experiments
        theta4a=unifrnd(a,b); % Reliazation of θ.
        
        % Generate the gaussian noise vector.
        w=mean_w + sqrt(var_w)*randn(n_obs,1);
        Cov_w=var_w*eye(n_obs);
        
        x = h*theta4a + w; % Generate synthetic data

        % Compute MMSE and the Squared Error for the estimation of the θ.
        MMSE_theta4a = mean_theta4a + var_theta4a*transpose(H)*inv(H*var_theta4a*transpose(H) + Cov_w)*(x-H*mean_theta4a);
        MSE_exp_per_n_obs(n_exp) = (theta4a-MMSE_theta4a).^2;
    end

    % Average over the 'Squared Errors' computed from the experiments for
    % this value of number of observations to take the MSE for this number.
    MSE(n_obs) = mean(MSE_exp_per_n_obs);

end

figure(figno)
clf
plot(1:n_observations_max, MSE, '.b','MarkerSize',15)
ylabel('Mean Squared Error (MSE)');
xlabel('Number of observations');
title(['Convergence of MSE for MMSE Estimator of unknown parameter θ (θ taken from Uniform Distribution U[',num2str(a),',',num2str(b),'] )']);
grid on;
figno=figno+1;


% b) Assume that θ is a random variable that follow an expotential distribution.
% Let θ follows exponential distribution with parameter λ=1/4, so that θ
% follows an exponential distribution with mean equal to μ. Declare the
% mean and the variance of the θ.
% Reminder: The unknown parametet θ we want to estimate from the data taken
% from the measurement system is a realization of this random variable.
lamba=1/4;
mean_theta4b=1/lamba;
var_theta4b=1/(lamba.^2);
% theta4b=exprnd(mean_theta4b)

N_total_experiments=500; % Total number of experiments to be made.
n_observations_max=50; % Maximum number of observations per experiment.
 
% Number of experiments to be made for each value of number of observations
% taken per experiment.
N_experiments=N_total_experiments/n_observations_max;

% Define MSE as a vector of length n_observations_max.
MSE=zeros(n_observations_max,1);

for n_obs=1:n_observations_max
    MSE_exp_per_n_obs=zeros(N_experiments,1);

    % Define the observation matrix H.
    H=h*ones(n_obs,1);

    for n_exp=1:N_experiments
        theta4b=exprnd(mean_theta4b); % Reliazation of θ.
        
        % Generate the gaussian noise vector.
        w=mean_w + sqrt(var_w)*randn(n_obs,1);
        Cov_w=var_w*eye(n_obs);
        
        x = h*theta4b + w; % Generate synthetic data

        % Compute MMSE and the Squared Error for the estimation of the θ.
        MMSE_theta4b = mean_theta4b + var_theta4b*transpose(H)*inv(H*var_theta4b*transpose(H) + Cov_w)*(x-H*mean_theta4b);
        MSE_exp_per_n_obs(n_exp) = (theta4b-MMSE_theta4b).^2;
    end

    % Average over the 'Squared Errors' computed from the experiments for
    % this value of number of observations to take the MSE for this number.
    MSE(n_obs) = mean(MSE_exp_per_n_obs);

end

figure(figno)
clf
plot(1:n_observations_max, MSE, '.b','MarkerSize',15)
ylabel('Mean Squared Error (MSE)');
xlabel('Number of observations');
title(['Convergence of MSE for MMSE Estimator of unknown parameter θ (θ taken from Exponential Distribution with λ=',num2str(lamba),')']);
grid on;
figno=figno+1;